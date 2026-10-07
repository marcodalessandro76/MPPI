"""
Tests of the Optics module. Most of them compare with the analytical susceptibilities of the classical anharmonic oscillator
(Boyd, Nonlinear Optics, Section 1.4), computed with the AnharmonicOscillator model.
"""
import os, io, contextlib
import numpy as np
import pytest
from mppi.Optics.Utils import fit_sum_frequencies, eval_sum_frequencies
from mppi import Optics as O
from mppi.Optics.Xn_frequency_mixing import field_amplitude
from mppi.Models.AnharmonicOscillator import AnharmonicOscillator
from mppi.Utilities.Constants import HaToeV

T0 = 5.          # switch on time of the fields (au), a nonzero value tests the phases
WP = 0.05        # pump frequency (Ha)
FREQS = np.array([0.32, 0.47, 0.61])
TIME = np.arange(0, 1800, 0.25)

def rel(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b))

def quiet():
    return contextlib.redirect_stdout(io.StringIO())

def test_fit_sum_frequencies():
    t = np.linspace(0,50,2000)
    Omegas = {1:1.0,2:2.0}
    y = 0.3 + 2.0*np.sin(1.0*t+0.4) + 0.5*np.sin(2.0*t-1.0)
    res, B0, residual = fit_sum_frequencies(t,y,Omegas)
    assert np.isclose(B0,0.3)
    assert np.isclose(res[1]['A'],2.0) and np.isclose(res[1]['phi'],0.4)
    assert np.isclose(res[2]['A'],0.5) and np.isclose(res[2]['phi'],-1.0)
    assert residual < 1e-8
    assert np.allclose(eval_sum_frequencies(t,res,B0),y)

def test_fit_sum_frequencies_underdetermined():
    # less samples than unknowns: lstsq returns empty residuals
    t = np.array([0.,1.,2.])
    res, B0, residual = fit_sum_frequencies(t,np.sin(t),{1:1.0,2:2.0})
    assert residual >= 0.

def test_xn_from_file_rejects_delta_field(ref_dir):
    # the Xn classes need sine-shaped fields, the LiF reference run uses a delta-shaped one
    file = os.path.join(ref_dir,'nl_results','LiF-delta_pulse','ndb.Nonlinear')
    for cls in [O.Xn_single_frequency, O.Xn_frequency_mixing]:
        with pytest.raises(ValueError):
            cls.from_file(file,verbose=False)

@pytest.fixture(scope='module')
def osc2():
    return AnharmonicOscillator(omega0=0.5, gamma=0.02, a=0.1)

@pytest.fixture(scope='module')
def osc3():
    return AnharmonicOscillator(omega0=0.5, gamma=0.02, b=0.1)

def test_field_amplitude():
    # E0 sin(w(t-t0)) = E(w) exp(-iwt) + E(-w) exp(iwt) with E(-w) = E(w)^*
    E0, w, t0, t = 2., 0.3, 1.5, np.linspace(0, 20, 7)
    Ep, Em = field_amplitude(E0, w, t0, 1), field_amplitude(E0, w, t0, -1)
    assert np.isclose(Em, np.conj(Ep))
    assert np.allclose(Ep*np.exp(-1j*w*t) + Em*np.exp(1j*w*t), E0*np.sin(w*(t-t0)))
    assert np.isclose(field_amplitude(E0, w, t0, -2), np.conj(Ep)**2)
    assert field_amplitude(E0, w, t0, 0) == 1.

def test_linear_response():
    osc = AnharmonicOscillator(omega0=0.5, gamma=0.01)
    time = np.arange(0, 3000, 0.5)
    data = osc.delta_response(time, amplitude=1e-3, initial_time=10.)
    with quiet():
        energy, eps = O.Linear_Response(time, data.Polarization[0], data.Efield[0], eta=0.)
    w = energy/HaToeV
    sel = (w > 0.2) & (w < 0.8)
    assert rel(eps[sel]-1, 4*np.pi*osc.chi1(w[sel])) < 2e-3

def test_single_frequency(osc2):
    data = osc2.sine_response(TIME, FREQS, amplitude=1e-5, initial_time=T0, substeps=4)
    with quiet():
        chi = O.Xn_single_frequency(data, X_order=2, verbose=False).compute_Xn()[0]
    assert rel(chi[1], osc2.chi1(FREQS)) < 1e-4
    assert rel(chi[2], osc2.chi2(FREQS, FREQS)) < 1e-4
    # optical rectification, degeneracy factor D=2
    assert rel(chi[0], 2*osc2.chi2(FREQS, -FREQS)) < 1e-4

def test_frequency_mixing_second_order(osc2):
    pump = {'frequency': WP, 'amplitude': 1e-4, 'initial_time': T0}
    data = osc2.sine_response(TIME, FREQS, amplitude=1e-5, initial_time=T0, pump=pump, substeps=4)
    with quiet():
        chi = O.Xn_frequency_mixing(data, X_order=(1, 1), verbose=False).compute_Xn()[0]
    assert rel(chi[(1, 0)], osc2.chi1(FREQS)) < 1e-4
    assert rel(chi[(0, 1)], osc2.chi1(WP)) < 1e-4
    # sum and difference frequency generation, degeneracy factor D=2
    assert rel(chi[(1, 1)], 2*osc2.chi2(FREQS, WP)) < 5e-3
    assert rel(chi[(1, -1)], 2*osc2.chi2(FREQS, -WP)) < 5e-3

def test_frequency_mixing_third_order(osc3):
    EP = 2e-3
    pump = {'frequency': WP, 'amplitude': EP, 'initial_time': T0}
    data = osc3.sine_response(TIME, FREQS, amplitude=1e-5, initial_time=T0, pump=pump, substeps=4)
    probe = osc3.sine_response(TIME, FREQS, amplitude=1e-5, initial_time=T0, substeps=4)
    with quiet():
        chi = O.Xn_frequency_mixing(data, X_order=(1, 3), verbose=False).compute_Xn()[0]
        chi1 = O.Xn_single_frequency(probe, X_order=1, verbose=False).compute_Xn()[0][1]
    # degeneracy factor D=3
    assert rel(chi[(1, 2)], 3*osc3.chi3(FREQS, WP, WP)) < 1e-2
    assert rel(chi[(1, -2)], 3*osc3.chi3(FREQS, -WP, -WP)) < 1e-2
    # chi3(w;w,wP,-wP) from the pump induced change of the linear term, D=6
    x3 = (chi[(1, 0)]-chi1)/(EP**2/4)
    assert rel(x3, 6*osc3.chi3(FREQS, WP, -WP)) < 1e-2

def test_frequency_mixing_sampling_after_dephasing():
    # at wp=0.32 the closest harmonics are wp-3wP and 3wP, the estimated time window is shorter than the simulation
    # and used to start before the dephasing time, where the free oscillations at omega0 spoil the fit (error ~7 on (1,2))
    osc = AnharmonicOscillator(omega0=0.5, gamma=0.01, b=0.1)
    freqs = np.array([0.32, 0.30+7*0.45/135])
    EP = 2e-3
    pump = {'frequency': WP, 'amplitude': EP, 'initial_time': T0}
    data = osc.sine_response(np.arange(0, 3000, 0.25), freqs, amplitude=1e-5, initial_time=T0, pump=pump, substeps=4)
    with quiet():
        xn = O.Xn_frequency_mixing(data, X_order=(1, 3), verbose=False)
        chi = xn.compute_Xn()[0]
    for ifreq in range(len(freqs)):
        assert xn._time_sampling(ifreq)[0]*xn.dt >= xn.deph - xn.dt
    assert rel(chi[(1, 2)], 3*osc.chi3(freqs, WP, WP)) < 2e-2

def test_harmonic_fits_are_cached(osc3, monkeypatch):
    import sys
    XFM, XSF = sys.modules['mppi.Optics.Xn_frequency_mixing'], sys.modules['mppi.Optics.Xn_single_frequency']
    calls = []
    def counting_fit(*args, **kwargs):
        calls.append(1)
        return fit_sum_frequencies(*args, **kwargs)
    monkeypatch.setattr(XFM, 'fit_sum_frequencies', counting_fit)
    monkeypatch.setattr(XSF, 'fit_sum_frequencies', counting_fit)
    pump = {'frequency': WP, 'amplitude': 2e-3, 'initial_time': T0}
    data = osc3.sine_response(TIME, FREQS, amplitude=1e-5, initial_time=T0, pump=pump, substeps=1)
    probe = osc3.sine_response(TIME, FREQS, amplitude=1e-5, initial_time=T0, substeps=1)
    with quiet():
        for xn in [O.Xn_frequency_mixing(data, X_order=(1, 3), verbose=False),
                   O.Xn_single_frequency(probe, X_order=3, verbose=False)]:
            calls.clear()
            chi = xn.compute_Xn()[0]
            xn.eval_Pw(), xn.eval_Ew(), xn.check_harmonic_reliability()
            assert len(calls) == 3*len(FREQS)   # one fit per frequency and direction
            # the stored fits are not used if the parameters of the analysis change
            xn.Trange_units = 'au'
            xn.Trange = [-1, TIME[-1]/2]
            xn.compute_Xn()
            assert len(calls) == 6*len(FREQS)

def test_sampling_warnings_summary(osc3):
    # wp close to 6wP: the optimal time window is longer than the simulation, the warning is printed once
    freqs = np.array([0.3033, 0.40, 0.3067])
    pump = {'frequency': WP, 'amplitude': 2e-3, 'initial_time': T0}
    data = osc3.sine_response(TIME, freqs, amplitude=1e-5, initial_time=T0, pump=pump, substeps=1)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        O.Xn_frequency_mixing(data, X_order=(1, 3), verbose=False).compute_Xn()
    lines = [l for l in out.getvalue().splitlines() if 'Warning' in l]
    assert len(lines) == 1
    assert 'optimal time window' in lines[0] and '(2 of 3 frequencies, indexes: 0, 2)' in lines[0]

def test_single_frequency_no_spurious_warning():
    # time grid and one of the frequencies of a yambo_nl run of the LiF project (200 fs, IO step 0.01 fs, 201 frequencies
    # in 10-20 eV): in floating point (t[-1]-T)+T > t[-1], which raised a spurious warning on the length of the time sampling
    from mppi.Utilities.Constants import FsToAu
    time = np.arange(20001)*0.4134137333656136
    freq = np.array([0.5667806500473366])
    T, tend = 2*np.pi/freq[0]/FsToAu, time[-1]/FsToAu
    assert (tend-T)+T > tend
    osc = AnharmonicOscillator(omega0=0.5, gamma=0.02)
    data = osc.sine_response(time, freq, amplitude=1e-5, initial_time=T0, substeps=2)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        chi = O.Xn_single_frequency(data, X_order=1, verbose=False).compute_Xn()[0]
    assert 'Warning' not in out.getvalue()
    assert rel(chi[1], osc.chi1(freq)) < 1e-4
