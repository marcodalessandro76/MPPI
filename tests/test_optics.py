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
    # time grid of a typical yambo_nl run (200 fs, IO step 0.01 fs) and a frequency of a grid of 201 points in 10-20 eV:
    # in floating point (t[-1]-T)+T > t[-1], which raised a spurious warning on the length of the time sampling
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

def test_lorentzian_broadening():
    from mppi.Optics.Utils import lorentzian_broadening
    osc, delta = AnharmonicOscillator(omega0=0.5, gamma=0.01), 0.01
    # on a wide grid the convolution gives chi(w + i delta) far from the edges
    w = np.linspace(0., 2., 2001)
    sel = (w > 0.3) & (w < 0.75)
    assert rel(lorentzian_broadening(w, osc.chi1(w), delta, renormalize=False)[sel], osc.chi1(w[sel]+1j*delta)) < 1e-3
    # on a narrow window the error comes from the missing tails, reduced by the normalization of the kernel
    w = np.linspace(0.30, 0.75, 136)
    err = [rel(lorentzian_broadening(w, osc.chi1(w), delta, renormalize=r), osc.chi1(w+1j*delta)) for r in [True, False]]
    assert err[0] < 4e-2 and err[1] < 7e-2
    # the last axis runs over the frequencies
    chi = np.array([osc.chi1(w), 2*osc.chi1(w)])
    assert np.allclose(lorentzian_broadening(w, chi, delta)[1], 2*lorentzian_broadening(w, osc.chi1(w), delta))
    with pytest.raises(ValueError):
        lorentzian_broadening(w[::-1], osc.chi1(w), delta)

def test_compute_Xn_broadening(osc3):
    # keys linear in the probe, compared with the analytic susceptibilities at the complex probe frequency
    w, delta, EP = np.linspace(0.30, 0.75, 61), 0.01, 2e-3
    pump = {'frequency': WP, 'amplitude': EP, 'initial_time': T0}
    data = osc3.sine_response(TIME, w, amplitude=1e-5, initial_time=T0, pump=pump, substeps=2)
    probe = osc3.sine_response(TIME, w, amplitude=1e-5, initial_time=T0, substeps=2)
    with quiet():
        mix = O.Xn_frequency_mixing(data, X_order=(1, 3), verbose=False)
        single = O.Xn_single_frequency(probe, X_order=3, verbose=False)
        chi, chi_b = mix.compute_Xn()[0], mix.compute_Xn(broadening=delta*HaToeV)[0]
        sf, sf_b = single.compute_Xn()[0], single.compute_Xn(broadening=delta*HaToeV)[0]
    # the Lorentzian tails outside the window give an error ~(delta/pi)(1/d1+1/d2) ~ 3% on the whole window
    z = w + 1j*delta
    assert rel(sf_b[1], osc3.chi1(z)) < 5e-2
    assert rel(chi_b[(1, 2)], 3*osc3.chi3(z, WP, WP)) < 5e-2
    assert rel(chi_b[(1, -2)], 3*osc3.chi3(z, -WP, -WP)) < 5e-2
    assert rel((chi_b[(1, 0)]-sf_b[1])/(EP**2/4), 6*osc3.chi3(z, WP, -WP)) < 5e-2
    # the other keys are not broadened
    for key in [0, 2, 3]:
        assert np.array_equal(sf_b[key], sf[key])
    for key in [(0, 0), (0, 1), (0, 2), (0, 3)]:
        assert np.array_equal(chi_b[key], chi[key])

# Transient response: chi^neq from the exact dynamics of the oscillator (delta probe + pump) vs Transient_Response with
# the analytic susceptibilities (keys of Xn_frequency_mixing, which include the degeneracy factors)
from scipy.integrate import trapezoid
from mppi.Optics.TransientResponse import envelope_weight
from mppi.Optics.Utils import time_weighted
from mppi.Utilities.Constants import FsToAu

TR_TIME = np.arange(0, 2800, 0.25)
TR_GRID = np.arange(0.15, 0.85 + WP/32, WP/16)   # step wP/16 (Ha)

def chi_neq_exact(osc, pumps, probe_times, w, amplitude=1e-6):
    # (P_Pp - P_P)/E_p at the frequencies w for each pump of the list pumps, minus the linear response of the same
    # integrator
    def chi(data, k):
        t, dP, tp = data.IO_TIME_points, data.Polarization[k][0], data.Efield[k]['initial_time']
        s = t >= tp
        return (dP[s]*np.exp(1j*np.outer(w, t[s]-tp))).sum(1)*(t[1]-t[0])/amplitude
    eq = osc.delta_probe_response(TR_TIME, lambda t: 0., [probe_times[0]], amplitude=amplitude, substeps=2)
    chi1 = chi(eq, 0)
    res = []
    for pump in pumps:
        data = osc.delta_probe_response(TR_TIME, pump, probe_times, amplitude=amplitude, substeps=2)
        res.append([chi(data, k) - chi1 for k in range(len(probe_times))])
    return res

def tr_keys(osc, E0P):
    g = TR_GRID
    if osc.a == 0.:
        return {(1, 0): osc.chi1(g) + E0P**2/4*6*osc.chi3(g, WP, -WP), (1, 2): 3*osc.chi3(g, WP, WP),
                (1, -2): 3*osc.chi3(g, -WP, -WP)}
    return {(1, 1): 2*osc.chi2(g, WP), (1, -1): 2*osc.chi2(g, -WP)}

def test_envelope_weight():
    tau, sigma, t = 1.5, 3., np.array([-1., 0., 1e-8, 2., 5.])
    f2 = np.exp(-tau**2/sigma**2)
    assert np.allclose(envelope_weight(tau, t, sigma, 2, 'probe time'), f2)
    W = envelope_weight(tau, t, sigma, 2)
    assert np.allclose(W[:3], f2)
    tt = np.linspace(0, 5, 20001)
    assert np.isclose(W[-1], trapezoid(np.exp(-(tau+tt)**2/sigma**2), tt)/5, rtol=1e-6)
    W1 = envelope_weight(tau, 5., sigma, 1)
    assert np.isclose(W1, trapezoid(np.exp(-(tau+tt)**2/(2*sigma**2)), tt)/5, rtol=1e-6)
    assert np.allclose(envelope_weight(tau, t, 1e6, 2), 1.)

def test_time_weighted_identity():
    osc = AnharmonicOscillator(omega0=0.5, gamma=0.01)
    chi = osc.chi1(TR_GRID)
    assert np.allclose(time_weighted(TR_GRID, chi, lambda t: np.ones_like(t)), chi)

def test_transient_response_monochromatic():
    tps = [1000., 1031.25]
    # third order, centrosymmetric oscillator
    osc3, E0P = AnharmonicOscillator(omega0=0.5, gamma=0.01, b=0.1), 4e-3
    with quiet():
        tr = O.Transient_Response(TR_GRID*HaToeV, osc3.chi1(TR_GRID), tr_keys(osc3, E0P), E0P, WP*HaToeV)
    w, dchi = tr.dchi(np.array(tps)/FsToAu, E0P)
    exact = chi_neq_exact(osc3, [lambda t: E0P*np.sin(WP*t)], tps, w/HaToeV)[0]
    for k in range(len(tps)):
        assert rel(dchi[k], exact[k]) < 1e-2
    # second order, non-centrosymmetric oscillator: odd part in the pump field
    osc2, E0P = AnharmonicOscillator(omega0=0.5, gamma=0.01, a=0.1), 1e-3
    with quiet():
        tr = O.Transient_Response(TR_GRID*HaToeV, osc2.chi1(TR_GRID), tr_keys(osc2, E0P), E0P, WP*HaToeV)
    w, dchi = tr.dchi(np.array(tps)/FsToAu, E0P)
    plus, minus = chi_neq_exact(osc2, [lambda t: E0P*np.sin(WP*t), lambda t: -E0P*np.sin(WP*t)], tps, w/HaToeV)
    for k in range(len(tps)):
        assert rel(dchi[k], (plus[k]-minus[k])/2) < 1e-2

def test_transient_response_gaussian():
    # pump width sigma = 2/gamma: the envelope at the probe time fails away from the pump maximum, the accumulated
    # envelope reproduces the exact third order response
    osc3, E0P, sigma, tP = AnharmonicOscillator(omega0=0.5, gamma=0.01, b=0.1), 4e-3, 200., 1000.
    taus = np.array([-150., 0., 150.])
    with quiet():
        tr = O.Transient_Response(TR_GRID*HaToeV, osc3.chi1(TR_GRID), tr_keys(osc3, E0P), E0P, WP*HaToeV)
    pump = lambda t: E0P*np.exp(-(t-tP)**2/(2*sigma**2))*np.sin(WP*(t-tP))
    w, acc = tr.dchi(taus/FsToAu, E0P, pump='accumulated', sigma=sigma/FsToAu)
    w, pt = tr.dchi(taus/FsToAu, E0P, pump='probe time', sigma=sigma/FsToAu)
    exact = chi_neq_exact(osc3, [pump], tP+taus, w/HaToeV)[0]
    for k in range(len(taus)):
        assert rel(acc[k], exact[k]) < 0.08
        assert rel(pt[k], exact[k]) > 2*rel(acc[k], exact[k])
    # for a long pulse both the models give the monochromatic result with the envelope at the delay
    w, mono = tr.dchi(taus/FsToAu, E0P)
    for model in ['probe time', 'accumulated']:
        assert np.allclose(tr.dchi(taus/FsToAu, E0P, pump=model, sigma=1e6)[1], mono)

def test_transient_response_options():
    osc3, osc2, E0P = AnharmonicOscillator(b=0.1), AnharmonicOscillator(a=0.1), 1e-3
    keys = {**tr_keys(osc3, E0P), **tr_keys(osc2, E0P)}
    chi1 = osc3.chi1(TR_GRID)
    with quiet():
        tr = O.Transient_Response(TR_GRID*HaToeV, chi1, keys, E0P, WP*HaToeV)
        tr_c = O.Transient_Response(TR_GRID*HaToeV, chi1, keys, E0P, WP*HaToeV, centrosymmetric=True)
    assert tr.orders == [2, 3] and tr_c.orders == [3]
    assert set(tr.components()) == {'chi_0', 'chi_+', 'chi_-', 'chi2_+', 'chi2_-'}
    assert set(tr_c.components()) == {'chi_0', 'chi_+', 'chi_-'}
    # the output grid is at 2wP from the edges and chi_+(w) = chi_(1,2)(w-2wP)
    assert np.allclose(tr.w, TR_GRID[32:-32]*HaToeV)
    assert np.allclose(tr.components()['chi_+'], 3*osc3.chi3(TR_GRID[:-64], WP, WP))
    assert np.allclose(tr.components()['chi_0'], 6*osc3.chi3(TR_GRID[32:-32], WP, -WP))
    # the total is the sum of the orders, the centrosymmetric option drops the second order terms
    taus = [0.3, 1.1]
    total = tr.dchi(taus, E0P)[1]
    assert np.allclose(total, tr.dchi(taus, E0P, orders=[2])[1] + tr.dchi(taus, E0P, orders=[3])[1])
    assert np.allclose(tr_c.dchi(taus, E0P)[1], tr.dchi(taus, E0P, orders=[3])[1])
    with pytest.raises(ValueError):
        tr_c.dchi(taus, E0P, orders=[2])
    with pytest.raises(ValueError):
        tr_c.dchi(taus, E0P, pump='accumulated')          # sigma is needed
    # the accumulated envelope is defined only for the third order terms
    with pytest.raises(ValueError):
        tr.dchi(taus, E0P, pump='accumulated', sigma=3.)
    assert np.allclose(tr.dchi(taus, E0P, pump='accumulated', sigma=3., orders=[3])[1],
                       tr_c.dchi(taus, E0P, pump='accumulated', sigma=3.)[1])
    with pytest.raises(ValueError):                        # the step does not divide the pump frequency
        O.Transient_Response(TR_GRID*HaToeV, chi1, keys, E0P, 1.1*WP*HaToeV, verbose=False)
    with pytest.raises(ValueError):                        # no keys of the transient response
        O.Transient_Response(TR_GRID*HaToeV, chi1, {(1, 0): chi1}, E0P, WP*HaToeV, verbose=False)

def test_transient_response_from_Xn(osc3):
    # components from the keys extracted from the simulated runs (probe step wP/4)
    EP = 2e-3
    w = np.arange(0.30, 0.75, WP/4)
    pump = {'frequency': WP, 'amplitude': EP, 'initial_time': T0}
    data = osc3.sine_response(TIME, w, amplitude=1e-5, initial_time=T0, pump=pump, substeps=4)
    probe = osc3.sine_response(TIME, w, amplitude=1e-5, initial_time=T0, substeps=4)
    with quiet():
        mix = O.Xn_frequency_mixing(data, X_order=(1, 3), verbose=False)
        single = O.Xn_single_frequency(probe, X_order=1, verbose=False)
        tr = O.Transient_Response.from_Xn(mix, single, centrosymmetric=True)
    assert np.isclose(tr.E0P, EP) and np.isclose(tr.wP, WP*HaToeV) and tr.orders == [3]
    wo, comp = tr.w/HaToeV, tr.components()
    assert rel(comp['chi_0'], 6*osc3.chi3(wo, WP, -WP)) < 1e-2
    assert rel(comp['chi_+'], 3*osc3.chi3(wo-2*WP, WP, WP)) < 1e-2
    assert rel(comp['chi_-'], 3*osc3.chi3(wo+2*WP, -WP, -WP)) < 1e-2
