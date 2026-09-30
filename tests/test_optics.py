import os
import numpy as np
import pytest
from mppi.Optics.Utils import fit_sum_frequencies, eval_sum_frequencies

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

def test_xn_single_frequency_from_file(ref_dir):
    from mppi.Optics import Xn_single_frequency
    folder = os.path.join(ref_dir,'nl_results','LiF-sine_pulse-time_200fs-step_10as')
    if not os.path.isdir(folder): pytest.skip('nl_results reference data not available')
    xn = Xn_single_frequency.from_file(os.path.join(folder,'ndb.Nonlinear'),verbose=False)
    assert xn.X_order == 3
