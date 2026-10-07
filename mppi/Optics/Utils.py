"""
Here we collect some static functions used by the classes of the optics module.
"""
import numpy as np
    
def fit_sum_frequencies(t, y, Omegas_dict, rcond=None):
    """
    Fit the data (t,y) with a function of the form:
        f(t) = B0 + sum_{k} A_k sin(Omega_k t + phi_k)

    Args:
        t, y (:py:class:`numpy.ndarray`): data
        Omegas_dict (:py:class:`dict`): {k: Omega}
        rcond (:py:class:`float`, optional): lstsq parameter

    Returns:
        results_dict: {k: {'Omega', 'A', 'phi', 'alpha', 'beta'}}
        :py:class:`numpy.float64`: constant offset
        :py:class:`numpy.float64`: residuals
    """

    keys = list(Omegas_dict.keys())
    Omegas = np.array([Omegas_dict[k] for k in keys])
    n_terms = len(Omegas)

    X_cols = [np.ones_like(t)]
    for O in Omegas:
        X_cols.append(np.sin(O * t))
    for O in Omegas:
        X_cols.append(np.cos(O * t))
    X = np.column_stack(X_cols)

    coeffs, residuals, _, _ = np.linalg.lstsq(X, y, rcond=rcond)
    # lstsq returns an empty residuals array if the system is rank deficient or underdetermined,
    # in this case the residual is computed explicitly
    if residuals.size == 0:
        residuals = np.array([np.sum((y - X @ coeffs)**2)])
    B0 = coeffs[0]
    results_dict = {}

    for i, key in enumerate(keys):
        O = Omegas[i]

        alpha = coeffs[1 + i]
        beta  = coeffs[1 + i + n_terms]

        A = np.sqrt(alpha**2 + beta**2)
        phi = np.arctan2(beta, alpha)

        results_dict[key] = {
            "Omega": O,
            "A": A,
            "phi": phi,
        }

    return results_dict, B0, np.sqrt(residuals)[0]

def lorentzian_broadening(freqs, chi, delta, renormalize=True):
    r"""
    Increase the broadening of a response function sampled on a grid of frequencies by a Lorentzian convolution.
    A causal response function is analytic in the upper half of the complex frequency plane, so

    .. math::
        \chi(\omega+i\Delta) = \frac{1}{\pi}\int d\omega'\,\frac{\Delta}{(\omega-\omega')^2+\Delta^2}\,\chi(\omega')

    i.e. the convolution evaluates the response at the complex frequency :math:`\omega+i\Delta`. If the frequency enters
    once in each resonance denominator, as in the linear susceptibility and in the susceptibilities linear in the probe
    field of a frequency mixing analysis, this adds :math:`\Delta` to the width of all the resonances that involve it. The
    convolution is not a physical broadening for the harmonics :math:`n\geq2` of a monochromatic field
    (:math:`\omega\to\omega+i\Delta` adds :math:`n\Delta` to the width of the n-photon resonance) and it is not defined for
    the terms that contain both :math:`E(\omega)` and :math:`E(-\omega)` (e.g. the optical rectification).

    The integral is computed with the trapezoidal rule on the grid, so the step of the grid must be much smaller than
    delta. The Lorentzian tails outside the frequency window are missing, with an error of order
    :math:`(\Delta/\pi)(1/d_1+1/d_2)`, where :math:`d_1` and :math:`d_2` are the distances from the edges of the window
    (so it affects the whole window, not only the frequencies close to its edges).
    With renormalize=True the kernel is normalized to one on the grid, which replaces the missing tails with the values
    inside the window (it overestimates the peaks close to the edges if the response decreases outside the window).

    Args:
        freqs (:py:class:`numpy.ndarray`): increasing frequencies of the grid
        chi (:py:class:`numpy.ndarray`): values of the response function on the grid (the last axis runs over the frequencies)
        delta (:py:class:`float`): broadening, in the same units of freqs
        renormalize (:py:class:`bool`): if True the kernel is normalized to one on the grid. Default is True

    Returns:
        :py:class:`numpy.ndarray`: the broadened response function on the same grid
    """
    w = np.asarray(freqs, dtype=float)
    if len(w) < 2:
        raise ValueError('The Lorentzian broadening needs at least two frequencies')
    if np.any(np.diff(w) <= 0.):
        raise ValueError('The frequencies of the Lorentzian broadening must be increasing')
    # trapezoidal weights, valid also for a non uniform grid
    weights = np.empty(len(w))
    weights[1:-1] = (w[2:] - w[:-2]) / 2.
    weights[0], weights[-1] = (w[1] - w[0]) / 2., (w[-1] - w[-2]) / 2.
    kernel = delta / np.pi / ((w[:, None] - w[None, :])**2 + delta**2) * weights[None, :]
    if renormalize:
        kernel /= kernel.sum(axis=1)[:, None]
    return np.asarray(chi) @ kernel.T

def print_sampling_warnings(warnings, nfreqs, max_shown=10):
    """
    Print a summary of the warnings raised in the choice of the time sampling of the harmonic analysis. Each warning is
    printed once, together with the number and the indexes of the frequencies of the external field that raised it.

    Args:
        warnings (:py:class:`dict`): {ifreq: list of warning messages}
        nfreqs (:py:class:`int`): total number of frequencies of the external field
        max_shown (:py:class:`int`): maximum number of frequency indexes printed for each warning
    """
    summary = {}
    for ifreq in sorted(warnings):
        for msg in warnings[ifreq]:
            summary.setdefault(msg, []).append(ifreq)
    for msg, ifreqs in summary.items():
        shown = ', '.join(str(i) for i in ifreqs[:max_shown])
        if len(ifreqs) > max_shown:
            shown += ', ...'
        print(f'Warning: {msg} ({len(ifreqs)} of {nfreqs} frequencies, indexes: {shown})')

def eval_sum_frequencies(t, results_dict, B0):
    """
    Evaluate:
        y(t) = B0 + sum_{key} A_key sin(Omega_key t + phi_key)
    
    given the results of fit_sum_frequencies.

    Args:
        t (:numpy.ndarray): time array
        results_dict (:py:class:`dict`): output of fit_sum_frequencies
        B0 (:py:class:`float`): constant offset

    Returns:
        y (:numpy.ndarray): evaluated function
    """
    
    y = B0 * np.ones_like(t)
    for key, res in results_dict.items():
        Omega = res["Omega"]
        A = res["A"]
        phi = res["phi"]

        y += A * np.sin(Omega * t + phi)

    return y
