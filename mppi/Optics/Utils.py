"""
Here we collect some static functions used by the classes of the optics module.
"""
import numpy as np

# def fit_sum_frequencies(t, y, Omegas, mapping=None, rcond=None):
#     """
#     Fit the data (t,y) with a function of the form:
#         f(t) = B0 + sum_k A_k sin(Omega_k t + phi_k)

#     Args:
#         t, y (:py:class:`numpy.ndarray`): data
#         Omegas (:py:class:`numpy.ndarray`): array of frequencies (positive, no duplicates)
#         mapping (:py:class:`dict`, optional): Omega -> contributions info
#         rcond (:py:class:`float`, optional): lstsq parameter

#     Returns:
#         results (list of dict, if mapping is provided) OR
#         (A, phi, B0, residuals) otherwise
#     """
    
#     Omegas = np.asarray(Omegas)
#     n_terms = len(Omegas)

#     X_cols = [np.ones_like(t)]
#     for O in Omegas:
#         X_cols.append(np.sin(O * t))
#     for O in Omegas:
#         X_cols.append(np.cos(O * t))

#     X = np.column_stack(X_cols)
#     coeffs, residuals, _, _ = np.linalg.lstsq(X, y, rcond=rcond)

#     B0 = coeffs[0]
#     A = np.zeros(n_terms)
#     phi = np.zeros(n_terms)

#     results = []

#     for i, O in enumerate(Omegas):
#         alpha = coeffs[1 + i]
#         beta  = coeffs[1 + i + n_terms]

#         A[i] = np.sqrt(alpha**2 + beta**2)
#         phi[i] = np.arctan2(beta, alpha)

#         if mapping is not None:
#             key = np.round(O, 12)
#             info = mapping.get(key, {"contributions": None})

#             results.append({
#                 "Omega": O,
#                 "A": A[i],
#                 "phi": phi[i],
#                 "contributions": info["contributions"]
#             })

#     if mapping is not None:
#         return results, B0, np.sqrt(residuals)
#     else:
#         return A, phi, B0, np.sqrt(residuals)
    
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

# def eval_sum_frequencies(t, A, phi, B0, Omegas):
#     """
#     Evaluate the sum of frequencies function 
#         y(t) = B0 + sum_k A_k sin(Omega_k t + phi_k)

#     at given time points.

#     Args:
#         t (:py:class:`numpy.ndarray`): array with the time values
#         A (:py:class:`numpy.ndarray`): amplitudes
#         phi (:py:class:`numpy.ndarray`): phases
#         B0 (:py:class:`float`): constant offset
#         omega (:py:class:`float`): angular frequency of the sine functions

#     Returns:
#         :py:class:`numpy.ndarray`: evaluated function values
#     """
#     n_frequencies = len(A)
#     y = B0 * np.ones_like(t)
#     for n in range(n_frequencies):
#         y += A[n] * np.sin(Omegas[n] * t + phi[n])
#     return y

def eval_sum_frequencies(t, results_dict, B0):
    """
    Evaluate:
        y(t) = B0 + sum_{key} A_key sin(Omega_key t + phi_key)

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