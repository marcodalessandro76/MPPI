"""
This module defines the tools to extract the non-linear optical properties of the system from the real-time polarization induced by the sum of two monocromatic external fields,
that represent the pump and the probe in a typical pump-probe experiment. 
The module can be loaded in the notebook as follows

>>> from mppi.Optics import Xn_frequency_mixing

or the class can be imported directly as

>>> from mppi.Optics.Xn_frequency_mixing import Xn_frequency_mixing

>>> Xn_frequency_mixing(data,Trange=[100,-1],X_order=2)

"""
import numpy as np
from mppi.Utilities import Constants as C
from mppi.Utilities import Utils as U
from mppi.Parsers import YamboNLDBParser


def fit_frequency_mixing(
    t, y,
    omega1, omega2,
    mixing_terms=None,
    rcond=None
):
    """
    Fit:
        f(t) = B0 + sum_k A_k sin(Omega_k t + phi_k)
    
    where Omega_k are mixed frequencies:
        Omega_k = a*omega1 + b*omega2
    
    Args:
        t, y: data
        omega1, omega2: base frequencies
        mixing_terms: list of (a, b) tuples defining frequencies
                      Omega = a*omega1 + b*omega2
                      Default includes desired mixing terms
        rcond: lstsq parameter
    
    Returns:
        Omegas: array of fitted frequencies
        A: amplitudes
        phi: phases
        B0: offset
        residuals
    """
    
    # ---- default: i tuoi termini ----
    if mixing_terms is None:
        mixing_terms = [
            (1, 0),   # omega1
            (1, 1),   # omega1 + omega2
            (1, -1),  # omega1 - omega2
            (1, 2),   # omega1 + 2 omega2
            (1, -2),  # omega1 - 2 omega2
        ]
    
    # ---- costruzione frequenze ----
    Omegas = np.array([a*omega1 + b*omega2 for (a, b) in mixing_terms])
    
    # ---- matrice ----
    X_cols = [np.ones_like(t)]  # offset
    
    for Omega in Omegas:
        X_cols.append(np.sin(Omega * t))
    for Omega in Omegas:
        X_cols.append(np.cos(Omega * t))
    
    X = np.column_stack(X_cols)
    
    # ---- fit ----
    coeffs, residuals, _, _ = np.linalg.lstsq(X, y, rcond=rcond)
    
    # ---- estrazione ----
    B0 = coeffs[0]
    
    n_terms = len(Omegas)
    
    A = np.zeros(n_terms)
    phi = np.zeros(n_terms)
    
    for i in range(n_terms):
        alpha = coeffs[1 + i]
        beta  = coeffs[1 + i + n_terms]
        
        A[i] = np.sqrt(alpha**2 + beta**2)
        phi[i] = np.arctan2(beta, alpha)
    
    return Omegas, A, phi, B0, np.sqrt(residuals)


class Xn_frequency_mixing():
    """
    Class to extract the non-linear susceptibility from the polarization induced by thesum of two monocromatic external fields, that represent the pump and the probe in a typical pump-probe experiment. 
    
    Args:
        data (:py:class:`YamboNLDBParser`) : data parsed from the nlndb.Nonlineardatabase database
        X_order (:py:class:`int`) : order of the non-linear susceptibility to be extracted. The zero-th term corresponds to the constant offset of the polarization harmonic expansion. 
            Default is 2, which means that the non-linear susceptibility up to the second order is extracted
        Trange (:py:class:`list`) : list with the time range in which the polarization is sampled to compute the sine fit 
            Default is [-1,-1], which means that a time sampling of one period of the external field is used
            (Tstart = time[-1] - Tperiod, Tend = time[-1])
            Further details are provided in the documentation of the method :py:meth:`set_time_sampling`
        Trange_units (:py:class:`str`) : set the units to time sampling. Default is 'fs' and the other possible choice is 'au'
        tol (:py:class:`float`) : tolerance for the fit of the sine function. Default is 1e-10
        inactive_harmonics (:py:class:`list` or None) : list with the harmonics to exclude from the fit 
            Example: [2] excludes the 2nd harmonic. Default is None, which means that all the harmonics up to X_order are included in the fit
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal
    
    Attributes:
        time (:py:class:`numpy.ndarray`) : array with the time values in au
        Pol (:py:class:`numpy.ndarray`) : array with the polarization of the system at each time step in au
        fields (:py:class:`list`) : list of dictionaries with the information on the external fields 
        nfreqs (:py:class:`int`) : number of frequencies of the external fields
        fields_freqs (:py:class:`numpy.ndarray`) : array with the frequencies of the external fields in Hartree
        field_time0 (:py:class:`numpy.ndarray`) : array with the switch on time of the external fields in au
        damp (:py:class:`float`) : damping factor in Hartree
        deph (:py:class:`float`) : dephasing time in au, defined as 12/damp. If the time sampling interval starts before this value
            a warning is raised since the fit of the sine function can be not accurate.
        dt (:py:class:`float`) : time sampling interval in au
    """

    def __init__(self,data,X_order=2,Trange=[-1, -1],Trange_units='fs',tol=1e-10,inactive_harmonics=None,verbose=True):
        self.time = data.IO_TIME_points
        self.pol = np.array(data.Polarization) 
        self.fields = data.Efield # in this case the variables Efield2 and Efield_general contains the same information to be checked in general
        self.nfreqs = data.n_frequencies
        self.fields_freqs = np.array([e['freq_range'][0] for e in self.fields]) 
        self.damp = data.NL_damping
        self.deph = 12/self.damp
        self.X_order = X_order
        self.dt = self.time[1]-self.time[0]
        self.Trange = Trange
        self.inactive_harmonics = inactive_harmonics
        self.Trange_units = Trange_units
        self.tol = tol
        EFIELDS = ['SIN','SOFTSIN']
        if self.fields[0]['name'] not in EFIELDS:
            raise ValueError(f'Invalid electric field for frequency mixing analysis. Expected one of: {EFIELDS}')
        if verbose:
            self.get_info()