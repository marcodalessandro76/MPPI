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
from mppi.Optics.Utils import fit_sum_frequencies, eval_sum_frequencies
from mppi.Parsers import YamboNLDBParser


# def fit_frequency_mixing(t, y, omega1, omega2, mixing_terms=None, rcond=None):
#     """
#     Fit:
#         f(t) = B0 + sum_k A_k sin(Omega_k t + phi_k)
    
#     where Omega_k are mixed frequencies:
#         Omega_k = a*omega1 + b*omega2
    
#     Args:
#         t, y (:py:class:`numpy.ndarray`): data
#         omega1, omega2 (:py:class:`float`): base frequencies
#         mixing_terms (:py:class:`list` or None): list of (a, b) tuples defining frequencies
#                       Omega = a*omega1 + b*omega2
#                       Default includes desired mixing terms
#         rcond (:py:class:`float` or None): lstsq parameter
    
#     Returns:
#         Omegas (:py:class:`numpy.ndarray`): array of fitted frequencies
#         A (:py:class:`numpy.ndarray`): amplitudes
#         phi (:py:class:`numpy.ndarray`): phases
#         B0 (:py:class:`float`): offset
#         residuals (:py:class:`numpy.ndarray`): residuals
#     """
    
#     if mixing_terms is None:
#         mixing_terms = [
#             (1, 0),   # omega1
#             (1, 1),   # omega1 + omega2
#             (1, -1),  # omega1 - omega2
#             (1, 2),   # omega1 + 2 omega2
#             (1, -2),  # omega1 - 2 omega2
#         ]
    
#     Omegas = np.array([a*omega1 + b*omega2 for (a, b) in mixing_terms])
    
#     X_cols = [np.ones_like(t)]  # offset
#     for Omega in Omegas:
#         X_cols.append(np.sin(Omega * t))
#     for Omega in Omegas:
#         X_cols.append(np.cos(Omega * t))
#     X = np.column_stack(X_cols)
#     coeffs, residuals, _, _ = np.linalg.lstsq(X, y, rcond=rcond)
    
#     B0 = coeffs[0]
#     n_terms = len(Omegas)
#     A = np.zeros(n_terms)
#     phi = np.zeros(n_terms)

#     for i in range(n_terms):
#         alpha = coeffs[1 + i]
#         beta  = coeffs[1 + i + n_terms]
        
#         A[i] = np.sqrt(alpha**2 + beta**2)
#         phi[i] = np.arctan2(beta, alpha)
    
#     return Omegas, A, phi, B0, np.sqrt(residuals)

# def generate_frequencies(omega1, omega2,max_order_E2=3,include_pure_E1=True,include_pure_E2=True,include_mixing=True,return_mapping=False,tol=1e-8):
#     """
#     Generate positive frequencies for weak-strong field expansion. The rationale behind the choice of producing positive frequencies only is that the
#     only the positive frequency of the probe field is considered in the expansion of the polarization, so if for instance omega1-m*omega2 is negative, 
#     it means that there is a positve contribution at m*omega2-omega1 which is not explicitly produced in the combination of the coefficients but is selected
#     by the check on the positive frequencies performed here.

#     Args:
#         omega1 (:py:class:`float`): frequency of the first field (pump)
#         omega2 (:py:class:`float`): frequency of the second field (probe)
#         max_order_E2 (:py:class:`int`): maximum order for the second field
#         include_pure_E1 (:py:class:`bool`): whether to include pure E1 terms
#         include_pure_E2 (:py:class:`bool`): whether to include pure E2 terms
#         include_mixing (:py:class:`bool`): whether to include mixing terms
#         return_mapping (:py:class:`bool`): whether to return the mapping
#         tol (:py:class:`float`): tolerance for frequency comparison

#     Returns:
#         Omegas (:py:class:`numpy.ndarray`): array of generated frequencies
#         mapping (:py:class:`dict`, optional): Omega -> list of contributions (n, m)
#     """
    
#     freq_map = {}

#     def add_freq(n, m):
#         Omega = abs(n*omega1 + m*omega2)
#         if Omega < tol:
#             return   
#         key = np.round(Omega,6)
#         if key not in freq_map:
#             freq_map[key] = {"Omega": Omega,"contributions": []}
        
#         freq_map[key]["contributions"].append({"n": n,"m": m})

#     if include_pure_E1:
#         add_freq(1, 0)
#     if include_pure_E2:
#         for m in range(1, max_order_E2 + 1):
#             add_freq(0, m)
#     if include_mixing:
#         for m in range(1, max_order_E2 + 1):
#             add_freq(1,  m)
#             add_freq(1, -m)

#     Omegas = np.array(sorted(freq_map.keys()))

#     if return_mapping:
#         return Omegas, freq_map
#     else:
#         return Omegas

def generate_frequencies(omega1, omega2,max_order_E2=3,include_pure_E1=True,include_pure_E2=True,include_mixing=True,tol=1e-8):
    """
    Generate positive frequencies as a dictionary:
        (n, m) -> Omega = |n*omega1 + m*omega2|

    Only positive frequencies are kept. The rationale behind the choice of producing positive frequencies only is that the
    only the positive frequency of the probe field is considered in the expansion of the polarization, so if for instance omega1-m*omega2 is negative, 
    it means that there is a positve contribution at m*omega2-omega1 which is not explicitly produced in the combination of the coefficients but is selected
    by the check on the positive frequencies performed here.

    Args:
        omega1 (:py:class:`float`): frequency of the first field (pump)
        omega2 (:py:class:`float`): frequency of the second field (probe)
        max_order_E2 (:py:class:`int`): maximum order for the second field
        include_pure_E1 (:py:class:`bool`): whether to include pure E1 terms
        include_pure_E2 (:py:class:`bool`): whether to include pure E2 terms
        include_mixing (:py:class:`bool`): whether to include mixing terms
        return_mapping (:py:class:`bool`): whether to return the mapping
        tol (:py:class:`float`): tolerance for frequency comparison

    Returns:
        Omegas_dict: dict with keys (n,m) and values Omega
    """
    
    Omegas_dict = {}

    def add_freq(n, m):
        Omega = n*omega1 + m*omega2
        
        if abs(Omega) < tol:
            return
        
        if Omega < 0:
            n, m = -n, -m
            Omega = -Omega
        key = (n, m)
        for (n2, m2), O2 in Omegas_dict.items():
            if np.isclose(Omega, O2, rtol=1e-6, atol=tol):
                return
        
        Omegas_dict[key] = Omega

    if include_pure_E1:
        add_freq(1, 0)
    if include_pure_E2:
        for m in range(1, max_order_E2 + 1):
            add_freq(0, m)
    if include_mixing:
        for m in range(1, max_order_E2 + 1):
            add_freq(1,  m)
            add_freq(1, -m)

    return Omegas_dict

# def eval_frequency_mixing(t, A, phi, B0, omega):
#     """
#     Evaluate the multiple harmonics function at given time points.

#     Args:
#         t (:py:class:`numpy.ndarray`): array with the time values
#         A (:py:class:`numpy.ndarray`): amplitudes
#         phi (:py:class:`numpy.ndarray`): phases
#         B0 (:py:class:`float`): constant offset
#         omega (:py:class:`float`): angular frequency of the sine functions

#     Returns:
#         :py:class:`numpy.ndarray`: evaluated function values
#     """
#     n_harmonics = len(A)
#     y = B0 * np.ones_like(t)
#     for n in range(n_harmonics):
#         y += A[n] * np.sin((n + 1) * omega * t + phi[n])
#     return y


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
        probe_freqs (:py:class:`numpy.ndarray`) : array with the frequencies of the probes in Hartree
        pump_freq (:py:class:`float`) : frequency of the pump in Hartree
        field_time0 (:py:class:`numpy.ndarray`) : array with the switch on time of the external fields in au
        damp (:py:class:`float`) : damping factor in Hartree
        deph (:py:class:`float`) : dephasing time in au, defined as 12/damp. If the time sampling interval starts before this value
            a warning is raised since the fit of the sine function can be not accurate.
        dt (:py:class:`float`) : time sampling interval in au
    """

    def __init__(self,data,X_order=2,Trange=[-1, -1],Trange_units='fs',tol=1e-10,inactive_harmonics=None,verbose=True):
        self.time = data.IO_TIME_points
        self.pol = np.array(data.Polarization) 
        self.probes = data.Efield
        self.pump = data.Efield_general[1]
        self.nfreqs = data.n_frequencies
        self.probe_freqs = np.array([e['freq_range'][0] for e in self.probes]) 
        self.pump_freq = self.pump['freq_range'][0]
        self.damp = data.NL_damping
        self.deph = 12/self.damp
        self.X_order = X_order
        self.dt = self.time[1]-self.time[0]
        self.Trange = Trange
        self.inactive_harmonics = inactive_harmonics
        self.Trange_units = Trange_units
        self.tol = tol
        EFIELDS = ['SIN','SOFTSIN']
        if self.probes[0]['name'] not in EFIELDS:
            raise ValueError(f'Invalid probe field for frequency mixing analysis. Expected one of: {EFIELDS}')
        if self.pump['name'] not in EFIELDS:
            raise ValueError(f'Invalid pump field for frequency mixing analysis. Expected one of: {EFIELDS}')
        if len(data.Efield_general) > 2 and data.Efield_general[2]['name'] != "none":
            raise ValueError("Only mixing of two fields implemented.")
        if verbose: 
            self.get_info()

    @classmethod
    def from_file(cls, file, verbose = True):
        """
        Initialize the class from a nlndb.Nonlineardatabase database file. The file is parsed through the :py:class:`YamboNLDBParser` class.    
        Args:
            file (:py:class:`string`) : name of the nlndb.Nonlineardatabase database file
            verbose (:py:class:`boolean`) : define the amount of information provided on terminal
        """
        data = YamboNLDBParser(file)
        return cls(data,verbose)
    
    def get_info(self):
        """
        Provide information on the keys structure of the instance of the class
        """
        print('Time range of the simulation (in au):',self.time[0],'-',self.time[-1])
        print('Time range of the simulation (in fs):',self.time[0]/C.FsToAu,'-',self.time[-1]/C.FsToAu)
        print('Type of the probe field:',self.probes[0]['name'])
        print('Type of the pump field:',self.pump['name'])
        print('Order of the non-linear susceptibility extracted:',self.X_order)
        print('Number of frequencies:',self.nfreqs)
        print('Frequency range of the probe (in Hartree):',self.probe_freqs[0],'-',self.probe_freqs[-1])
        print('Frequency range of the probe (in eV):',self.probe_freqs[0]*C.HaToeV,'-',self.probe_freqs[-1]*C.HaToeV)
        print('Frequency of the pump :',self.pump_freq,' Hartree -',self.pump_freq*C.HaToeV,' eV')
        print('Damping factor:',self.damp,'Hartree',self.damp*C.HaToeV,'eV')
        print('Dephasing time:',self.deph,'au',self.deph/C.FsToAu,'fs')