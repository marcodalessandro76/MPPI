"""
This module defines the tools to extract the non-linear optical properties of the system from the real-time polarization.
The module can be loaded in the notebook as follows

>>> from mppi.Utilities import NLoptics as NL

>>> LR.Nonlinear_Response

"""
import numpy as np
from mppi.Utilities import Constants as C
from mppi.Parsers import YamboNLDBParser

def fit_single_frequency(t, y, omega, rcond=None):
    """
    Fit the data to a sine function A*sin(omega*t+phi) at a given frequency omega. 
    
    Args:
        t (:py:class:`numpy.ndarray`): array with the time values   
        y (:py:class:`numpy.ndarray`): array with the values of the function 
        omega (:py:class:`float`): angular frequency of the sine function
        rcond (:py:class:`float`): relative condition number for the least squares fit. Default is None.
    
        Returns:
        :py:class:`tuple` : tuple with the values of the amplitude A and the phase phi of the fitted sine function
    """ 

    X = np.column_stack([
        np.sin(omega * t),
        np.cos(omega * t)])
    
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=rcond)
    alpha, beta = coeffs
    
    A = np.sqrt(alpha**2 + beta**2)
    phi = np.arctan2(beta, alpha)
    
    return A, phi

def fit_two_harmonics(t, y, omega, rcond=None):
    """
    Fit the data to the function A1*sin(omega*t+phi1) + A2*sin(2*omega*t+phi2) at a given frequency omega.
    
    Args:
        t (:py:class:`numpy.ndarray`): array with the time values
        y (:py:class:`numpy.ndarray`): array with the values of the function
        omega (:py:class:`float`): angular frequency of the sine functions
        rcond (:py:class:`float`): relative condition number for the least squares fit. Default is None.
    
    Returns:
        :py:class:`tuple` : tuple with the values of the amplitudes A1, A2 and the phases phi1, phi2 of the fitted sine functions
    """
    X = np.column_stack([
        np.sin(omega * t),
        np.cos(omega * t),
        np.sin(2 * omega * t),
        np.cos(2 * omega * t)])
    
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    alpha1, beta1, alpha2, beta2 = coeffs
    
    # Ricostruzione parametri
    A1 = np.sqrt(alpha1**2 + beta1**2)
    phi1 = np.arctan2(beta1, alpha1)
    
    A2 = np.sqrt(alpha2**2 + beta2**2)
    phi2 = np.arctan2(beta2, alpha2)
    
    return A1, phi1, A2, phi2

def fit_multiple_harmonics(t, y, omega, n_harmonics = 1, rcond=None):
    """
    Fit the data to the function sum_{n=1}^{n_harmonics} A_n*sin(n*omega*t+phi_n) at a given frequency omega.
    
    Args:
        t (:py:class:`numpy.ndarray`): array with the time values
        y (:py:class:`numpy.ndarray`): array with the values of the function
        omega (:py:class:`float`): angular frequency of the sine functions
        n_harmonics (:py:class:`int`): number of harmonics to fit. Default is 1
        rcond (:py:class:`float`): relative condition number for the least squares fit. Default is None
    
    Returns:
        :py:class:`tuple` : tuple with the values of the amplitudes A_n and the phases phi_n of the fitted sine functions. 
        Each component of the tuple is an array with n_harmonics values corresponding to the n harmonics.
    """
    X = np.column_stack([
        np.sin(n * omega * t) for n in range(1, n_harmonics + 1)] +
        [np.cos(n * omega * t) for n in range(1, n_harmonics + 1)])
    
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=rcond)
    
    A = np.zeros(n_harmonics)
    phi = np.zeros(n_harmonics)
    
    for n in range(n_harmonics):
        alpha = coeffs[n]
        beta = coeffs[n + n_harmonics]
        
        A[n] = np.sqrt(alpha**2 + beta**2)
        phi[n] = np.arctan2(beta, alpha)
    
    return A, phi   

class Xn_from_sine():
    """
    Class to extract the non-linear susceptibility (up to the second order) from the polarization induced by a sine-shaped external field.

    Args:
        data (:py:class:`YamboNLDBParser`) : data parsed from the nlndb.Nonlineardatabase database
        X_order (:py:class:`int`) : order of the non-linear susceptibility to be extracted. Can be 1 or 2. Default is 2.
        Trange (:py:class:`list`) : list with the time range in which the polarization is sampled to compute the sine fit. 
            Default is [-1,-1], which means that a time sampling of one period of the external field is used
            (Tstart = time[-1] - Tperiod, Tend = time[-1])
            Further details are provided in the documentation of the method :py:meth:`set_time_sampling`
        Trange_units (:py:class:`str`) : set the units to time sampling. Default is 'fs' and the other possible choice is 'au'.
        tol (:py:class:`float`) : tolerance for the fit of the sine function. Default is 1e-10
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

    def __init__(self,data,X_order=2,Trange=[-1, -1],Trange_units='fs',tol=1e-10,verbose=True):
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
        self.Trange_units = Trange_units
        self.tol = tol
        EFIELDS = ['SIN','SOFTSIN']
        if self.fields[0]['name'] not in EFIELDS:
            raise ValueError(f'Invalid electric field for frequency mixing analysis. Expected one of: {EFIELDS}')
        if verbose:
            self.get_info()
    
    @classmethod
    def from_file(cls, file, verbose = True):
        """
        Initialize the class from a nlndb.Nonlineardatabase database file. The file is parsed through the :py:class:`YamboNLDBParser` class.    
        Args:
            file (:py:class:`string`) : name of the nlndb.Nonlineardatabase database file
            verbose (:py:class:`boolean`) : define the amount of information provided on terminal"""
        data = YamboNLDBParser(file)
        return cls(data,verbose)
    
    def get_info(self):
        """
        Provide information on the keys structure of the instance of the class
        """
        print('Time range of the simulation (in au):',self.time[0],'-',self.time[-1])
        print('Time range of the simulation (in fs):',self.time[0]/C.FsToAu,'-',self.time[-1]/C.FsToAu)
        print('Type of the external field:',self.fields[0]['name'])
        print('Number of frequencies:',self.nfreqs)
        print('Frequency range of the field (in Hartree):',self.fields_freqs[0],'-',self.fields_freqs[-1])
        print('Frequency range of the field (in eV):',self.fields_freqs[0]*C.HaToeV,'-',self.fields_freqs[-1]*C.HaToeV)
        print('Damping factor:',self.damp,'Hartree',self.damp*C.HaToeV,'eV')
        print('Dephasing time:',self.deph,'au',self.deph/C.FsToAu,'fs')
    
    def set_time_sampling(self,ifreq):
        """
        Define the time interval in which the polarization is sampled to compute the sine fit.
        The default time range is one period of the external field (T_start = time[-1] - T_period, T_end = time[-1]). 
        Instead if positive values of Trange are provided the sampling time range is set to that values. Checks are performed to 
        ensure that the time range is consistent with the time sampling of the polarization. Lastly, if the time sampling interval 
        starts before the dephasing time a warning is raised since the fit of the sine function can be not accurate.
        Args:
            ifreq (:py:class:`int`) : index of the frequency of the external field
        
        Returns:
            :py:class:`tuple` : tuple with the indexes of the start and end time of the sampling interval
        """
        if self.Trange_units == 'fs':
            time = self.time/C.FsToAu
            dt = self.dt/C.FsToAu
            Tperiod = 2.0*np.pi/self.fields_freqs[ifreq]/C.FsToAu
            deph = self.deph/C.FsToAu
        elif self.Trange_units == 'au':
            time = self.time
            dt = self.dt
            Tperiod = 2.0*np.pi/self.fields_freqs[ifreq]
            deph = self.deph
        else:
            raise ValueError("Invalid time units. Please use 'fs' or 'au'.")
        
        if self.Trange[0] < 0.:     
            Tstart = time[-1] - Tperiod
        else:
            Tstart = self.Trange[0]
        if Tstart < deph:
            print('Warning: the time sampling starts before the dephasing time. The fit of the sine function can be not accurate.') 
        if Tstart > time[-1]:
            Tstart = time[-1] - Tperiod
            print('Warning: the time sampling starts after the end of the time range. Tstart is set to time[-1]-Tperiod.') 
        
        if self.Trange[1] < 0.:   
            Tend = time[-1]   
        else:            
            Tend = self.Trange[1]
        if Tend > time[-1]:
            Tend = time[-1]
            print('Warning: the time sampling ends after the end of the time range. Tend is set to time[-1].') 
        if Tend < (Tstart + Tperiod):
            Tend = Tstart + Tperiod
            print('Warning: the time sampling ends before it starts. Tend is set to Tstart + Tperiod.')
        
        iTstart = int(np.round( Tstart / dt))
        iTend = int(np.round( Tend / dt)) 
        
        return iTstart, iTend
    
    def perform_harm_analysis(self,ifreq):
        """
        Perform harmonic analysis of the polarization using the fit_multiple_harmonics function.

        Args:
            ifreq (:py:class:`int`) : index of the frequency of the external field
        
        Returns:
            :py:class:`tuple` : tuple with the values of the amplitude A and the phase phi of the fitted sine function.
                Each component of the tuple is an array with 3 values corresponding to the 3 cartesian directions.
        """
        X_order = self.X_order
        iTstart, iTend = self.set_time_sampling(ifreq)
        omega = self.fields_freqs[ifreq]

        A, phi = np.zeros((X_order,3)), np.zeros((X_order,3))
        for idir in range(3):
            A[:,idir], phi[:,idir] = fit_multiple_harmonics(self.time[iTstart:iTend], self.pol[ifreq,idir,iTstart:iTend], omega, X_order, rcond=self.tol)
        return A, phi
    
    def eval_pol_freq(self):
        """
        Compute the polarization in the frequency domain at the (multiples) of the frequencies of the external fields according to the value of X_order
        
        Returns:
            :py:class:`numpy.ndarray` : array with the polarization in the frequency domain at the (multiple) frequencies of the external fields according 
                to the value of X_order. The shape of the array is (X_order,nfreqs,3) where 3 is the number of cartesian directions 
        """
        X_order = self.X_order
        pol_w = np.zeros((X_order,self.nfreqs,3),dtype=complex)
        for ifreq in range(self.nfreqs):
            A, phi = self.perform_harm_analysis(ifreq)
            for n_harm in range(1,X_order+1):
                for idir in range(3):
                    pol_w[n_harm-1,ifreq,idir] = 1j*A[n_harm-1,idir]/2.0*np.exp(-1j*phi[n_harm-1,idir])
            
        return pol_w
    
    def eval_Efield_freq(self):
        """
        Evaluate the external fields in the frequency domain for all the values of the self.fields_freqs array. 
        
        Returns:
            :py:class:`numpy.ndarray` : array with the external fields in the frequency domain for all the values of the self.fields_freqs array 
                  
        """
        efield_w = np.zeros(self.nfreqs,dtype=complex)
        for ifreq in range(self.nfreqs):
            E0 = self.fields[ifreq]['amplitude']
            t0 = self.fields[ifreq]['initial_time']
            omega = self.fields_freqs[ifreq]
            efield_w[ifreq] = 1j*E0/2.0*np.exp(1j*t0*omega)
        return efield_w
    
    def compute_Xn(self):
        """
        Compute the non-linear susceptibility of the system at the frequencies of the external fields. 
        
        Returns:
            :py:class:`numpy.ndarray` : array with the non-linear susceptibility of the system at the frequencies of the external fields. The shape of the array is (3,nfreqs) where 3 is the number of cartesian directions and nfreqs is the number of frequencies of the external fields.
        """
        pol_w = self.eval_pol_freq()
        efield_w = self.eval_Efield_freq()
        Xn = np.zeros((3,self.nfreqs),dtype=complex)
        for ifreq in range(self.nfreqs):
            for idir in range(3):
                Xn[idir,ifreq] = pol_w[idir,ifreq]/efield_w[ifreq]
        return Xn

