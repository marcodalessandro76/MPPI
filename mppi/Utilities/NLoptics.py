"""
This module defines the tools to extract the non-linear optical properties of the system from the real-time polarization.
The module can be loaded in the notebook as follows

>>> from mppi.Utilities import NLoptics as NL

>>> LR.Nonlinear_Response

"""
import numpy as np
from mppi.Utilities import Constants as C
from mppi.Parsers import YamboNLDBParser

def fit_sine_with_offset(t, y, omega, rcond=None):
    """
    Fit a sine function with an offset C0+A*sin(omega*t+phi) to the data y(t) at a given frequency omega. 
    
    Args:
        t (:py:class:`numpy.ndarray`): array with the time values   
        y (:py:class:`numpy.ndarray`): array with the values of the function 
        omega (:py:class:`float`): frequency of the sine function to fit in rad
        rcond (:py:class:`float`): relative condition number for the least squares fit. Default is None.
    
        Returns:
        :py:class:`tuple` : tuple with the values of the offset C0, the amplitude A and the phase phi of the fitted sine function
    """ 

    X = np.column_stack([
        np.ones_like(t),         
        np.sin(omega * t),
        np.cos(omega * t)])
    
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=rcond)
    C0, alpha, beta = coeffs
    
    A = np.sqrt(alpha**2 + beta**2)
    phi = np.arctan2(beta, alpha)
    
    return C0, A, phi

class Xn_from_sine():
    """
    Class to extract the non-linear susceptibility (up to the second order) from the polarization induced by a sine-shaped external field.

    Args:
        data (:py:class:`YamboNLDBParser`) : data parsed from the nlndb.Nonlineardatabase database
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
        Trange (:py:class:`list`) : list with the time range in which the polarization is sampled to compute the sine fit. 
            Default is [-1,-1], which means that a time sampling of one period of the external field is used
            (Tstart = time[-1] - Tperiod, Tend = time[-1])
            Further details are provided in the documentation of the method :py:meth:`set_time_sampling`
        Trange_units (:py:class:`str`) : set the units to time sampling. Default is 'fs' and the other possible choice is 'au'
        tol (:py:class:`float`) : tolerance for the fit of the sine function. Default is 1e-10
    """

    def __init__(self,data,Trange=[-1, -1],Trange_units='fs',tol=1e-10,verbose=True):
        self.time = data.IO_TIME_points
        self.pol = np.array(data.Polarization) 
        self.fields = data.Efield # in this case the variables Efield2 and Efield_general contains the same information to be checked in general
        self.nfreqs = data.n_frequencies
        self.fields_freqs = np.array([e['freq_range'][0] for e in self.fields]) 
        #self.field_time0 = np.array([e['initial_time'] for e in self.fields]) 
        self.damp = data.NL_damping
        self.deph = 12/self.damp
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
    
    def eval_harmonic_fit(self,ifreq):
        """
        Evaluate the amplitude and phase of the polarization at the frequency of the external field using the 
        fit_sine_with_offset function.

        Args:
            ifreq (:py:class:`int`) : index of the frequency of the external field
        
        Returns:
            :py:class:`tuple` : tuple with the values of the offset C0, the amplitude A and the phase phi of the fitted sine function.
                Each component of the tuple is an array with 3 values corresponding to the 3 cartesian directions.
        """
        iTstart, iTend = self.set_time_sampling(ifreq)
        omega = self.fields_freqs[ifreq]

        C0, A, phi = np.zeros(3), np.zeros(3), np.zeros(3) 
        for idir in range(3):
             C0[idir], A[idir], phi[idir] = fit_sine_with_offset(self.time[iTstart:iTend], self.pol[ifreq,idir,iTstart:iTend], omega, rcond=self.tol)
        return C0, A, phi
    
    def eval_harm_pol_w(self):
        """
        Compute the polarization in the frequency domain at the frequencies of the external fields using the fit of the sine function. 
        
        Returns:
            :py:class:`numpy.ndarray` : array with the polarization in the frequency domain at the frequencies of the external fields. 
                The shape of the array is (nfreqs,3) where nfreqs is the number of frequencies of the external fields and 3 is the number of cartesian directions.
        """
        pol_w = np.zeros((self.nfreqs,3),dtype=complex)
        for ifreq in range(self.nfreqs):
            _, A, phi = self.eval_harmonic_fit(ifreq)
            for idir in range(3):
                pol_w[ifreq,idir] = 1j*A[idir]/2.0*np.exp(-1j*phi[idir])
        return pol_w
    
    def eval_Efield_w(self):
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

