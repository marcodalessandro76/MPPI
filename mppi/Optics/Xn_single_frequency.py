"""
This module defines the tools to extract the non-linear optical properties of the system from the real-time polarization induced by a sine-shaped (monochromatic) external field. 
The module can be loaded in the notebook as follows

>>> from mppi.Optics import Xn_single_frequency

or the class can be imported directly as

>>> from mppi.Optics.Xn_single_frequency import Xn_single_frequency

>>> Xn_single_frequency(data,X_order=2,Trange=[100,-1],Trange_units='fs')

"""
import numpy as np
from mppi.Utilities import Constants as C
from mppi.Utilities import Utils as U
from mppi.Optics.Utils import fit_sum_frequencies, eval_sum_frequencies
from mppi.Parsers import YamboNLDBParser

# def generate_frequencies(omega,n_harmonics=3,inactive_harmonics=None,return_mapping=False,tol=1e-8):
#     """
#     Generate positive frequencies for harmonic expansion with a single monochromatic field.

#     Args:
#         omega (:py:class:`float`): frequency of the first field (pump)
#         n_harmonics (:py:class:`int`): number of harmonics to include in the expansion. Default is 3
#         inactive_harmonics (:py:class:`list` or None): list with the harmonics to exclude from the expansion 
#             Example: [2] excludes the 2nd harmonic. Default is None
#         return_mapping (:py:class:`bool`): whether to return the mapping
#         tol (:py:class:`float`): tolerance for frequency comparison

#     Returns:
#         Omegas (:py:class:`numpy.ndarray`): array of generated frequencies
#         mapping (:py:class:`dict`, optional): Omega -> list of contributions (n, m)
#     """
    
#     freq_map = {}

#     def add_freq(n):
#         Omega = abs(n * omega)
#         if Omega < tol:
#             return   
#         key = np.round(Omega,12)
#         if key not in freq_map:
#             freq_map[key] = {"Omega": Omega,"contributions": []}
#         freq_map[key]["contributions"].append({"n": n})

#     all_harmonics = set(range(1, n_harmonics + 1))    
#     inactive_harmonics = set(inactive_harmonics or [])
#     active_harmonics = sorted(all_harmonics - inactive_harmonics)

#     for n in active_harmonics:
#         add_freq(n)
    
#     Omegas = np.array(sorted(freq_map.keys()))

#     if return_mapping:
#         return Omegas, freq_map
#     else:
#         return Omegas
    
def generate_frequencies(omega,n_harmonics=3,inactive_harmonics=None,tol=1e-8):
    """
    Generate positive frequencies for harmonic expansion with a single monochromatic field.

    Args:
        omega (:py:class:`float`): frequency of the first field (pump)
        n_harmonics (:py:class:`int`): number of harmonics to include in the expansion. Default is 3
        inactive_harmonics (:py:class:`list` or None): list with the harmonics to exclude from the expansion 
            Example: [2] excludes the 2nd harmonic. Default is None
        return_mapping (:py:class:`bool`): whether to return the mapping
        tol (:py:class:`float`): tolerance for frequency comparison

    Returns:
        Omegas_dict: {n: Omega_n = n * omega}
    """

    Omegas_dict = {}

    all_harmonics = set(range(1, n_harmonics + 1))
    inactive_harmonics = set(inactive_harmonics or [])
    active_harmonics = sorted(all_harmonics - inactive_harmonics)

    for n in active_harmonics:
        Omega = n * omega

        if abs(Omega) < tol:
            continue

        Omegas_dict[n] = Omega

    return Omegas_dict

# def fit_multiple_harmonics(t, y, omega, n_harmonics=2, inactive_harmonics=None, rcond=None):
#     """
#     Fit the data to:
#         f(t) = B0 + sum_{n=1}^{n_harmonics} A_n*sin(n*omega*t+phi_n)
    
#     at a given frequency omega, allowing exclusion of some selected harmonics.    
    
#     Args:
#         t (:py:class:`numpy.ndarray`): array with the time values
#         y (:py:class:`numpy.ndarray`): array with the values of the function
#         omega (:py:class:`float`): angular frequency of the sine functions
#         n_harmonics (:py:class:`int`): number of harmonics to fit. Default is 2
#         inactive_harmonics (:py:class:`list` or None): harmonics to exclude.
#                                            Example: [2] excludes the 2nd harmonic.
#         rcond (:py:class:`float`): lstsq parameter

#     Returns:
#         :py:class:`tuple` : tuple with the output, namely
#             - A (:py:class:`numpy.ndarray`): amplitudes
#             - phi (:py:class:`numpy.ndarray`): phases
#             - B0 (:py:class:`float`): constant offset
#             - residuals (:py:class:`numpy.ndarray`): euclidean norms of the squared errors of the fit

#         Harmonics that were excluded have A=0 and phi=0.
#     """
    
#     all_harmonics = set(range(1, n_harmonics + 1))    
    
#     inactive_harmonics = set(inactive_harmonics or [])
#     active_harmonics = sorted(all_harmonics - inactive_harmonics)
    
#     X = np.column_stack(
#         [np.ones_like(t)] + 
#         [np.sin(n * omega * t) for n in active_harmonics] + 
#         [np.cos(n * omega * t) for n in active_harmonics]
#     )
    
#     coeffs, residuals, _, _ = np.linalg.lstsq(X, y, rcond=rcond)
    
#     B0 = coeffs[0]
#     A = np.zeros(n_harmonics)
#     phi = np.zeros(n_harmonics)
    
#     n_active = len(active_harmonics)
    
#     for i, n in enumerate(active_harmonics):
#         alpha = coeffs[1 + i]
#         beta  = coeffs[1 + i + n_active]
        
#         A[n-1] = np.sqrt(alpha**2 + beta**2)
#         phi[n-1] = np.arctan2(beta, alpha)
    
#     return A, phi, B0, np.sqrt(residuals)

# def eval_multiple_harmonics(t, A, phi, B0, omega):
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


class Xn_single_frequency():
    """
    Class to extract the non-linear susceptibility from the polarization induced by a sine-shaped external field.

    Args:
        data (:py:class:`YamboNLDBParser`) : data parsed from the nlndb.Nonlineardatabase database
        X_order (:py:class:`int`) : order of the non-linear susceptibility to be extracted. The zero-th term corresponds to the constant offset of the polarization harmonic expansion. 
            Default is 3, which means that the non-linear susceptibility up to the third order is extracted
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

    def __init__(self,data,X_order=3,Trange=[-1, -1],Trange_units='fs',tol=1e-10,inactive_harmonics=None,verbose=True):
        self.time = data.IO_TIME_points
        self.pol = np.array(data.Polarization) 
        self.fields = data.Efield
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
        if len(data.Efield_general)>1 and data.Efield_general[1]["name"] != 'none':
                raise ValueError("This analysis is for one monochromatic field only.")
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
        print('Type of the external field:',self.fields[0]['name'])
        print('Order of the non-linear susceptibility extracted:',self.X_order)
        print('Number of frequencies:',self.nfreqs)
        print('Frequency range of the field (in Hartree):',self.fields_freqs[0],'-',self.fields_freqs[-1])
        print('Frequency range of the field (in eV):',self.fields_freqs[0]*C.HaToeV,'-',self.fields_freqs[-1]*C.HaToeV)
        print('Damping factor:',self.damp,'Hartree',self.damp*C.HaToeV,'eV')
        print('Dephasing time:',self.deph,'au',self.deph/C.FsToAu,'fs')
    
    def set_time_sampling(self,ifreq):
        """
        Define the time interval in which the polarization is sampled to compute the harmonic fit.
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
    
    def perform_harmonic_analysis(self,ifreq):
        """
        Perform harmonic analysis of the polarization using the fit_sum_frequencies function.

        Args:
            ifreq (:py:class:`int`) : index of the frequency of the external field
        
        Returns:
            :py:class:`tuple` : tuple with the values of the amplitude A and the phase phi of the fitted sine function.
                Each component of the tuple is an array with 3 values corresponding to the 3 cartesian directions.
        """
        X_order = self.X_order
        inactive_harmonics = self.inactive_harmonics
        iTstart, iTend = self.set_time_sampling(ifreq)
        t = self.time[iTstart:iTend]
        omega = self.fields_freqs[ifreq]

        A, phi, B0, residuals = np.zeros((X_order,3)), np.zeros((X_order,3)), np.zeros(3),  np.zeros(3)
        results_xyz = []
        Omegas = generate_frequencies(omega,n_harmonics=X_order,inactive_harmonics=inactive_harmonics)
        for idir in range(3):
            y = self.pol[ifreq,idir,iTstart:iTend]
            #A[:,idir], phi[:,idir], B0[idir], residuals[idir] = fit_sum_frequencies(t,y,Omegas,rcond=self.tol)
            results_xyz.append(fit_sum_frequencies(t,y,Omegas,rcond=self.tol))
            #A[:,idir], phi[:,idir], B0[idir], residuals[idir] = fit_multiple_harmonics(t,y,omega,n_harmonics=X_order,inactive_harmonics=inactive_harmonics,rcond=self.tol)
        #return A, phi, B0, residuals
        return results_xyz
    
    def check_harmonic_reliability(self,plot_ifreq=None,plot_dir=0):
        """
        Check the reliability of multiple harmonics fit by comparing the residuals with the amplitude of the fitted sine function. 
        If the residuals are larger than one tenth of the amplitude a warning is raised since the fit can be not accurate.
        
        Args:
            plot_ifreq (:py:class:`int` or None): index of the frequency of the external field to be plotted. If None, no plot is produced. Default is None
            plot_dir (:py:class:`int`): index of the cartesian direction to be plotted. Default is 0, which corresponds to the x direction

        """
        for ifreq in range(self.nfreqs):
            A, _, _, residuals = self.perform_harmonic_analysis(ifreq)
            for idir in range(3):
                if residuals[idir] > max(np.abs((A[:,idir])))/10.0:
                    print(f'Warning: Fit for frequency {ifreq} in direction {idir} is not accurate since the residuals are larger than the amplitude.')
                    print(f'Amplitudes (for each harmonic): {A[:,idir]}',f'Residuals of the fit: {residuals[idir]}')
        
        if plot_ifreq is not None:
            time = self.time
            pol = self.pol[plot_ifreq,plot_dir]
            iTstart, iTend = self.set_time_sampling(plot_ifreq)
            A, phi, B0, _ = self.perform_harmonic_analysis(plot_ifreq)
            Omegas = generate_frequencies(self.fields_freqs[plot_ifreq],n_harmonics=self.X_order,inactive_harmonics=self.inactive_harmonics)
            pol_fit = eval_sum_frequencies(time, A[:,plot_dir], phi[:,plot_dir], B0[plot_dir], Omegas)
            #pol_fit = eval_multiple_harmonics(time, A[:,plot_dir], phi[:,plot_dir], B0[plot_dir], self.fields_freqs[plot_ifreq])
            Tperiod = 2.0*np.pi/self.fields_freqs[plot_ifreq]
            Tmin,Tmax = np.max([time[iTstart]-5*Tperiod,time[0]])/C.FsToAu, np.min([time[iTend]+5*Tperiod,time[-1]])/C.FsToAu
            U.Plot_Array(time/C.FsToAu, pol_fit,xlim=(Tmin,Tmax), label='Harm fit',data2=pol,label2='Pol',figsize=(6,3))
            diffe = pol-pol_fit
            U.Plot_Array(time/C.FsToAu, diffe, label='Difference fit-pol',figsize=(6,3)) 

    def eval_Pw(self,plot=False):
        """
        Compute the polarization in the frequency domain at the (multiples) of the frequencies of the external fields according to the value of X_order. 
        Pw has ``X_order+1`` components, where the zero-th term corresponds to the constant offset of the polarization harmonic expansion. The first one is the
        response at the frequency of the external field, the second one is the response at the second harmonic and so on.
        
        Args:
            plot (:py:class:`bool`, optional): If True, plot the polarization in the frequency domain. Default is False

        Returns:
            :py:class:`numpy.ndarray` : array with the polarization in the frequency domain at the (multiple) frequencies of the external fields according 
                to the value of X_order. The shape of the array is (X_order,nfreqs,3) where 3 is the number of cartesian directions 
        """
        X_order = self.X_order
        pol_w = np.zeros((X_order+1,self.nfreqs,3),dtype=complex)
        for ifreq in range(self.nfreqs):
            A, phi, B0, _ = self.perform_harmonic_analysis(ifreq)
            pol_w[0,ifreq,:] = B0
            for n_harm in range(X_order):
                pol_w[n_harm+1,ifreq,:] = 1j*A[n_harm,:]/2.0*np.exp(-1j*phi[n_harm,:])
        
        if plot:
            energy = self.fields_freqs * C.HaToeV
            for n_harm in range(self.X_order+1):
                U.Plot_ComplexArray(energy, pol_w[n_harm], label=f'P_{n_harm}')
        
        return pol_w
    
    def eval_Ew_power(self):
        """
        Evaluate the n_harmonic power of the external fields in the frequency domain for all the values of the self.fields_freqs array. 
        The choice of the field factor for the zero-th order is done in agreement with the one of the YamboPy implementation of the non-linear susceptibility.
        
        Returns:
            :py:class:`numpy.ndarray` : array with the harmonic powers of the external fields in the frequency domain for all the values of the self.fields_freqs array.
                The shape of the array is (X_order+1,nfreqs)
                  
        """
        X_order = self.X_order
        efield_w_power = np.zeros((X_order+1,self.nfreqs),dtype=complex)
        for n_harm in range(X_order+1):
            for ifreq in range(self.nfreqs):
                E0 = self.fields[ifreq]['amplitude']
                t0 = self.fields[ifreq]['initial_time']
                omega = self.fields_freqs[ifreq]
                if n_harm == 0: 
                    efield_w_power[n_harm,ifreq] = -E0**2/4.0#*np.exp(1j*t0*omega)
                else:
                    efield_w_power[n_harm,ifreq] =np.power(1j*E0/2.0*np.exp(1j*t0*omega), n_harm)
        return efield_w_power
    
    def compute_Xn(self,set_units_of_measure=False,plot=False):
        """
        Compute the non-linear susceptibility of the system at the frequencies of the external fields. Analougly to the polarization, the non-linear susceptibility has ``X_order+1`` components,
        where the zero-th term corresponds to the constant offset of the polarization harmonic expansion. The first one is the response at the frequency of the external field, the second one is the response at the second harmonic and so on. 
        The non-linear susceptibility is computed as Xn = P(n*omega)/E^n. 
        
        Args:
            set_units_of_measure (:py:class:`bool`, optional): If True, set the units of measure of the non-linear susceptibility. Default is False. 
                The conversion factor is set in agreement with the YamboPy implementation of the non-linear susceptibility
            plot (:py:class:`bool`, optional): If True, plot the Xn in the frequency domain. Default is False

        Returns:
            :py:class:`numpy.ndarray` : array with the non-linear susceptibility of the system at the frequencies of the external fields. 
                The shape of the array is (X_order,nfreqs,3) where 3 is the number of cartesian directions
        """
        X_order = self.X_order
        pol_w = self.eval_Pw()
        Xn = np.zeros((X_order+1,self.nfreqs,3),dtype=complex)
        efield_w_power = self.eval_Ew_power()
        
        for n_harm in range(X_order+1):
            for idir in range(3):
                Xn[n_harm,:,idir] = pol_w[n_harm,:,idir] / efield_w_power[n_harm]
            if set_units_of_measure:
                Xn[n_harm] *= self.set_units_of_measure(n_harm)

        if plot:
            energy = self.fields_freqs * C.HaToeV
            for n_harm in range(self.X_order+1):
                U.Plot_ComplexArray(energy, Xn[n_harm], label=f'X_{n_harm}')
        
        return Xn
    
    def set_units_of_measure(self,n_harm):
        """
        Set the units of measure of the non-linear susceptibility in agreement with the YamboPy implementation. 

        Args:
            n_harm (:py:class:`int`) : order of the non-linear susceptibility.
        """
        ratio = C.SVCMm1_to_VMm1 / C.AU_toVMm1
        if n_harm == 0: 
            return np.power(ratio,1,dtype=np.float64)
        else:
            return np.power(ratio,n_harm-1,dtype=np.float64)

