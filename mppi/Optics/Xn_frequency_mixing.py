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

def generate_frequencies(omega1, omega2,max_order_E1=1,max_order_E2=3,include_pure_E1=True,include_pure_E2=True,include_mixing=True,tol=1e-8):
    """
    Generate positive frequencies as a dictionary:
        (n, m) -> Omega = |n*omega1 + m*omega2|

    Only positive frequencies are kept. The rationale behind the choice of producing positive frequencies only is that the
    only the positive frequency of the probe field is considered in the expansion of the polarization, so if for instance omega1-m*omega2 is negative, 
    it means that there is a positve contribution at m*omega2-omega1 which is not explicitly produced in the combination of the coefficients but is selected
    by the check on the positive frequencies performed here.

    Args:
        omega1 (:py:class:`float`): frequency of the second field (probe)
        omega2 (:py:class:`float`): frequency of the first field (pump)
        max_order_E1 (:py:class:`int`): maximum order for the probe. Default is 1 
        max_order_E2 (:py:class:`int`): maximum order for the pump. Default is 3
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
        for n in range(1, max_order_E1 + 1):
            add_freq(n, 0)
    if include_pure_E2:
        for m in range(1, max_order_E2 + 1):
            add_freq(0, m)
    if include_mixing:
        for n in range(1, max_order_E1 + 1):
            for m in range(1, max_order_E2 + 1):
                add_freq(n, m)
                add_freq(n, -m)

    return Omegas_dict

import numpy as np

def estimate_time_window(Omegas, safety_factor=8):
    """
    Estimate optimal time window for resolving frequencies. Assume that Omegas are expressed in Hartree.

    Args:
        Omegas (:py:class:`dict`): dictionary with keys (n, m) and values Omega
        safety_factor (:py:class:`float`): how conservative you want to be

    Returns:
        Time_window (:py:class:`float`): estimated optimal time window for resolving frequencies in au
    """
    Omegas_values= np.array([Omegas[key] for key in Omegas])
    
    diffs = []
    for i in range(len(Omegas_values)):
        for j in range(i+1, len(Omegas_values)):
            diffs.append(abs(Omegas_values[i] - Omegas_values[j]))
    
    delta_min = np.min(diffs)
    
    Time_window = safety_factor * 2.0*np.pi / delta_min
    
    return Time_window


class Xn_frequency_mixing():
    """
    Class to extract the non-linear susceptibility from the polarization induced by the sum of two monocromatic external fields, that 
    represent the pump and the probe in a typical pump-probe experiment. 

    Args:
        data (:py:class:`YamboNLDBParser`) : data parsed from the nlndb.Nonlineardatabase database
        X_order (:py:class:`inttuple`) : tuple with the (probe,pump) orders of the non-linear susceptibility to be extracted. 
            Default is (1, 3)
        Trange (:py:class:`list`) : list with the time range in which the polarization is sampled to compute the sine fit 
            Default is [-1,-1], which means that a time sampling of one period of the external field is used
            (Tstart = time[-1] - Tperiod, Tend = time[-1])
            Further details are provided in the documentation of the method :py:meth:`set_time_sampling`
        Trange_units (:py:class:`str`) : set the units to time sampling. Default is 'fs' and the other possible choice is 'au'
        tol (:py:class:`float`) : tolerance for the fit of the sine function. Default is 1e-10
        ref_pol (:py:class:`numpy.ndarray` or None) : array with the reference polarization to be subtracted from the original one. Default is None
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal
    
    Attributes:
        time (:py:class:`numpy.ndarray`) : array with the time values in au
        pol (:py:class:`numpy.ndarray`) : array with the polarization of the system at each time step in au
        probes (:py:class:`list`) : list of dictionaries with the information on the probe field for each frequency
        pump (:py:class:`dict`) : dictionary with the information on the pump field
        nfreqs (:py:class:`int`) : number of frequencies of the probe field
        probe_freqs (:py:class:`numpy.ndarray`) : array with the frequencies of the probes in Hartree
        pump_freq (:py:class:`float`) : frequency of the pump in Hartree
        damp (:py:class:`float`) : damping factor in Hartree
        deph (:py:class:`float`) : dephasing time in au, defined as 12/damp. If the time sampling interval starts before this value
            a warning is raised since the fit of the sine function can be not accurate.
        dt (:py:class:`float`) : time sampling interval in au
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal
    """

    def __init__(self,data,X_order=(1, 3),Trange=[-1, -1],Trange_units='fs',tol=1e-10,ref_pol=None,verbose=True):
        self.time = data.IO_TIME_points
        self.pol = np.array(data.Polarization)
        if ref_pol is not None:
            self.pol -= ref_pol
            print('Reference polarization found. The difference between the original polarization and the reference one is computed.')
        self.probes = data.Efield
        self.pump = data.Efield_general[1]
        self.nfreqs = data.n_frequencies
        self.probe_freqs = np.array([e['freq_range'][0] for e in self.probes]) 
        self.pump_freq = self.pump['freq_range'][0]
        self.damp = data.NL_damping
        self.deph = 6/self.damp
        self.X_order = X_order
        self.dt = self.time[1]-self.time[0]
        self.Trange = Trange
        self.Trange_units = Trange_units
        self.tol = tol
        self.verbose = verbose
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
        print('Order of the non-linear susceptibility extracted:',self.X_order,'(probe,pump)')
        print('Number of frequencies:',self.nfreqs)
        print('Frequency range of the probe (in Hartree):',self.probe_freqs[0],'-',self.probe_freqs[-1])
        print('Frequency range of the probe (in eV):',self.probe_freqs[0]*C.HaToeV,'-',self.probe_freqs[-1]*C.HaToeV)
        print('Frequency of the pump :',self.pump_freq,' Hartree -',self.pump_freq*C.HaToeV,' eV')
        print('Damping factor:',self.damp,'Hartree',self.damp*C.HaToeV,'eV')
        print('Dephasing time:',self.deph,'au',self.deph/C.FsToAu,'fs')
    
    def set_time_sampling(self,ifreq):
        """
        Define the time interval in which the polarization is sampled to compute the harmonic fit.
        For frequency mixing analysis, the polarization is not periodic so the default time range T is chosen in order to be able
        to resolve the frequency of the harmonic of the fields according to Dw = 2*pi/T, where Dw is the frequency resolution and is
        computed as the minimum value of the frequencies of the harmonics of the fields.
        So if Trange=[-1, -1] (T_start = time[-1] - T, T_end = time[-1]).
        Instead if positive values of Trange are provided the sampling time range is set to that values. Checks are performed to 
        ensure that the time range is consistent with the time sampling of the polarization. Lastly, if the time sampling interval 
        starts before the dephasing time a warning is raised since the fit of the sine function can be not accurate.
        
        Args:
            ifreq (:py:class:`int`) : index of the frequency of the external field
        
        Returns:
            :py:class:`tuple` : tuple with the indexes of the start and end time of the sampling interval
        """

        Omegas = generate_frequencies(self.probe_freqs[ifreq], self.pump_freq, max_order_E1=self.X_order[0], max_order_E2=self.X_order[1])
        Time_window = estimate_time_window(Omegas)
        sim_time = self.time[-1] - self.time[0]
        if self.verbose:
            print(f'Estimated optimal time window for resolving frequency {ifreq}: {Time_window:.2f} au - {Time_window/C.FsToAu:.2f} fs')
        
        if self.Trange_units == 'fs':
            time = self.time/C.FsToAu
            dt = self.dt/C.FsToAu
            deph = self.deph/C.FsToAu
            Time_window = Time_window/C.FsToAu
            sim_time = sim_time/C.FsToAu
        elif self.Trange_units == 'au':
            time = self.time
            dt = self.dt
            deph = self.deph
        else:
            raise ValueError("Invalid time units. Please use 'fs' or 'au'.")
                
        if self.Trange[0] < 0.:
            if Time_window < sim_time:   
                Tstart = time[-1] - Time_window
            else:
                Tstart = deph
                print(f'Warning: the optimal time window for resolving the frequencies {ifreq} is larger than the simulation time')
        else:
            Tstart = self.Trange[0]
        
        if Tstart < deph:
            print('Warning: the time sampling starts before the dephasing time. The fit can be not accurate.') 
        if Tstart > time[-1]:
            raise ValueError('The time sampling starts after the end of the time range') 
        
        if self.Trange[1] < 0.:   
            Tend = time[-1]   
        else:            
            Tend = self.Trange[1]
        if Tend > time[-1]:
            raise ValueError('The time sampling ends after the end of the time range') 
        if Tend < (Tstart + Time_window):
            print('Warning: the time sampling is shorter than the estimated optimal time window.')
        
        iTstart = int(np.round( Tstart / dt))
        iTend = int(np.round( Tend / dt)) 
        if self.verbose:
            print(f'Time sampling for frequency {ifreq}: {Tstart:.2f} - {Tend:.2f} {self.Trange_units} - indexes: {iTstart} - {iTend}')
        return iTstart, iTend
    
    def perform_harmonic_analysis(self,ifreq):
        """
        Perform harmonic analysis of the polarization using the fit_sum_frequencies function.

        Args:
            ifreq (:py:class:`int`) : index of the frequency of the external field
        
        Returns:
            :py:class:`list` : list with the results for each  cartesian direction.
                Each element of the list is a tuple with A,phi (organized in a dict in which the keys are the harmonic indices), B0 and
                residuals
        """
        X_order = self.X_order
        iTstart, iTend = self.set_time_sampling(ifreq)
        t = self.time[iTstart:iTend]
        omega1 = self.probe_freqs[ifreq]
        omega2 = self.pump_freq
        
        results_xyz = []
        Omegas = generate_frequencies(omega1, omega2,max_order_E1=X_order[0],max_order_E2=X_order[1])
        for idir in range(3):
            y = self.pol[ifreq,idir,iTstart:iTend]
            results_xyz.append(fit_sum_frequencies(t,y,Omegas,rcond=self.tol))
        
        return results_xyz
    
    def check_harmonic_reliability(self,plot_ifreq=None,plot_dir=0):
        """
        Check the reliability of multiple harmonics fit by comparing the residuals with the amplitude of the fitted sine function. 
        If the residuals are larger than one tenth of the amplitude a warning is raised since the fit can be not accurate.
        if plot_ifreq is not None, the fit of the sine function is plotted together with the polarization for the specified frequency and cartesian direction.
        
        Args:
            plot_ifreq (:py:class:`int` or None): index of the frequency of the external field to be plotted. If None, no plot is produced. Default is None
            plot_dir (:py:class:`int`): index of the cartesian direction to be plotted. Default is 0, which corresponds to the x direction

        """
        for ifreq in range(self.nfreqs):
            results_xyz = self.perform_harmonic_analysis(ifreq)
            for idir in range(3):
                A = np.array([res["A"] for res in results_xyz[idir][0].values()])
                A_max = max(A)
                residuals = results_xyz[idir][2]
                ratio = residuals / A_max if A_max > 0. else 0.
                if ratio > 0.1:
                    print(f'Warning: Fit for frequency {ifreq} in direction {idir} is not accurate since the residuals are comparable to the amplitude.')
                    print(f'Amplitudes (for each harmonic): {A}',f'Residuals of the fit: {residuals}')
        
        if plot_ifreq is not None:
            time = self.time
            pol = self.pol[plot_ifreq,plot_dir]
            iTstart, _ = self.set_time_sampling(plot_ifreq)
            print(iTstart)
            Omegas = generate_frequencies(self.probe_freqs[ifreq], self.pump_freq, max_order_E1=self.X_order[0], max_order_E2=self.X_order[1])
            omega_min = min(Omegas.values())
            Tperiod_max = 2.0*np.pi/omega_min
            results = self.perform_harmonic_analysis(plot_ifreq)[plot_dir]
            pol_fit = eval_sum_frequencies(time, results[0], results[1])
            Tmin,Tmax = time[iTstart]/C.FsToAu, np.min([time[iTstart]+5*Tperiod_max,time[-1]])/C.FsToAu
            print('Tmin,Tmax for the plot:',iTstart,Tmin,Tmax)
            U.Plot_Array(time/C.FsToAu, pol_fit,xlim=(Tmin,Tmax), label='Harm fit',data2=pol,label2='Pol',figsize=(6,3))
            diffe = pol-pol_fit
            U.Plot_Array(time/C.FsToAu, diffe, label='Difference fit-pol',figsize=(6,3)) 
    
    def eval_Pw(self,plot=False,plot_dir=0):
        """
        Compute the polarization in the frequency domain at the (multiples) of the frequencies of the external fields according to the value of the 
        'X_order' and 'inactive_harmonics' variables. Pw at the harmonic omega is computed as 
            P(omega) = 1j*A(omega)/2.0*exp(-1j*phi(omega))

        For each cartesian direction Pw is a dict with the keys of the generate_frequencies function.  
        
        Args:
            plot (:py:class:`bool`, optional): If True, plot the polarization in the frequency domain. Default is False
            plot_dir (:py:class:`int`, optional): Index of the cartesian direction to be plotted. Default is 0, which corresponds to the x direction

        Returns:
            :py:class:`numpy.list` : each element contains a cartesian direction. For each direction the function returns a dict
                with the keys of the generate_frequencies function and the values of the polarization in the frequency domain 
        """

        Pw_xyz = []

        for idir in range(3):       
            Pw = {}
            Pw[(0, 0)] = np.zeros(self.nfreqs, dtype=complex) 
            for ifreq in range(self.nfreqs):
                results = self.perform_harmonic_analysis(ifreq)[idir]
                Pw[(0, 0)][ifreq] = results[1] # constant offset
                for key, res in results[0].items():
                    if key not in Pw:
                        Pw[key] = np.zeros(self.nfreqs, dtype=complex)
                    Pw[key][ifreq] = 1j*res["A"] / 2.0 * np.exp(-1j * res["phi"])
            Pw_xyz.append(Pw)
        
        if plot:
            dir = ['x','y','z']
            energy = self.probe_freqs * C.HaToeV
            harmonics = list(Pw_xyz[plot_dir].keys())
            for n_harm in harmonics:
                U.Plot_ComplexArray(energy, Pw_xyz[plot_dir][n_harm], label=f'P_{dir[plot_dir]}[{n_harm}]')
        
        return Pw_xyz
    
    def eval_Ew(self):
        """
        Evaluate the product of the (n_harmonic powers of) the pump and probe fields in the frequency domain for all the values of the self.fields_freqs array. 
        The choice of the field factors is done in agreement with the one of the YamboPy implementation for the frequency mixing of the non-linear susceptibility.
        
        Returns:
            :py:class:`numpy.dict` : dict with the harmonic powers of the external fields in the frequency domain for all the values of the self.fields_freqs array
                  
        """
        Ew = {}
        Pw_x = self.eval_Pw()[0] 
        E0_pump = self.pump['amplitude']
        t0_pump = self.pump['initial_time']
        omega_pump = self.pump_freq
        for ifreq in range(self.nfreqs):
            E0_probe = self.probes[ifreq]['amplitude']
            t0_probe = self.probes[ifreq]['initial_time']
            omega_probe = self.probe_freqs[ifreq]

            for key in Pw_x:
                if key not in Ew:   
                    Ew[key] = np.zeros(self.nfreqs, dtype=complex)
                n_p,n_P = np.abs(key)
                if n_p == 0:
                    Ew[key][ifreq] = 1.0
                else:
                    Ew[key][ifreq] = np.power(1j * E0_probe / 2.0 * np.exp(1j * t0_probe * omega_probe),n_p)
                if n_P == 0:
                    Ew[key][ifreq] *= 1.0
                else:
                    Ew[key][ifreq] *= np.power(1j * E0_pump / 2.0 * np.exp(1j * t0_pump * omega_pump),n_P)

        return Ew
    
    def compute_Xn(self,set_units_of_measure=False,plot=False,plot_dir=0):
        """
        Compute the non-linear susceptibility of the system at (multiple of the) frequencies of the external fields. 
        For each cartesian direction Xn is a dict with the keys of the generate_frequencies function.
        
        Args:
            set_units_of_measure (:py:class:`bool`, optional): If True, set the units of measure of the non-linear susceptibility. Default is False. 
                The conversion factor is set in agreement with the YamboPy implementation of the non-linear susceptibility
            plot (:py:class:`bool`, optional): If True, plot the Xn in the frequency domain. Default is False
            plot_dir (:py:class:`int`, optional): Index of the cartesian direction to be plotted. Default is 0, which corresponds to the x direction

        Returns:
            :py:class:`numpy.list` : each element contains a cartesian direction. For each direction the function returns a dict
                with the keys of the generate_frequencies function and the values of the Xn in the frequency domain 
        """

        Xn_xyz = []
        Pw_xyz = self.eval_Pw()
        Ew = self.eval_Ew()
        for idir in range(3):   
            Xn = {}
            for key, res in Pw_xyz[idir].items():
                n_p,n_P = np.abs(key)
                if key not in Xn:
                    Xn[key] = np.zeros(self.nfreqs, dtype=complex)
                Xn[key] = Pw_xyz[idir][key] / Ew[key]
                if set_units_of_measure:
                    Xn[key] *= self.set_units_of_measure(n_p+n_P)
            Xn_xyz.append(Xn)

        if plot:
            dir = ['x','y','z']
            energy = self.probe_freqs * C.HaToeV
            harmonics = list(Pw_xyz[plot_dir].keys())
            for n_harm in harmonics:
                U.Plot_ComplexArray(energy, Xn_xyz[plot_dir][n_harm], label=f'X_{dir[plot_dir]}[{n_harm}]')
        
        return Xn_xyz
    
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


