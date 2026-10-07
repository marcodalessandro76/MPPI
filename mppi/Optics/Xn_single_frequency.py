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
from mppi.Optics.Utils import fit_sum_frequencies, eval_sum_frequencies, print_sampling_warnings, lorentzian_broadening
from mppi.Parsers import YamboNLDBParser

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


class Xn_single_frequency():
    r"""
    Class to extract the non-linear susceptibilities from the polarization induced by a sine-shaped external field.

    Convention for the susceptibilities. Fields and polarization are expanded as :math:`E(t) = \sum_n E(\omega_n)e^{-i\omega_n t}`, with
    :math:`E(-\omega) = E(\omega)^*`, as in R. W. Boyd, *Nonlinear Optics*, 4th ed. (2020), Section 1.3. The susceptibility of the n-th harmonic
    is the ratio between the component of the polarization at :math:`n\omega` and :math:`E(\omega)^n`, i.e. the susceptibility
    :math:`\chi^{(n)}(n\omega;\omega,\dots,\omega)` of Boyd (degeneracy factor D=1). The zero-th order is the ratio between the static
    polarization and :math:`|E(\omega)|^2`, so it includes the degeneracy factor D=2 of Boyd, Eq. (1.3.19):
    :math:`\chi_0 = 2\chi^{(2)}(0;\omega,-\omega)`. Note that the YamboPy implementation divides the zero-th order by
    :math:`E(\omega)^2 = -|E(\omega)|^2` instead, so its value has the opposite sign.

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
        self.verbose = verbose
        self._fits = {} # results of the harmonic analysis, see _harmonic_analysis
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
        data = YamboNLDBParser(file,verbose=verbose)
        return cls(data,verbose=verbose)
    
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
        The default time range is the last period of the external field (T_end = time[-1], T_start = T_end - T_period).
        Instead if positive values of Trange are provided the sampling time range is set to that values. Checks are performed to
        ensure that the time range is consistent with the time sampling of the polarization (with a tolerance of half a time
        step) and that it contains at least one period of the field. Lastly, if the time sampling interval starts before the
        dephasing time a warning is raised since the fit of the sine function can be not accurate.
        
        Args:
            ifreq (:py:class:`int`) : index of the frequency of the external field

        Returns:
            :py:class:`tuple` : tuple with the indexes of the start and end time of the sampling interval
        """
        iTstart, iTend, warnings = self._time_sampling(ifreq)
        print_sampling_warnings({ifreq: warnings}, self.nfreqs)
        return iTstart, iTend

    def _time_sampling(self,ifreq):
        """
        Compute the time sampling interval of :py:meth:`set_time_sampling`.

        Returns:
            :py:class:`tuple` : the indexes of the start and end time of the sampling interval and the list of the warnings
        """
        warnings = []
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
        tol = dt/2. # tolerance on the comparison between times

        if self.Trange[1] < 0.:
            Tend = time[-1]
        else:
            Tend = self.Trange[1]
        if Tend > time[-1] + tol:
            Tend = time[-1]
            warnings.append('the time sampling ends after the end of the time range. Tend is set to time[-1]')

        if self.Trange[0] < 0.:
            Tstart = Tend - Tperiod
        else:
            Tstart = self.Trange[0]
        if Tstart > time[-1] - tol:
            Tstart = time[-1] - Tperiod
            warnings.append('the time sampling starts after the end of the time range. Tstart is set to time[-1]-Tperiod')
        if Tstart < deph - tol:
            warnings.append('the time sampling starts before the dephasing time. The fit of the sine function can be not accurate')

        if Tend < Tstart + Tperiod - tol:
            Tend = Tstart + Tperiod
            warnings.append('the time sampling is shorter than one period of the field. Tend is set to Tstart + Tperiod')

        iTstart = int(np.round((Tstart - time[0]) / dt))
        iTend = int(np.round((Tend - time[0]) / dt))

        return iTstart, iTend, warnings

    def perform_harmonic_analysis(self,ifreq):
        """
        Perform harmonic analysis of the polarization using the fit_sum_frequencies function.
        The results are stored, so the fit of each frequency is performed only once (as long as X_order, Trange,
        Trange_units, tol and inactive_harmonics are not changed).

        Args:
            ifreq (:py:class:`int`) : index of the frequency of the external field

        Returns:
            :py:class:`list` : list with the results for each  cartesian direction.
                Each element of the list is a tuple with A,phi (organized in a dict in which the keys are the harmonic indices), B0 and
                residuals
        """
        results_xyz, warnings = self._harmonic_analysis(ifreq)
        print_sampling_warnings({ifreq: warnings}, self.nfreqs)
        return results_xyz

    def _harmonic_analysis(self,ifreq):
        """
        Perform (or read from the stored results) the harmonic analysis of :py:meth:`perform_harmonic_analysis`.

        Returns:
            :py:class:`tuple` : the results of the fit for each cartesian direction and the list of the warnings raised
                by the choice of the time sampling
        """
        inactive = tuple(sorted(self.inactive_harmonics)) if self.inactive_harmonics is not None else None
        key = (ifreq, self.X_order, tuple(self.Trange), self.Trange_units, self.tol, inactive)
        if key not in self._fits:
            iTstart, iTend, warnings = self._time_sampling(ifreq)
            t = self.time[iTstart:iTend]
            Omegas = generate_frequencies(self.fields_freqs[ifreq],n_harmonics=self.X_order,inactive_harmonics=self.inactive_harmonics)
            results_xyz = []
            for idir in range(3):
                y = self.pol[ifreq,idir,iTstart:iTend]
                results_xyz.append(fit_sum_frequencies(t,y,Omegas,rcond=self.tol))
            self._fits[key] = (results_xyz, warnings)
        return self._fits[key]

    def _print_sampling_warnings(self):
        """
        Print a summary of the warnings raised by the time sampling of all the frequencies.
        """
        warnings = {ifreq: self._harmonic_analysis(ifreq)[1] for ifreq in range(self.nfreqs)}
        print_sampling_warnings(warnings, self.nfreqs)

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
            results_xyz = self._harmonic_analysis(ifreq)[0]
            for idir in range(3):
                A = np.array([res["A"] for res in results_xyz[idir][0].values()])
                A_max = max(A)
                residuals = results_xyz[idir][2]
                ratio = residuals / A_max if A_max > 0. else 0.
                if ratio > 0.1:
                    print(f'Warning: Fit for frequency {ifreq} in direction {idir} is not accurate since the residuals are comparable to the amplitude.')
                    print(f'Amplitudes (for each harmonic): {A}',f'Residuals of the fit: {residuals}')
        self._print_sampling_warnings()

        if plot_ifreq is not None:
            time = self.time
            pol = self.pol[plot_ifreq,plot_dir]
            iTstart, iTend = self._time_sampling(plot_ifreq)[:2]
            results = self._harmonic_analysis(plot_ifreq)[0][plot_dir]
            pol_fit = eval_sum_frequencies(time, results[0], results[1])
            Tperiod = 2.0*np.pi/self.fields_freqs[plot_ifreq]
            Tmin,Tmax = np.max([time[iTstart]-5*Tperiod,time[0]])/C.FsToAu, np.min([time[iTend]+5*Tperiod,time[-1]])/C.FsToAu
            U.Plot_Array(time/C.FsToAu, pol_fit,xlim=(Tmin,Tmax), label='Harm fit',data2=pol,label2='Pol',figsize=(6,3))
            diffe = pol-pol_fit
            U.Plot_Array(time/C.FsToAu, diffe, label='Difference fit-pol',figsize=(6,3)) 

    def eval_Pw(self,plot=False,plot_dir=0):
        """
        Compute the polarization in the frequency domain at the (multiples) of the frequencies of the external fields according to the value of the 
        'X_order' and 'inactive_harmonics' variables. Pw at the harmonic omega is computed as::

            P(omega) = 1j*A(omega)/2.0*exp(-1j*phi(omega))

        For each cartesian direction Pw is a dict with the keys of the generate_frequencies function.  
        
        Args:
            plot (:py:class:`bool`, optional): If True, plot the polarization in the frequency domain. Default is False
            plot_dir (:py:class:`int`, optional): Index of the cartesian direction to be plotted. Default is 0, which corresponds to the x direction

        Returns:
            :py:class:`numpy.list` : each element contains a cartesian direction. For each direction the function returns a dict
                with the keys of the generate_frequencies function and the values of the polarization in the frequency domain 
        """

        Pw_xyz = self._eval_Pw()
        self._print_sampling_warnings()

        if plot:
            dir = ['x','y','z']
            energy = self.fields_freqs * C.HaToeV
            harmonics = list(Pw_xyz[plot_dir].keys())
            for n_harm in harmonics:
                U.Plot_ComplexArray(energy, Pw_xyz[plot_dir][n_harm], label=f'P_{dir[plot_dir]}[{n_harm}]')
        
        return Pw_xyz

    def _eval_Pw(self):
        """
        Compute the polarization in the frequency domain of :py:meth:`eval_Pw`, without printing the warnings.
        """
        Pw_xyz = []
        for idir in range(3):
            Pw = {}
            Pw[0] = np.zeros(self.nfreqs, dtype=complex)
            for ifreq in range(self.nfreqs):
                results = self._harmonic_analysis(ifreq)[0][idir]
                Pw[0][ifreq] = results[1] # constant offset
                for key, res in results[0].items():
                    if key not in Pw:
                        Pw[key] = np.zeros(self.nfreqs, dtype=complex)
                    Pw[key][ifreq] = 1j*res["A"] / 2.0 * np.exp(-1j * res["phi"])
            Pw_xyz.append(Pw)
        return Pw_xyz

    def eval_Ew(self):
        r"""
        Evaluate the (n_harmonic powers of) the external fields in the frequency domain for all the values of the self.fields_freqs array.
        The field :math:`E_0\sin(\omega(t-t_0))` has the complex amplitudes :math:`E(\omega) = (iE_0/2)e^{i\omega t_0}` and
        :math:`E(-\omega) = E(\omega)^*` (Boyd, Nonlinear Optics, Eq. (1.3.6)). The n-th harmonic is associated to :math:`E(\omega)^n`,
        while the zero-th order (optical rectification) is associated to :math:`E(\omega)E(-\omega) = |E(\omega)|^2 = E_0^2/4`.
        
        Returns:
            :py:class:`numpy.dict` : dict with the harmonic powers of the external fields in the frequency domain for all the values of the self.fields_freqs array
                  
        """
        return self._eval_Ew(self._eval_Pw()[0].keys())

    def _eval_Ew(self,keys):
        """
        Compute the powers of the field amplitudes of :py:meth:`eval_Ew` for the given keys.
        """
        Ew = {}
        for ifreq in range(self.nfreqs):
            E0 = self.fields[ifreq]['amplitude']
            t0 = self.fields[ifreq]['initial_time']
            omega = self.fields_freqs[ifreq]

            for key in keys:
                if key not in Ew:   
                    Ew[key] = np.zeros(self.nfreqs, dtype=complex)
                n = abs(key)
                if n == 0:
                    Ew[key][ifreq] = E0**2 / 4.0
                else:
                    Ew[key][ifreq] = np.power(1j * E0 / 2.0 * np.exp(1j * t0 * omega),n)

        return Ew

    def compute_Xn(self,set_units_of_measure=False,plot=False,plot_dir=0,broadening=None,renormalize=True):
        """
        Compute the non-linear susceptibility of the system at (multiple of the) frequencies of the external fields.
        For each cartesian direction Xn is a dict with the keys of the generate_frequencies function.

        The broadening of the linear susceptibility (key 1) can be increased a posteriori with a Lorentzian convolution on the
        grid of the frequencies of the fields, which gives :math:`\\chi^{(1)}(\\omega+i\\Delta)` (see
        :py:func:`mppi.Optics.Utils.lorentzian_broadening`). It requires a grid step much smaller than the broadening. The
        other keys are not broadened: for the harmonics n >= 2 the convolution would add n times the broadening to the width
        of the n-photon resonance, and it is not defined for the zero-th order.

        Args:
            set_units_of_measure (:py:class:`bool`, optional): If True, set the units of measure of the non-linear susceptibility. Default is False.
                The conversion factor is set in agreement with the YamboPy implementation of the non-linear susceptibility
            plot (:py:class:`bool`, optional): If True, plot the Xn in the frequency domain. Default is False
            plot_dir (:py:class:`int`, optional): Index of the cartesian direction to be plotted. Default is 0, which corresponds to the x direction
            broadening (:py:class:`float`, optional): Lorentzian broadening (eV) added to the linear susceptibility. Default is None (no broadening)
            renormalize (:py:class:`bool`, optional): normalize the Lorentzian kernel to one on the grid of the frequencies, see
                :py:func:`mppi.Optics.Utils.lorentzian_broadening`. Default is True

        Returns:
            :py:class:`numpy.list` : each element contains a cartesian direction. For each direction the function returns a dict
                with the keys of the generate_frequencies function and the values of the Xn in the frequency domain
        """

        Xn_xyz = []
        Pw_xyz = self._eval_Pw()
        Ew = self._eval_Ew(Pw_xyz[0].keys())
        self._print_sampling_warnings()
        for idir in range(3):
            Xn = {}
            for key, res in Pw_xyz[idir].items():
                n = abs(key)
                if key not in Xn:
                    Xn[key] = np.zeros(self.nfreqs, dtype=complex)
                Xn[key] = Pw_xyz[idir][key] / Ew[key]
                if set_units_of_measure:
                    Xn[key] *= self.set_units_of_measure(n)
                if broadening is not None and key == 1:
                    Xn[key] = lorentzian_broadening(self.fields_freqs, Xn[key], broadening/C.HaToeV, renormalize=renormalize)
            Xn_xyz.append(Xn)
        if broadening is not None and self.verbose:
            print(f'Lorentzian broadening of {broadening} eV applied to the key 1, the other keys are not broadened')

        if plot:
            dir = ['x','y','z']
            energy = self.fields_freqs * C.HaToeV
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

