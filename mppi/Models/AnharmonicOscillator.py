r"""
This module implements the classical anharmonic oscillator model of the non-linear optical response
described in R. W. Boyd, *Nonlinear Optics*, 4th edition (Academic Press, 2020), Section 1.4.

The model provides both the analytical expression of the linear and non-linear susceptibilities and the
numerical solution of the equation of motion in the time domain. The resulting polarization is stored in an
object with the same structure of the :py:class:`YamboNLDBParser` class, so it can be analyzed with the tools of the
Optics module (:py:class:`Xn_single_frequency`, :py:class:`Xn_frequency_mixing`, :py:func:`Linear_Response`)
and the results can be compared with the exact susceptibilities.

All the quantities are expressed in atomic units (:math:`e = m = \hbar = 1`, the charge of the electron is -1)
and the polarization is related to the field as :math:`P = \chi E` (gaussian-like convention, no :math:`\epsilon_0`,
as in the analysis of the yambo_nl results). The module can be loaded in the notebook as follows

>>> from mppi.Models import AnharmonicOscillator as AO

>>> osc = AO.AnharmonicOscillator(omega0=0.5,gamma=0.01,a=0.1)

"""
import numpy as np

class NLData():
    """
    Container of the results of a computation of the time dependent polarization. Its attributes have the
    same names of the ones of the :py:class:`YamboNLDBParser` class used by the Optics module.

    Attributes:
        IO_TIME_points (:py:class:`numpy.ndarray`) : time points in au
        Polarization (:py:class:`list`) : for each run, array with shape (3, time) with the polarization
        Efield (:py:class:`list`) : for each run, dictionary with the parameters of the first field
        Efield2 (:py:class:`list`) : for each run, dictionary with the parameters of the second field
        Efield_general (:py:class:`list`) : dictionaries with the parameters of the two external fields
        n_frequencies (:py:class:`int`) : number of runs (frequencies of the first field)
        n_runs (:py:class:`int`) : number of runs
        NL_damping (:py:class:`float`) : damping of the oscillator (in Hartree), that sets the decay rate of the
            transient of the polarization
    """
    def __init__(self,time,polarization,efields,efields2,damping):
        self.IO_TIME_points = time
        self.Polarization = polarization
        self.Efield = efields
        self.Efield2 = efields2
        self.Efield_general = [efields[0],efields2[0]]
        self.n_frequencies = len(efields)
        self.n_runs = len(efields)
        self.N_ext_fields = 1 if efields2[0]['name'] == 'none' else 2
        self.NL_damping = damping

    def get_time(self,convert_to_fs=True):
        """
        Return the time points (in fs or au)
        """
        from mppi.Utilities.Constants import FsToAu
        return self.IO_TIME_points/FsToAu if convert_to_fs else self.IO_TIME_points

def build_field(name='SIN',frequency=0.,amplitude=0.,initial_time=0.):
    """
    Build a dictionary with the parameters of an external field, with the keys used by the :py:class:`YamboNLDBParser` class.
    The field is directed along x.

    Args:
        name (:py:class:`str`) : type of the field, 'SIN', 'DELTA' or 'none'
        frequency (:py:class:`float`) : frequency of the field (Hartree)
        amplitude (:py:class:`float`) : amplitude of the field (au)
        initial_time (:py:class:`float`) : switch on time of the field (au)

    Returns:
        :py:class:`dict` : dictionary with the parameters of the field
    """
    return {'name':name,'versor':np.array([1.,0.,0.]),'amplitude':amplitude,'initial_time':initial_time,
            'freq_range':np.array([frequency,frequency]),'intensity':0.}

class AnharmonicOscillator():
    r"""
    Classical anharmonic oscillator (Boyd, Section 1.4). The displacement :math:`x` of the electron obeys

    .. math::
        \ddot{x} + 2\gamma\dot{x} + \omega_0^2 x + a x^2 - b x^3 = -E(t)

    and the polarization is :math:`P = -N x`. The term :math:`a` (non-centrosymmetric medium) produces the second order
    response, while the term :math:`b` (centrosymmetric medium) produces the third order one.

    The susceptibilities are defined according to Boyd, Eqs. (1.3.12) and (1.3.20): the fields are written as
    :math:`E(t) = \sum_n E(\omega_n) e^{-i\omega_n t}` with :math:`E(-\omega) = E(\omega)^*`, and the component of the
    polarization at :math:`\omega_\sigma = \omega_1 + \dots + \omega_n` is

    .. math::
        P(\omega_\sigma) = D\,\chi^{(n)}(\omega_\sigma;\omega_1,\dots,\omega_n)\,E(\omega_1)\cdots E(\omega_n)

    where the degeneracy factor :math:`D` is the number of distinct permutations of the frequencies of the fields.

    Args:
        omega0 (:py:class:`float`) : resonance frequency (Hartree)
        gamma (:py:class:`float`) : damping (Hartree). The full width at half maximum of the absorption is 2*gamma
        a (:py:class:`float`) : coefficient of the quadratic term of the restoring force
        b (:py:class:`float`) : coefficient of the cubic term of the restoring force
        density (:py:class:`float`) : number of oscillators per unit volume (au)
    """

    def __init__(self,omega0=0.5,gamma=0.01,a=0.,b=0.,density=1.):
        self.omega0 = omega0
        self.gamma = gamma
        self.a = a
        self.b = b
        self.density = density

    def D(self,omega):
        r"""
        Resonance denominator :math:`D(\omega) = \omega_0^2 - \omega^2 - 2i\omega\gamma` (Boyd, Eq. (1.4.10))
        """
        omega = np.asarray(omega)
        return self.omega0**2 - omega**2 - 2j*omega*self.gamma

    def chi1(self,omega):
        r"""
        Linear susceptibility :math:`\chi^{(1)}(\omega) = N/D(\omega)` (Boyd, Eq. (1.4.17a))
        """
        return self.density/self.D(omega)

    def chi2(self,omega1,omega2):
        r"""
        Second order susceptibility :math:`\chi^{(2)}(\omega_1+\omega_2;\omega_1,\omega_2)` (Boyd, Eqs. (1.4.20),
        (1.4.24), (1.4.26) and (1.4.27)). The frequencies are signed, for instance chi2(w1,-w2) is the difference
        frequency generation and chi2(w,-w) the optical rectification.
        """
        omega1, omega2 = np.asarray(omega1), np.asarray(omega2)
        return self.density*self.a/(self.D(omega1+omega2)*self.D(omega1)*self.D(omega2))

    def chi3(self,omega1,omega2,omega3):
        r"""
        Third order susceptibility :math:`\chi^{(3)}(\omega_1+\omega_2+\omega_3;\omega_1,\omega_2,\omega_3)` of the
        centrosymmetric oscillator (Boyd, Eq. (1.4.52) for the xxxx component). The frequencies are signed.
        The expression is valid for a = 0, otherwise the cascaded second order processes give a further third order
        contribution.
        """
        if self.a != 0.:
            raise ValueError('chi3 is implemented only for the centrosymmetric oscillator (a = 0)')
        omega1, omega2, omega3 = np.asarray(omega1), np.asarray(omega2), np.asarray(omega3)
        return self.density*self.b/(self.D(omega1+omega2+omega3)*self.D(omega1)*self.D(omega2)*self.D(omega3))

    def solve(self,time,field,x0=0.,v0=0.,substeps=10,kick_index=None,kick_velocity=0.):
        """
        Solve the equation of motion with a fourth order Runge-Kutta method on the (uniform) time grid. The
        equation can be solved for several runs at the same time.

        Args:
            time (:py:class:`numpy.ndarray`) : uniform time grid (au)
            field (:py:class:`function`) : function of the time that returns the field for each run
            x0, v0 (:py:class:`float` or :py:class:`numpy.ndarray`) : initial position and velocity
            substeps (:py:class:`int`) : number of integration steps for each time step of the grid
            kick_index (:py:class:`numpy.ndarray` or None) : for each run, index of the point of the time grid at which
                the velocity is changed by kick_velocity (-1 for no kick), i.e. a delta-shaped field
                :math:`E_0\\delta(t-t_k)` with kick_velocity = :math:`-E_0`. Default is None (no kick)
            kick_velocity (:py:class:`float`) : change of the velocity of the kick

        Returns:
            :py:class:`numpy.ndarray` : the polarization -N*x, with shape (runs, time)
        """
        def acceleration(t,x,v):
            return -2.*self.gamma*v - self.omega0**2*x - self.a*x**2 + self.b*x**3 - field(t)
        nruns = np.shape(field(time[0]))
        x = np.zeros(nruns) + x0
        v = np.zeros(nruns) + v0
        h = (time[1]-time[0])/substeps
        pos = np.zeros(nruns+(len(time),))
        pos[...,0] = x
        t = time[0]
        for it in range(1,len(time)):
            if kick_index is not None:
                v = v + kick_velocity*(np.asarray(kick_index) == it-1)
            for _ in range(substeps):
                k1x, k1v = v, acceleration(t,x,v)
                k2x, k2v = v+0.5*h*k1v, acceleration(t+0.5*h,x+0.5*h*k1x,v+0.5*h*k1v)
                k3x, k3v = v+0.5*h*k2v, acceleration(t+0.5*h,x+0.5*h*k2x,v+0.5*h*k2v)
                k4x, k4v = v+h*k3v, acceleration(t+h,x+h*k3x,v+h*k3v)
                x = x + h/6.*(k1x+2*k2x+2*k3x+k4x)
                v = v + h/6.*(k1v+2*k2v+2*k3v+k4v)
                t = t + h
            t = time[it]
            pos[...,it] = x
        return -self.density*pos

    def delta_response(self,time,amplitude=1e-3,initial_time=0.,substeps=10):
        """
        Compute the polarization induced by a delta-shaped field :math:`E(t) = E_0\\delta(t-t_0)`, i.e. a kick of
        the velocity of the electron at the initial time of the field. This is the field used by yambo_nl to
        compute the linear response.

        Args:
            time (:py:class:`numpy.ndarray`) : uniform time grid (au)
            amplitude (:py:class:`float`) : amplitude :math:`E_0` of the field (au)
            initial_time (:py:class:`float`) : time of the kick (au), rounded to the closest point of the grid
            substeps (:py:class:`int`) : number of integration steps for each time step of the grid

        Returns:
            :py:class:`NLData` : the time dependent polarization and the parameters of the field
        """
        i0 = int(np.argmin(np.abs(time-initial_time)))
        pol = np.zeros((3,len(time)))
        pol[0,i0:] = self.solve(time[i0:],lambda t : np.zeros(1),v0=-amplitude,substeps=substeps)[0]
        efield = build_field('DELTA',amplitude=amplitude,initial_time=time[i0])
        return NLData(time,[pol],[efield],[build_field('none')],self.gamma)

    def delta_probe_response(self,time,pump,probe_times,amplitude=1e-6,substeps=10):
        r"""
        Compute the polarization induced by a delta-shaped probe :math:`E_0\delta(t-t_p)` in presence of a pump field,
        for each probe time: the difference :math:`P_{Pp}-P_P` between the polarization with probe and pump and the one
        with the pump only, which defines the transient response :math:`\chi^{neq} = (P_{Pp}-P_P)/E_p` of a pump-probe
        experiment. The :py:func:`Linear_Response` of each run gives :math:`1+4\pi\chi^{neq}`.

        Args:
            time (:py:class:`numpy.ndarray`) : uniform time grid (au)
            pump (:py:class:`function`) : function of the time (au) that returns the pump field (au)
            probe_times (:py:class:`float` or :py:class:`numpy.ndarray`) : times of the probe (au), rounded to the
                closest points of the grid
            amplitude (:py:class:`float`) : amplitude :math:`E_0` of the probe (au)
            substeps (:py:class:`int`) : number of integration steps for each time step of the grid

        Returns:
            :py:class:`NLData` : one run for each probe time, with the polarization :math:`P_{Pp}-P_P` and the
            parameters of the (delta-shaped) probe. The pump is not described in the result
        """
        idx = [int(np.argmin(np.abs(time-tp))) for tp in np.atleast_1d(probe_times)]
        nruns = len(idx)+1
        kick_index = np.array([-1]+idx)  # the first run has the pump only
        pol_x = self.solve(time,lambda t : pump(t)*np.ones(nruns),substeps=substeps,kick_index=kick_index,
                           kick_velocity=-amplitude)
        polarization = []
        for p in pol_x[1:]:
            pol = np.zeros((3,len(time)))
            pol[0] = p - pol_x[0]
            polarization.append(pol)
        efields = [build_field('DELTA',amplitude=amplitude,initial_time=time[i]) for i in idx]
        return NLData(time,polarization,efields,[build_field('none')]*len(idx),self.gamma)

    def sine_response(self,time,frequencies,amplitude=1e-3,initial_time=0.,pump=None,substeps=10):
        r"""
        Compute the polarization induced by a sine-shaped field :math:`E(t) = E_0 \sin(\omega(t-t_0))` switched on
        at :math:`t_0`, for each value of the frequency. If the pump is provided a second sine-shaped field (with
        the same structure) is added to the first one, as in a pump and probe computation.

        Args:
            time (:py:class:`numpy.ndarray`) : uniform time grid (au)
            frequencies (:py:class:`numpy.ndarray`) : frequencies of the first field (Hartree)
            amplitude (:py:class:`float`) : amplitude of the first field (au)
            initial_time (:py:class:`float`) : switch on time of the first field (au)
            pump (:py:class:`dict` or None) : dictionary with the keys 'frequency', 'amplitude' and 'initial_time'
                of the second field
            substeps (:py:class:`int`) : number of integration steps for each time step of the grid

        Returns:
            :py:class:`NLData` : the time dependent polarization and the parameters of the fields
        """
        freqs = np.atleast_1d(np.asarray(frequencies,dtype=float))
        def sine(t,w,E0,t0):
            return E0*np.sin(w*(t-t0))*(t >= t0)
        if pump is None:
            field = lambda t : sine(t,freqs,amplitude,initial_time)
            field2 = build_field('none')
        else:
            wP, EP, tP = pump['frequency'], pump['amplitude'], pump.get('initial_time',0.)
            field = lambda t : sine(t,freqs,amplitude,initial_time) + sine(t,wP,EP,tP)
            field2 = build_field('SIN',wP,EP,tP)
        pol_x = self.solve(time,field,substeps=substeps)
        polarization = []
        for p in pol_x:
            pol = np.zeros((3,len(time)))
            pol[0] = p
            polarization.append(pol)
        efields = [build_field('SIN',w,amplitude,initial_time) for w in freqs]
        return NLData(time,polarization,efields,[field2]*len(freqs),self.gamma)
