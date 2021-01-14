"""
This module collects some tools to deal with Gaussian pulses.
This class of pulses is described by a sinusoidal function, whose frequency
defines the `energy` :math:`E=\hbar\omega` of the pulse, times a gaussian envelope.
The time-dependent amplitude of the pulse is parametrized as

.. math::
   :label: pulseamplitude

   {\cal G}(t) = A sin(\omega t)e^{-(t-(t_0+t_{start}))^2/2w^2}

where ``A`` is the amplitude of the pulse. The :math:`t_start` parameter shifts the pulse
along the time interval and the value of :math:`t_0` sets the position of the maximum of the
gaussian envelope. Here :math:`t0` is expressed in function of the width parameter`
as :math:`t_0=3w`.

The width parameter determines the value of the FWHM of the pulse according to the
relation

.. math::
   :label: fwhm_def

   FWHM = 2 w \sqrt{2ln(2)}

The Gaussian pulse contributes to the interaction Hamiltonian :math:`H_I` through the
dipole interaction. The amplitude parameter of :math:`H_I` is expressed by the
Rabi coupling

.. math::
   :label: rabi_coupl

   \Omega_r = \\frac{Ad}{\hbar}

where `d` is the value of the dipole parallel to the field. The Rabi coupling
determines the :math:`\\theta` parameter of the field, that represents the `area`
of the envelope of the pulse

.. math::
    \\theta = \int_{-\infty}^{\infty}\Omega_r e^{-t^2/2w^2} = \sqrt{2\pi}w\Omega_r

According to the physics of the optical obsorption in the two level systems, it is
useful to determine the value of the Rabi coupling for which the :math:`\\theta`
parameter has some specific values, like :math:`\pi` or :math:`\pi/2`. This
functionality is emplemented in the module.

"""
import numpy as np
from mppi import Utilities as U

def gaussianPulse(time, energy = 1.5, amplitude = 1, width = 100, fwhm = None,t_start = 0,
        envelope_only = False, THz_pulse = False, verbose = True, **kwargs):
    """
    Build a Gaussian pulse.

    Args:
        time : the value of the time variable (in fs), it can be a single :py:class:`float`
            or an :py:class:`array`
        energy (:py:class:`float`) : the energy of the pulse (in eV)
        amplitude (:py:class:`float`) : the amplitude of the pulse
        width (:py:class:`float`) : the width parameter of the Gaussian (in fs)
        fwhm (:py:class:`float`) : if not None set the FWHM of the pulse (in fs). The width
            is set to :math:`fwhm/(2\sqrt{2ln(2)})`
        t_start (:py:class:`float`) : time shift for the origin of the pulse (in fs)
        envelope_only (:py:class:`bool`) : if True the sinusodial oscillating factor
            is not considered
        THz_pulse (:py:class:`bool`) : if True expresses the energy in meV and the
            time variable and the width parameter in ps
        verbose (:py:class:`bool`) : defines the amount of information provided
            on terminal

    Returns:
        :py:class:`float` or :py:class:`array` (depending on the time provided as input)
        with the value of the Gaussian pulse

    """
    if not THz_pulse:
        h_red = U.Planck_reduced_ev_ps*1e3 # in eV*fs
        timeUnit = 'fs'
        inverseTimeUnit = 'fs^-1'
        if verbose: print('time unit: fs - energy unit: eV')
    else:
        h_red = U.Planck_reduced_ev_ps*1e3 # in meV*ps
        timeUnit = 'ps'
        inverseTimeUnit = 'ps^-1'
        if verbose: print('time unit: ps - energy unit: meV')
    omega = energy/h_red
    if fwhm is None:
        fwhm = width*(2.*np.sqrt(2.*np.log(2.)))
    else:
        width = fwhm/(2.*np.sqrt(2.*np.log(2.)))
    if verbose :
        print('period of the oscillations',2.*np.pi/omega,timeUnit)
        print('width of the pulse',width,timeUnit)
        print('fwhm of the pulse',fwhm,timeUnit)
    exp_arg = (time-2*t_start-3.*width)**2/(2.*width**2)
    pulse = amplitude*np.exp(-exp_arg)
    if not envelope_only:
        pulse *= np.sin(omega*(time-t_start))
    return pulse

def pulseParametersFromIntensity(dipole, intensity, width = 100, fwhm = None,
        THz_pulse = False, verbose = True):
    """
    Compute the Rabi coupling frequency and the pulse area in function of the values of the
    transition dipole and of the field intensity.

    Args:
        dipole (:py:class:`array`) : array with the transition dipole (real and imaginary part),
            as provided by the :class:YamboDipolesParser class
        intensity (:py:class:`float`) : field intensity in kW/cm^2
        width (:py:class:`float`) : width parameter of the Gaussian pulse (in fs)
        fwhm (:py:class:`float`) : if not None set the FWHM of the pulse (in fs). The width
            is set to :math:`fwhm/(2\sqrt{2ln(2)})`
        THz_pulse (:py:class:`bool`) : if True expresses the width and the fwhm parameters in ps
        verbose (:py:class:`bool`) : sets the amount of information provided on terminal

    Returns:
        :py:class:`tuple` : a tuple (Omega_r,amplitude_vm,theta) with the Rabi coupling
            frequency (in fs^-1), the field amplitude (in V/m) and the pulse area

    """
    Z0 = U.vacuum_impedence
    if not THz_pulse:
        h_red = U.Planck_reduced_ev_ps*1e3 # in eV*fs
        timeUnit = 'fs'
        inverseTimeUnit = 'fs^-1'
        if verbose: print('time unit: fs')
    else:
        h_red = U.Planck_reduced_ev_ps # in eV*ps
        timeUnit = 'ps'
        inverseTimeUnit = 'ps^-1'
        if verbose: print('time unit: ps')
    if fwhm is not None:
        width = fwhm/(2.*np.sqrt(2.*np.log(2.)))
        if verbose: print('set width to',width,timeUnit)
    intensity = intensity*1e3*1e4 #W/m^2
    amplitude_vm = np.sqrt(Z0*intensity) #V/m
    amplitude = amplitude_vm*U.Bohr_radius #V/a0 in atomic units
    dip_mod = np.linalg.norm(dipole)
    Omega_r = dip_mod*amplitude/h_red

    theta = np.sqrt(2*np.pi)*width*Omega_r
    if verbose:
        print('coupling frequency (%s):'%inverseTimeUnit,Omega_r)
        print('field amplitude (V/m):',amplitude_vm)
        print('pulse area :',theta)
    return (Omega_r,amplitude_vm,theta)

def pulseParametersFromTheta(dipole, theta, width = 100, fwhm = None, THz_pulse = False, verbose=True):
    """
    Compute the field intensity and the Rabi coupling that correspond to the pulse area given as input.

    Args:
        dipole (:py:class:`array`) : array with the transition dipole (real and imaginary part),
            as provided by the :class:YamboDipolesParser class
        theta : pulse area
        width (:py:class:`float`) : width parameter of the Gaussian pulse (in fs)
        fwhm (:py:class:`float`) : if not None set the FWHM of the pulse (in fs). The width
            is set to :math:`fwhm/(2\sqrt{2ln(2)})`
        THz_pulse (:py:class:`bool`) : if True expresses the width and the fwhm parameters in ps
        verbose (:py:class:`bool`) : sets the amount of information provided on terminal

    Returns:
        :py:class:`tuple` : a tuple (Omega_r,amplitude_vm,intensity) with the Rabi coupling (in fs^-1),
            the field amplitude (in V/m) and the field intensity (in kW/cm^2) that correspond
            to the pulse area provided in input

    """
    Z0 = U.vacuum_impedence
    if not THz_pulse:
        h_red = U.Planck_reduced_ev_ps*1e3 # in eV*fs
        timeUnit = 'fs'
        inverseTimeUnit = 'fs^-1'
        if verbose: print('time unit: fs')
    else:
        h_red = U.Planck_reduced_ev_ps # in eV*ps
        timeUnit = 'ps'
        inverseTimeUnit = 'ps^-1'
        if verbose: print('time unit: ps')
    if fwhm is not None:
        width = fwhm/(2.*np.sqrt(2.*np.log(2.)))
        if verbose: print('set width to',width,timeUnit)
    dip_mod = np.linalg.norm(dipole)
    Omega_r = theta/(np.sqrt(2*np.pi)*width) #fs^-1 if THz_pulse is False
    amplitude = Omega_r*h_red/dip_mod #V/a0 in atomic units
    amplitude_vm = amplitude/U.Bohr_radius #V/m
    intensity = amplitude_vm**2/Z0 #W/m^2
    intensity = intensity*1e-3*1e-4 #kW/cm^2
    if verbose:
        print('coupling frequency (%s):'%inverseTimeUnit,Omega_r)
        print('field amplitude (V/m):',amplitude_vm)
        print('field intensity (kW/cm^2) :',intensity)
    return (Omega_r,amplitude_vm,intensity)

def evalPulseFourierTransform(time, pulse, THz_pulse = False, verbose = True):
    """
    Eval the Fourier transform of the pulse.

    Args:
        time (:py:class:`array`) : array with the values of the time variable (in fs).
        The function assumes that the dt of the array is constant
        pulse (:py:class:`array`) : array with the pulse
        THz_pulse (:py:class:`bool`) : if True assumes that the time variable is provided in
            ps and expresses the energy in meV
        verbose (:py:class:`bool`) : sets the amount of information provided on terminal

    Returns:
        :py:class:`tuple` : a tuple (energies,pulseFTmod) with the (positive range) time conjugate
        variable converted in eV (or in meV, depending on the ``THz_pulse`` option) and the modulus of
        the FT ot the pulse

    Note:
        To increase the resolution in the frequency space choose a long time interval

    """
    if not THz_pulse:
        h = U.Planck_ev_ps*1e3 # in eV*fs
        energyUnit = 'eV'
        if verbose: print('time unit: fs - energy unit: eV')
    else:
        h = U.Planck_ev_ps*1e3 # in meV*ps
        energyUnit = 'meV'
        if verbose: print('time unit: ps - energy unit: meV')
    dt = time[1]-time[0]
    N = len(time)
    freqs = np.fft.fftfreq(N,d=dt)
    energies = h*freqs[0:int(N/2)]
    if verbose:
        de = energies[1]-energies[0]
        print('energy resolution of the FT:',de,energyUnit)
        print('maximum energy:',energies[-1],energyUnit)
    pulseFT = np.fft.fft(pulse)[0:int(N/2)]
    pulseFTmod = np.sqrt(pulseFT.real**2+pulseFT.imag**2)
    if verbose: # compute the FWHM of the FT
        M = max(pulseFTmod)
        ind = np.where(pulseFTmod >= 0.5*M)[0]
        fwhm = energies[ind[-1]]-energies[ind[0]]
        print('FWHM of the FT of the pulse:',fwhm,energyUnit)
    return (energies,pulseFTmod)
