"""
This module collects some tools to deal with Gaussian pulses.
This class of pulses is described by a sinusoidal function, whose frequency
defines the `energy` :math:`E=\hbar\omega` of the pulse, times a gaussian envelope.
The time-dependent amplitude of the pulse is parametrized as

.. math::
   :label: pulseamplitude

   {\cal G}(t) = A sin(\omega t)e^{-(t-t_0)^2/2w^2}

where ``A`` is the amplitude of the pulse and the value of :math:`t_0` sets the
position of the maximum of the gaussian envelope. Here :math:`t_0` is expressed
in function of the width parameter as the multiple of the half period of the sine
function nearest to 3w. Due to this choice the maximum of the gaussian profile
matches with the maximum of the sine function. Finally, the origin of the time
variable can be set using the :math:`t_{start}` parameter.

The width parameter determines the value of the FWHM of the pulse according to the
relation

.. math::
   :label: fwhm_def

   FWHM = 2 w \sqrt{2ln(2)}

"""
import numpy as np
from mppi.Utilities import Constants as C

def gaussianPulse(time, energy = 1.5, amplitude = 1, width = 100, fwhm = None,t_start = 0,
        envelope_only = False, THz_pulse = False, verbose = True):
    """
    Build a Gaussian pulse.

    Args:
        time : the value of the time variable (in fs), it can be a single :py:class:`float`
            or an :py:class:`array`
        energy (:py:class:`float`) : the energy of the pulse (in eV)
        amplitude (:py:class:`float`) : the amplitude of the pulse (it can be a complex number)
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
        with the value of the Gaussian pulse. The pulse can be complex (depending on the
        amplitude provided as input)

    """
    if not THz_pulse:
        h_red = C.Planck_reduced_ev_ps*1e3 # in eV*fs
        timeUnit = 'fs'
        inverseTimeUnit = 'fs^-1'
        if verbose: print('time unit: fs - energy unit: eV')
    else:
        h_red = C.Planck_reduced_ev_ps*1e3 # in meV*ps
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
    t = time-t_start
    t0 = np.pi/omega*round(omega/np.pi*3*width)
    #t0 = 3.*width
    exp_arg = (t-t0)**2/(2.*width**2)
    pulse = amplitude*np.exp(-exp_arg)
    if not envelope_only:
        pulse *= np.sin(omega*t)
    return pulse

def doubleGaussianPulse(time, energy = 1.5, amplitude1 = 1, amplitude2 = 1, width1 = 100, width2 = 100,
                        fwhm1 = None, fwhm2 = None, t_start1 = 0, t_start2 = 600, envelope_only = False,
                        THz_pulse = False, verbose = True):
    """
    Build a sine oscillating pulse with e double Gaussian envelope.

    Args:
        time : the value of the time variable (in fs), it can be a single :py:class:`float`
            or an :py:class:`array`
        energy (:py:class:`float`) : the energy of the pulse (in eV)
        amplitude1 (:py:class:`float`) : the amplitude of the first pulse (it can be a complex number)
        amplitude2 (:py:class:`float`) : the amplitude of the second pulse (it can be a complex number)
        width1 (:py:class:`float`) : the width parameter of the first Gaussian (in fs)
        fwhm1 (:py:class:`float`) : if not None set the FWHM of the first pulse (in fs). The width1
            is set to :math:`fwhm1/(2\sqrt{2ln(2)})`
        width2 (:py:class:`float`) : the width parameter of the second Gaussian (in fs)
        fwhm2 (:py:class:`float`) : if not None set the FWHM of the second pulse (in fs). The width2
            is set to :math:`fwhm2/(2\sqrt{2ln(2)})`
        t_start1 (:py:class:`float`) : time shift for the origin of the first pulse (in fs)
        t_start2 (:py:class:`float`) : time shift for the origin of the second pulse (in fs)
        envelope_only (:py:class:`bool`) : if True the sinusodial oscillating factor
            is not considered
        THz_pulse (:py:class:`bool`) : if True expresses the energy in meV and the
            time variable and the width parameter in ps
        verbose (:py:class:`bool`) : defines the amount of information provided
            on terminal

    Returns:
        :py:class:`float` or :py:class:`array` (depending on the time provided as input)
        with the value of the double Gaussian pulse. The pulse can be complex (depending on the
        amplitude provided as input)

    """
    if not THz_pulse:
        h_red = C.Planck_reduced_ev_ps*1e3 # in eV*fs
        timeUnit = 'fs'
        inverseTimeUnit = 'fs^-1'
        if verbose: print('time unit: fs - energy unit: eV')
    else:
        h_red = C.Planck_reduced_ev_ps*1e3 # in meV*ps
        timeUnit = 'ps'
        inverseTimeUnit = 'ps^-1'
        if verbose: print('time unit: ps - energy unit: meV')
    omega = energy/h_red
    if fwhm1 is None:
        fwhm1 = width1*(2.*np.sqrt(2.*np.log(2.)))
    else:
        width1 = fwhm1/(2.*np.sqrt(2.*np.log(2.)))
    if fwhm2 is None:
        fwhm2 = width2*(2.*np.sqrt(2.*np.log(2.)))
    else:
        width2 = fwhm2/(2.*np.sqrt(2.*np.log(2.)))
    if verbose :
        print('period of the oscillations',2.*np.pi/omega,timeUnit)
        print('width of the first pulse',width1,timeUnit)
        print('fwhm of the first pulse',fwhm1,timeUnit)
        print('width of the second pulse',width2,timeUnit)
        print('fwhm of the second pulse',fwhm2,timeUnit)
    pulse1 = gaussianPulse(time, energy=energy, amplitude=amplitude1,width=width1,
            t_start= t_start1,envelope_only=True,THz_pulse=THz_pulse,verbose=False)
    pulse2 = gaussianPulse(time, energy=energy, amplitude=amplitude2,width=width2,
            t_start= t_start2,envelope_only=True,THz_pulse=THz_pulse,verbose=False)
    doublePulse = pulse1 + pulse2
    t = time-t_start1
    if not envelope_only:
        doublePulse *= np.sin(omega*t)
    return doublePulse

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
        h = C.Planck_ev_ps*1e3 # in eV*fs
        energyUnit = 'eV'
        if verbose: print('time unit: fs - energy unit: eV')
    else:
        h = C.Planck_ev_ps*1e3 # in meV*ps
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
