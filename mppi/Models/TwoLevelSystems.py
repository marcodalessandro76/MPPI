"""
This module collects some tools to deal with the physics of Two Level Systems (TLS).
This modeling is particularly useful to describe both electromagnetic transition
between two *optically active* bands and magnetic system in which a time-dependent
magnetic field induces oscillations between two (non degenerate) states with different
spin.

The feature of the model, as well as the derivation of the equation used in this
implementation of the physics of the TLS are described in this notebook_.

.. _notebook: tutorials/Model_TLS_optical_absorption.ipynb

"""

import numpy as np
from mppi.Utilities import Constants as C
from scipy.integrate import odeint

def pulseParametersFromIntensity(mu12, intensity, width = 100, fwhm = None,
        THz_pulse = False, verbose = True):
    """
    Compute the Rabi coupling frequency and the `area` of the Gaussian pulse in function of
    the values of the transition dipole :math:`\mu_{12}` and of the field intensity.

    Args:
        mu12 (:py:class:`complex`) : complex value of the transition dipole :math:`\mu_{12}` that enters
            in the definition of the Rabi coupling
        intensity (:py:class:`float`) : field intensity in kW/cm^2
        width (:py:class:`float`) : width parameter of the Gaussian pulse (in fs)
        fwhm (:py:class:`float`) : if not None set the FWHM of the pulse (in fs). The width
            is set to :math:`fwhm/(2\sqrt{2ln(2)})`
        THz_pulse (:py:class:`bool`) : if True expresses the width and the fwhm parameters in ps
        verbose (:py:class:`bool`) : sets the amount of information provided on terminal

    Return:
        :py:class:`dict` : a dictionary withe the (generally complex) Rabi coupling,
                    the field amplitude (in V/m) and the pulse area

    """
    Z0 = C.vacuum_impedence
    if not THz_pulse:
        h_red = C.Planck_reduced_ev_ps*1e3 # in eV*fs
        timeUnit = 'fs'
        inverseTimeUnit = 'fs^-1'
        if verbose: print('time unit: fs')
    else:
        h_red = C.Planck_reduced_ev_ps # in eV*ps
        timeUnit = 'ps'
        inverseTimeUnit = 'ps^-1'
        if verbose: print('time unit: ps')
    if fwhm is not None:
        width = fwhm/(2.*np.sqrt(2.*np.log(2.)))
        if verbose: print('the width parameter of the pulse is',width,timeUnit)
    intensity = intensity*1e3*1e4 #W/m^2
    amplitude_vm = np.sqrt(Z0*intensity) #V/m
    amplitude = amplitude_vm*C.Bohr_radius #V/a0 in atomic units
    Omega0 = mu12*amplitude/h_red

    theta = np.sqrt(2*np.pi)*width*abs(Omega0)
    if verbose:
        print('Rabi coupling (%s):'%inverseTimeUnit,Omega0)
        print('Rabi coupling (module) (%s):'%inverseTimeUnit,abs(Omega0))
        print('field amplitude (V/m):',amplitude_vm)
        print('pulse area :',theta)
    return dict(Omega0=Omega0,amplitude=amplitude_vm,theta=theta)

def pulseParametersFromTheta(mu12, theta, width = 100, fwhm = None, THz_pulse = False, verbose=True):
    """
    Compute the field intensity and the Rabi coupling that correspond to the `area` of the Gaussian
    pulse given as input.

    Args:
        mu12 (:py:class:`complex`) : complex value of the transition dipole :math:`\mu_{12}` that enters
            in the definition of the Rabi coupling
        theta : pulse area
        width (:py:class:`float`) : width parameter of the Gaussian pulse (in fs)
        fwhm (:py:class:`float`) : if not None set the FWHM of the pulse (in fs). The width
            is set to :math:`fwhm/(2\sqrt{2ln(2)})`
        THz_pulse (:py:class:`bool`) : if True expresses the width and the fwhm parameters in ps
        verbose (:py:class:`bool`) : sets the amount of information provided on terminal

    Returns:
        :py:class:`dict` : a dictionay with the (generally complex) Rabi coupling,
            the field amplitude (in V/m) and the field intensity (in kW/cm^2)

    """
    Z0 = C.vacuum_impedence
    if not THz_pulse:
        h_red = C.Planck_reduced_ev_ps*1e3 # in eV*fs
        timeUnit = 'fs'
        inverseTimeUnit = 'fs^-1'
        if verbose: print('time unit: fs')
    else:
        h_red = C.Planck_reduced_ev_ps # in eV*ps
        timeUnit = 'ps'
        inverseTimeUnit = 'ps^-1'
        if verbose: print('time unit: ps')
    if fwhm is not None:
        width = fwhm/(2.*np.sqrt(2.*np.log(2.)))
        if verbose: print('the width parameter of the pulse is',width,timeUnit)
    Omega0_abs = theta/(np.sqrt(2*np.pi)*width) #fs^-1 if THz_pulse is False
    amplitude = Omega0_abs*h_red/abs(mu12) #V/a0 in atomic units
    amplitude_vm = amplitude/C.Bohr_radius #V/m
    intensity = amplitude_vm**2/Z0 #W/m^2
    intensity = intensity*1e-3*1e-4 #kW/cm^2
    Omega0 = mu12*amplitude/h_red
    if verbose:
        print('Rabi coupling (%s):'%inverseTimeUnit,Omega0)
        print('Rabi coupling (module) (%s):'%inverseTimeUnit,abs(Omega0))
        print('field amplitude (V/m):',amplitude_vm)
        print('field intensity (kW/cm^2) :',intensity)
    return dict(Omega0=Omega0,amplitude=amplitude_vm,intensity=intensity)

def solveBlochEq(u0, time, Omega, delta = 0, eta = 0):
    """
    Solve the Bloch equations in the rotating frame.

    Args:
        x0 (:py:class:`array`) : initial condition
        time (:py:class:`array`) : array with the time values
        Omega (`function`) : (complex) envelope function
        delta (:py:class:`float`) : difference between pulse energy and the energy shift
            of the two levels of the system. It has to be provided in the inverse dimensions
            of the time variable
        eta (:py:class:`float`) : dephasing parameter for the polarization (u1 and u2 elements of
            the Bloch vector). It has to be provided in the inverse dimensions
            of the time variable

    Returns:
        (:py:class:`array`) : array with the solution of the Bloch equations.
                The index shape of the array is (3,len(time))

    """
    def Bloch_Eq(u, t, delta, Omega):
        dudt = [-delta*u[1]-np.real(Omega(t))*u[2]-eta*u[0],\
                 delta*u[0]+np.imag(Omega(t))*u[2]-eta*u[1],\
                 np.real(Omega(t))*u[0]-np.imag(Omega(t))*u[1]]
        return dudt
    u = odeint(Bloch_Eq, u0, time, args=(delta, Omega))
    return u.transpose()

def solveBlochEq_Real(u0, time, Omega, delta = 0, eta = 0):
    """
    Solve the Bloch equations in the rotating frame. This function assumes that
    the provided Rabi coupling is real and positive.

    Args:
        u0 (:py:class:`array`) : initial condition
        time (:py:class:`array`) : array with the time values
        Omega (`function`) : (real) envelope function built in terms of
            the module of the Rabi coupling
        delta (:py:class:`float`) : difference between pulse energy and the energy shift
            of the two levels of the system. It has to be provided in the inverse dimensions
            of the time variable
        eta (:py:class:`float`) : dephasing parameter for the polarization (u1 and u2 elements of
            the Bloch vector). It has to be provided in the inverse dimensions
            of the time variable

    Returns:
        (:py:class:`array`) : array with the solution of the Bloch equations.
                The index shape of the array is (3,len(time))

    """
    def Bloch_Eq_Real(u, t, delta, Omega):
        dudt = [-delta*u[1]-Omega(t)*u[2]-eta*u[0],\
                 delta*u[0]-eta*u[1],\
                 Omega(t)*u[0]]
        return dudt
    u = odeint(Bloch_Eq_Real, u0, time, args=(delta, Omega))
    return u.transpose()

def convertToRotatingFrame(omega, time, u, invert = False):
    """
    Convert the Bloch vector components to the rotating frame (from the u_i to the u'_i).
    If the ``invert`` argument is True perform the inverse change of basis (from the u'_i to the u_i).

    Args:
        omega (:py:class:`float`) : angular frequency of the pulse, expressed in the
            reciprocal unit of the time variable
        time (:py:class:`array`) : array with the time values
        u (:py:class:`array`) : array with the Bloch vector in the original frame. The
            function assumes that the array has two indices, the first is the component
            and the second one is the time index
        invert (:py:class:`bool`) : if True perform the inverse transformation

    Returns:
        (:py:class:`array`) : array with the transformed Bloch vector

    """
    uprime = np.zeros([3,len(time)])
    if invert: alpha = -1
    else: alpha = 1
    uprime[0,:] = np.cos(omega*time)*u[0,:] - alpha*np.sin(omega*time)*u[1,:]
    uprime[1,:] = alpha*np.sin(omega*time)*u[0,:] + np.cos(omega*time)*u[1,:]
    uprime[2,:] = u[2,:]
    return uprime
