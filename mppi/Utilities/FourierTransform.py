"""
This module contains a function for the computation of the one-dimensional fourier transform of a
real valued function. The module can be loaded in the notebook in one of the following way

>>> from mppi.Utilities import FourierTransform

>>>FourierTransform.eval_FT

or as

>>> import mppi.Utilities.FourierTransform

>>> FourierTransform.eval_FT

or to load directly some elements

>>> from mppi.Utilities.FourierTransform import eval_FT

>>> eval_FT

"""
from ast import If
from mppi.Utilities import Constants as C
import numpy as np

import os

def eval_FT(time,values,time_units = 'fs', verbose = True):
    """
    Compute the FT of the (real-valued) values array using the time array as x-values.
    The FT is provided only for positive energies since due to the real-valued of the
    values array the FT is assumed to be an even function.

    Args:
        time (:py:class:`np.array`) : x-values array (assumes a uniform sampling)
        values (:py:class:`np.array`) : values array (assumes real-value function)
        time_units (:py:class:`string`) : set the units to compute the energy array.
            Default is 'fs' and the other possible choice is 'ps'
        verbose (:py:class:`bool`) : defines the amount of information provided
            on terminal

    Return:
        :py:class:`dict`  : dictionary with the energy array (in eV), the real and imaginary part
            of the FT and the associated module

    """
    h = C.Planck_ev_ps
    if time_units == 'fs':
        h *= 1e3
        if verbose : print('Time values expressed in fs. Energy array is provided in eV.')
    elif time_units == 'ps':
        if verbose : print('Time values expressed in ps. Energy array is provided in eV.')
    else:
        print('Unkwon time units. FT has not been computed.')
        return 0
    dt = time[1]-time[0]
    Nt = len(time)
    freqs = np.fft.rfftfreq(Nt,d=dt)
    energy = h*freqs   #[0:int(N/2)]
    ft = np.fft.rfft(values)   #[0:int(N/2)]
    ft_abs = np.sqrt(ft.real**2+ft.imag**2)
    if verbose :
        print('Maximum energy value =',energy[-1])
        print('Energy sampling step =',energy[1]-energy[0])
    return {'energy':energy,'ft_real':ft.real,'ft_im':ft.imag,'ft_abs':ft_abs}

def set_time_sampling_pars(e_max,de,time_units='fs',verbose=True):
    """
    Compute the time sampling parameters (time step and number of points) to obtain a FT with a given maximum 
    energy and energy sampling step.
    The function assumes a uniform time sampling and that the FT is computed with the numpy.fft.rfft function, the
    number of sampling points in the time domain is chosen as even.
    """
    if time_units == 'fs':
        h = C.Planck_ev_fs
    elif time_units == 'ps':
        h = C.Planck_ev_ps
    else:
        print('Unknown time units. Time sampling parameters have not been computed.')
        return 0
    dt = h/(2*e_max)
    Nt = int(h/(dt*de))
    if Nt%2 != 0:
        Nt += 1
    if verbose: 
        print('Time step =',dt,)
        print('Number of time points =',Nt)
    return dt, Nt

def eval_energy_array(Nt,dt=1.,time_units='fs',verbose=True):
    """
    Compute the energy array associated to a time sampling with Nt points and time step dt.
    The energy is given by energy = h*freqs, where freqs is the array of frequencies 
    is computedf by the np.fft.rfftfreq function.

    Args:
        Nt (:py:class:`int`) : number of sampling points in the time domain
        dt (:py:class:`float`) : time step of the sampling in the time domain.Default is 1.0.
        time_units (:py:class:`string`) : set the units to time sampling.  
            Default is 'fs' and the other possible choice is 'au' (atomic units)

    Return:
        :py:class:`np.array` : energy array associated to the time sampling. 
            The energy is given in eV if time_units is 'fs' and in Hartree if time_units is 'au'.

    """
    if time_units == 'fs':
        h = C.Planck_ev_fs
    elif time_units == 'au':
        h = C.Planck_au
    else:
        print('Unknown time units. Energy array has not been computed.')
        return 0
    freqs = np.fft.rfftfreq(Nt,d=dt)
    energy = h*freqs
    if time_units == 'fs':
        energy = C.Planck_ev_fs*freqs
    elif time_units == 'au':
        energy = C.Planck_au*freqs
    if verbose :
        print('Maximum energy value =',energy[-1])
        print('Energy sampling step =',energy[1]-energy[0])
    return energy