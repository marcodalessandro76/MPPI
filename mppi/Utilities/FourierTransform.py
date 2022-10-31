"""
This module a function for the computation of the one-dimensional fourier transform of a
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
    from mppi.Utilities import Constants as C
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
    N = len(time)
    freqs = np.fft.fftfreq(N,d=dt)
    energy = h*freqs[0:int(N/2)]
    ft = np.fft.fft(values)[0:int(N/2)]
    ft_abs = np.sqrt(ft.real**2+ft.imag**2)
    if verbose :
        print('Maximum energy value =',energy[-1])
        print('Energy sampling step =',energy[1]-energy[0])
    return {'energy':energy,'ft_real':ft.real,'ft_im':ft.imag,'ft_abs':ft_abs}
