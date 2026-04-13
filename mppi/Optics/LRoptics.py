"""
This module defines the tools to extract the linear optical properties of the system from the real-time polarization induced by a delta-shaped external field.
The module can be loaded in the notebook as follows

>>> from mppi.Optics import LRoptics as LR

>>> LR.Linear_Response

or the method can be imported directly as

>>> from mppi.Optics.LRoptics import Linear_Response

"""
import numpy as np
from mppi.Utilities import FourierTransform as FT
from mppi.Utilities import Constants as C
from mppi.Utilities.Utils import damp_ft

def eval_Efield_w(energy,efield):
    """
    Compute the external field in the frequency domain.
    Only the case of a delta-shaped field is implemented, for which the FT is given by E(w)=E0*exp(-i*w*t_initial), 
    where E0 is the amplitude of the field and t_initial is the switch on time of the field.
    The exponent corresponds to exp(i*E*t_0/hbar), where hbar = 1 in atomic units, so w has to be expressed in Hartree and
    t_initial in atomic units as well. 

    Args:
        freqs (:py:class:`numpy.ndarray`): array with the energy values expressed in Hartree
        efield (:py:class:`dict`): dictionary with the information on the external field used in the calculation.
    
    Return:
        :py:class:`numpy.ndarray`: array with the external field in the frequency domain
    """
    if efield['name'] == 'DELTA':
        efield_w=efield['amplitude']*np.exp(1j*energy[:]*efield['initial_time'])
    else:
        print('Fields different from Delta function not implemented yet')
        return 0


    return efield_w

def Linear_Response(time, pol,efield, pol_ref=None, damp_type="LORENTZIAN", eta=0.1,time_units='au'):
    """
    Compute the linear response of the system to an external delta-shaped field from the polarization in the time domain.
    The function subtract the reference polarization (if provided) and applies a damping function, if eta is not zero.
    The FT is computed throughthe numpy.fft.rfft function assuming  that the time array is uniformly sampled.

    Args:
        time (:py:class:`numpy.ndarray`) : array with the time values
        pol (:py:class:`numpy.ndarray`) : array with the polarization of the system at each time step
        efield (:py:class:`dict`) : dictionary with the information on the external field used in the calculation. 
        pol_ref (:py:class:`numpy.ndarray`) : array with the reference polarization of the system at each time step. Default is None.
        damp_type (:py:class:`str`) : The type of damping function to apply. Can be "LORENTZIAN" or "GAUSSIAN". Default is "LORENTZIAN".
        eta (:py:class:`float`) : The damping factor to apply expressed in eV. Default is 0.1.
        time_units (:py:class:`str`) : set the units to time sampling. Default is 'au' and the other possible choice is 'fs'.
    
    Returns:
         :py:class:`tuple`  : tuple with the energy array (in eV) and the complex dielectric function along the field direction
    """
    
    if efield["name"] != "DELTA":
        print("Linear response implemented only for Delta function external fields ")
        return 0
    
    t = time.copy()
    if time_units == 'fs':
        print('Time units = fs. Rescaled to au')
        t *= C.FsToAu
    elif time_units != 'au':
        print("Invalid time units. Please use 'fs' or 'au'.")
        return 0

    Nt = t.size
    dt = t[1]-t[0]
    energy = FT.eval_energy_array(Nt,dt,time_units='au',verbose=False) 
    t_initial = efield['initial_time']

    if pol_ref is not None:
        pol=pol-pol_ref
    pol_damped=damp_ft(pol,t,t_initial,damp_type=damp_type,eta=eta,time_units='au')
    P_w = 2*np.conjugate(np.fft.rfft(pol_damped))
    efield_w = eval_Efield_w(energy,efield)    
    pol_w_along_E=np.zeros(energy.size,dtype=complex)
    for i_d in range(3):
        pol_w_along_E[:]+=P_w[i_d,:]*efield["versor"][i_d]   
    eps=1.0+4.0*np.pi*pol_w_along_E/efield_w

    return C.HaToeV*energy, eps                  