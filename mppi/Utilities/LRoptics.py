"""
This module defines the tools to extract the linear optical properties of the system from the real-time polarization.
The module can be loaded in the notebook as follows

>>> from mppi.Utilities import LRoptics as LR

>>> LR.Linear_Response

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
    
#####################################################
##### TO BE REMOVED #################################

import sys
from matplotlib import pyplot as plt
ha2ev = 27.211386

def Fourier_Interpolation(ft, fw, time,freqs,mode="T2W"):
    #
    # I assume constant time-step and regular frequency grid
    # otherwise this subroutine needs to be generalized
    #
    t_step=time[1]-time[0]
    f_step=freqs[1]-freqs[0]
    #
    if mode.upper() == "T2W":
        fw[:,:]=0.0+0.0j
        for i_w in range(fw.shape[1]):
            for i_c in range(ft.shape[0]):
                    fw[i_c,i_w] = fw[i_c,i_w]+np.sum(ft[i_c,:]*np.exp(1j*freqs[i_w]*time[:]))*t_step
    elif mode.upper() == "W2T":
        ft[:,:]=0.0+0.0j
        for i_w,i_t in zip(range(ft.shape[1]),range(fw.shape[1])):
            for i_c in range(ft.shape[0]):
                ft[i_c,i_t] = ft[i_c,i_t]+fw[i_c,i_w]*np.exp(-1j*freqs[i_w]*time[i_t])*f_step


def damp_it(ft, time, t_initial, damp_type="LORENTZIAN", damp_factor=0.1/ha2ev):
    ft_damped=np.empty_like(ft)

    if damp_type.upper() == "LORENTZIAN":
        ft_damped[:]=ft[:]*np.exp(-abs(time[:]-t_initial)*damp_factor)
    elif damp_type.upper() == "GAUSSIAN":
        ft_damped[:]=ft[:]*np.exp(-(time[:]-t_initial)**2*damp_factor**2)
    else:
        print("Wrong damping type ")
        sys.exit(0)

    return ft_damped

def get_Efield_w(freqs,efield):
    if efield["name"] == "DELTA":
        print('t_initial =',efield["initial_time"])
        efield_w=efield["amplitude"]*np.exp(1j*freqs[:]*efield["initial_time"])
    else:
        print("Fields different from Delta function not implemented yet")
        sys.exit(0)

    return efield_w

def YamboPy_Linear_Response(time, pol, efield, pol_ref=None, e_range=[0.0, 20.0], n_freqs=200,plot=True):
    if efield["name"] != "DELTA":
        print("Linear response implemented only for Delta function external fields ")
        sys.exit(0)

    # Pump and probe minus the pump if present
    if pol_ref is not None:
        pol=pol-pol_ref

    freqs=np.linspace(e_range[0],e_range[1],n_freqs)/ha2ev
    pol_w=np.zeros((pol.shape[0],n_freqs),dtype=complex)

    efield_w=get_Efield_w(freqs,efield)

    Fourier_Interpolation(pol,pol_w,time,freqs,mode="T2W")
    #
    # Polarization along the field direction
    #
    pol_w_along_E=np.zeros(n_freqs,dtype=complex)
    for i_d in range(3):
        pol_w_along_E[:]+=pol_w[i_d,:]*efield["versor"][i_d]
    #
    # EPS = 1+4.0*pi*Xhi(omega)
    #
    eps=np.zeros(n_freqs,dtype=complex)
    eps=1.0+4.0*np.pi*pol_w_along_E/efield_w

    fig, axs = plt.subplots(2)
    fig.suptitle('Dielectric constant along the field direction')
    axs[0].set_xlim(e_range)
    axs[0].plot(freqs*ha2ev,eps.real,label='Real part ')
    axs[0].legend()
    axs[0].set_xticks([])
    axs[1].set_xlim(e_range)
    axs[1].plot(freqs*ha2ev,eps.imag,label='Imag part ')
    axs[1].legend()
    axs[1].set_xlabel('eV')
    if plot:
        plt.show()