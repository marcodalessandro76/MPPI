"""
This module contains useful physical constant and conversion factors
"""
import numpy as np

HaToeV = 27.211386
"""
Energy conversion factor from Hartree to eV units
"""

Planck_au = 2.*np.pi
"""
Planck constant in Atomic units. (hbar = 1 in atomic units).
"""

Planck_ev_ps = 4.135667696e-3
"""
Planck constant in eV*ps.
"""

Planck_ev_fs = 4.135667696
"""
Planck constant in eV*fs.
"""

Planck_reduced_ev_ps = 6.58211957e-4
"""
Reduced Planck constant (h/2pi) in eV*ps.
"""

Planck_reduced_ev_fs = 6.58211957e-1
"""
Reduced Planck constant (h/2pi) in eV*fs.
"""

Light_speed_nm_fsm1 = 299.792458
"""
Light speed in nm/fs
"""

Light_speed_au = 137.03599911
"""
Light speed in atomic units (Bohr/fs)
"""

FsToAu = 41.341373
"""
Femto-second to atomic units conversion factor
"""

Bohr_radius = 5.291772e-11
"""
The Bohr radius in meter.
"""

Bohr_to_Angstrom = 0.529177249
"""
Conversion factor from Bohr (the length units in a.u.) to Angstrom.
"""

electron_charge = 1.60217662e-19
"""
The charge of the electron in Coulomb
"""

vacuum_impedence = 376.730313
"""
The impedence of free space in Ohm. The field amplitude E (in V/m) is related to
the field intensity P (in W/m^2) by the relation E = sqrt(Z0*P)
"""

SVCMm1_to_VMm1 =29.98*np.power(10,3,dtype=np.double)
"""
Conversion factor from [statV/cm] to [V/m]
"""

AU_toVMm1 =5.14220632*np.power(10,11,dtype=np.double)
"""
Conversion factor from atomic units of electric field to V/m
"""

Boltzmann = 8.617333262145e-5
"""
The Boltzmann constant in eV/K  
"""

high_sym_fcc = {'G':[0.,0.,0.],
                'X':[0.,0.,1.],
                'L':[0.5,0.5,0.5],
                'K':[0.,1.,1.],
                'W':[1.0,0.5,0.]}
"""
High symmetry points of the fcc lattice (expressed in cartesian coordinates in units of 2pi/alat)
"""

high_sym_fcc_crystal = {
                'G':[0.,0.,0.],
                'X':[0.,0.5,0.5],
                'L':[0.5,0.5,0.5],
                'K':[3./8.,3./4.,3./8.],
                'W':[1./4.,3./4.,1./2.]}
"""
High symmetry points of the fcc lattice (expressed in crystal coordinates)
"""
