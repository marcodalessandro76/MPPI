"""
This module contains useful physical constant and conversion factors
"""

HaToeV = 27.211386
"""
Energy conversion factor from Hartree to eV units
"""

Planck_ev_ps = 4.135667696e-3
"""
Planck constant in eV*ps.
"""

Planck_reduced_ev_ps = 6.58211957e-4
"""
Reduced Planck constant (h/2pi) in eV*ps.
"""


Bohr_radius = 5.291772e-11
"""
The Bohr radius in meter.
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

high_sym_fcc = {'G':[0.,0.,0.],
                'X':[0.,0.,1.],
                'L':[0.5,0.5,0.5],
                'K':[0.,1.,1.],
                'W':[1.0,0.5,0.]}
"""
High symmetry points of the fcc lattice (expressed in cartesian coordinates in units of 2pi/alat)
"""
