"""
This module contains useful functions for working with lattice related topics
like compute the volume of a (direct or reciprocal) lattice, perform the conversion from
cartesian to crystal units or generate the atomic positions for a given value of replica
of the lattice cell. This functions are also used by the PwParser and the YamboDftParser classes.
"""

import numpy as np
from itertools import product

def convert_to_crystal(lattice,vector_cartesian):
    """
    Convert the cartesian coordinates of a vector into crystal ones. It can be used
    for both direct lattice vectors (like the atomic position) and reciprocal lattice
    related quantities (like the k-points).

    Args:
        lattice (:py:class:`array`) : array with the (direct or reciprocal) lattice vectors.
            The i-th row represents the i-th lattice vectors in cartesian units
        vector_cartesian (:py:class:`array`) : array with the cartesian coordinates
            of the vector

    Returns:
        (:py:class:`array`) : array with the crystal coordinates of the vector

    """
    M = np.linalg.inv(lattice.transpose())
    return np.dot(M,vector_cartesian)

def convert_to_cartesian(lattice,vector_crystal):
    """
    Convert the crystal coordinates of a vector into cartesian ones. It can be used
    for both direct lattice vectors (like the atomic position) and reciprocal lattice
    related quantities (like the k-points).

    Args:
        lattice (:py:class:`array`) : array with the (direct or reciprocal) lattice vectors.
            The i-th row represents the i-th lattice vectors in cartesian units
        vector_crystal (:py:class:`array`) : array with the crystal coordinates
            of the vector

    Returns:
        (:py:class:`array`) : array with the cartesian coordinates of the vector

    """
    M = lattice.transpose()
    return np.dot(M,vector_cartesian)

def eval_lattice_volume(lattice):
    """
    Compute the volume of a lattice

    Args:
        lattice (:py:class:`array`) : array with the lattice vectors. The i-th
            row represents the i-th lattice vectors in cartesian units

    Returns:
        :py:class:`float` : lattice volume

    """
    a1,a2,a3 = np.array(lattice)
    return np.dot(a1,np.cross(a2,a3))

def get_lattice(lattice, alat, rescale = False):
    """
    Compute the lattice vectors. If rescale = True the vectors are expressed in units
    of the lattice constant.

    Args:
        lattice (:py:class:`array`) : array with the lattice vectors. The i-th
            row represents the i-th lattice vectors in cartesian units
        alat (:py:class:`float`) : lattice parameter
        rescale (:py:class:`bool`)  : if True express the lattice vectors in units alat

    Returns:
        :py:class:`array` : array with the lattice vectors a_i as rows

    """
    if not rescale:
        return lattice
    else:
        return lattice/alat

def get_reciprocal_lattice(lattice, alat, rescale = False):
    """
    Calculate the reciprocal lattice vectors. If rescale = False the vectors are normalized
    so that np.dot(a_i,b_j) = 2*np.pi*delta_ij, where a_i is a basis vector of the direct
    lattice. If rescale = True the reciprocal lattice vectors are expressed in units of
    2*np.pi/alat

    Args:
        lattice (:py:class:`array`) : array with the lattice vectors. The i-th
            row represents the i-th lattice vectors in cartesian units
        alat (:py:class:`float`) : lattice parameter
        rescale (:py:class:`bool`)  : if True express the reciprocal lattice vectors in units of 2*np.pi/alat

    Returns:
        :py:class:`array` : array with the reciprocal lattice vectors b_i as rows

    """
    a1,a2,a3 = lattice
    vol = eval_lattice_volume(lattice)
    if rescale: norm_factor = alat/vol
    else: norm_factor = 2.*np.pi/vol
    b1 = norm_factor*np.cross(a2,a3)
    b2 = norm_factor*np.cross(a3,a1)
    b3 = norm_factor*np.cross(a1,a2)
    return np.array([b1,b2,b3])

def unit_cell_3d(lattice, atom_pos, Nx, Ny, Nz):
    """
    Make arrays of x-, y- and z-positions of a lattice from the
    lattice vectors, the atom positions and the number of unit cells.

    Args:
        lattice (:py:class:`array`) : array with the lattice vectors. The i-th
            row represents the i-th lattice vectors in cartesian units
        atom_pos (:py:class:`list`) : each element of the list represents the position
            of an atom in the unit cell, written in the basis of the lattice vectors
        Nx (:py:class:`int`) : number of unit cells to be plotted in the x-direction
        Ny (:py:class:`int`) : number of unit cells to be plotted in the y-direction
        Nz (:py:class:`int`) : number of unit cells to be plotted in the z-direction

    Returns:
        :py:class:`array` : Array containing the coordinates of all atoms to be plotted
    """
    a,b,c = lattice
    latt_coord_x = []
    latt_coord_y = []
    latt_coord_z = []

    for atom in atom_pos:
        xpos = atom[0]*a[0] + atom[1]*b[0] + atom[2]*c[0]
        ypos = atom[0]*a[1] + atom[1]*b[1] + atom[2]*c[1]
        zpos = atom[0]*a[2] + atom[1]*b[2] + atom[2]*c[2]
        xpos_all = [xpos + n*a[0] + m*b[0] + k*c[0] for n, m, k in
                     product(range(Nx), range(Ny), range(Nz))]
        ypos_all = [ypos + n*a[1] + m*b[1] + k*c[1] for n, m, k in
                     product(range(Nx), range(Ny), range(Nz))]
        zpos_all = [zpos + n*a[2] + m*b[2] + k*c[2] for n, m, k in
                     product(range(Nx), range(Ny), range(Nz))]
        latt_coord_x.append(xpos_all)
        latt_coord_y.append(ypos_all)
        latt_coord_z.append(zpos_all)
    latt_coord_x = np.array(latt_coord_x).flatten()
    latt_coord_y = np.array(latt_coord_y).flatten()
    latt_coord_z = np.array(latt_coord_z).flatten()
    return latt_coord_x, latt_coord_y, latt_coord_z
