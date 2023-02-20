"""
This module contains useful functions for working with lattice related topics
like compute the volume of a (direct or reciprocal) lattice, perform the conversion from
cartesian to crystal units or generate the atomic positions for a given value of replica
of the lattice cell. This functions are also used by the PwParser and the YamboDftParser classes.
"""

import numpy as np

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
    return np.dot(M,vector_crystal)

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

def get_yambo_kpoints(kpoints, alat, use_scalar_alat = True):
    """
    Get the kpoints using cartesian coordinates in units of 2*np.pi/alat (with a vector alat).

    Args:
        kpoints (:py:class:`array`) : array with the kpoints.
        alat (:py:class:`array`) : vector lattice parameters
        use_scalar_alat (:py:class:`bool`)  : if True express the kpoints in units of 2*np.pi/alat[0]

    Returns:
        :py:class:`array` : array with the kpoints

    """
    if use_scalar_alat: kpoints = alat[0]/alat*kpoints
    return kpoints

def build_lattice(lattice_vectors, atom_pos, Nx, Ny, Nz):
    """
    Create the arrays of the x,y and z positions of each atom in the lattice.
    The function make usage of the lattice vectors, of the atoms positions and
    of the number of replica of the unit cell in the x,y and z directions.

    Args:
        lattice_vectors (:py:class:`array`) : array with the lattice vectors. The i-th
            row represents the i-th lattice vectors in cartesian units
        atom_pos (:py:class:`list`) : each element has the structure of the
            atomic_positions attribute of the pw xml data-file-schema
            ['atom_name',x,y,z] and specify the position of the atom in the unit cell
            in cartesian coordinate
        Nx (:py:class:`int`) : number of unit cells to be plotted in the x-direction
        Ny (:py:class:`int`) : number of unit cells to be plotted in the y-direction
        Nz (:py:class:`int`) : number of unit cells to be plotted in the z-direction

    Returns:
        :py:class:`list` : Each element has the structure
            ['atom_name',x,y,z] where x,y and z are arrays with corresponding coordinates
            of all the atoms of 'atom__name' type in the lattice

    """
    from itertools import product

    a,b,c = lattice_vectors
    latt_coord = []
    for atom in atom_pos:
        atom_name = atom[0]
        position = atom[1]
        xpos = np.array([position[0] + n*a[0] + m*b[0] + k*c[0] for n, m, k in
                     product(range(Nx), range(Ny), range(Nz))])
        ypos = np.array([position[1] + n*a[1] + m*b[1] + k*c[1] for n, m, k in
                     product(range(Nx), range(Ny), range(Nz))])
        zpos = np.array([position[2] + n*a[2] + m*b[2] + k*c[2] for n, m, k in
                     product(range(Nx), range(Ny), range(Nz))])
        latt_coord.append([atom_name,xpos,ypos,zpos])
    return latt_coord
