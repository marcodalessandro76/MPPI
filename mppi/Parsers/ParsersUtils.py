"""
Function used by the the PwParser and the YamboDftParser classes.
The elements of this module can be also directly accessed by the user to perform some
useful operations like compute the volume of a lattice, or perform the conversion from
cartesian to crystal units.
"""

import numpy as np
from mppi.Utilities.Constants import HaToeV

################################################################################
#
# def vec_in_list(veca,vec_list,atol=1e-6):
#     """
#     Check if a vector exists in a list of vectors
#     """
#     return np.array([ np.allclose(veca,vecb,rtol=atol,atol=atol) for vecb in vec_list ]).any()
#
# def expand_kpts(kpts,syms):
#     """
#     Take a list of qpoints and symmetry operations and return the full brillouin zone
#     with the corresponding index in the irreducible brillouin zone
#     """
#     full_kpts = []
#     print("nkpoints:", len(kpts))
#     for nk,k in enumerate(kpts):
#         for sym in syms:
#             full_kpts.append((nk,np.dot(sym,k)))
#
#     return full_kpts
#
# def car_red(car,lat): -> it is implented here as convert_to_crystal
#     """
#     Convert cartesian coordinates to reduced
#     """
#     return np.array([np.linalg.solve(np.array(lat).T,coord) for coord in car])
#
# def red_car(red,lat):
#     """
#     Convert reduced coordinates to cartesian
#     """
#     return np.array([coord[0]*lat[0]+coord[1]*lat[1]+coord[2]*lat[2] for coord in red])
#
# from latticedb of yambopy, note that the sym in the rec_lat is realized as the
# inverse transpose of the one in the direct lattice (check)
#
# def sym_rec(self):
#     """Convert cartesian transformations to reciprocal transformations"""
#     if not hasattr(self,"_sym_rec"):
#         sym_rec = np.zeros([self.nsym,3,3])
#         for n,s in enumerate(self.sym_car):
#             sym_rec[n] = np.linalg.inv(s).T
#         self._sym_rec = sym_rec
#     return self._sym_rec
#
# def sym_rec_red(self):
#     """Convert reduced transformations to reduced reciprocal transformations"""
#     if not hasattr(self,"_sym_rec_red"):
#         sym_rec_red = np.zeros([self.nsym,3,3],dtype=int)
#         for n,s in enumerate(self.sym_red):
#             sym_rec_red[n] = np.linalg.inv(s).T
#         self._sym_rec_red = sym_rec_red
#     return self._sym_rec_red
#
# @property
# def sym_rec(self):
#     """Convert cartesian transformations to reciprocal transformations"""
#     if not hasattr(self,"_sym_rec"):
#         sym_rec = np.zeros([self.nsym,3,3])
#         for n,s in enumerate(self.sym_car):
#             sym_rec[n] = np.linalg.inv(s).T
#         self._sym_rec = sym_rec
#     return self._sym_rec
#
# @property
# def time_rev_list(self):
#     """get a list of symmetries with time reversal"""
#     time_rev_list = [False]*self.nsym
#     for i in range(self.nsym):
#         time_rev_list[i] = ( i >= self.nsym/(self.time_rev+1) )
#     return time_rev_list
################################################################################

def get_variable_from_db(ndb_file,var_name):
    """
    Extract the value of a variable from a ndb database

    Args:
        ndb_file (:py:class:`string`) : the name of the database
        var_name (:py:class:`string`) : name of the variable

    Return:
        :py:class:`numpy.ndarray`  : array with the values of the variable
    """
    from netCDF4 import Dataset as Ds
    db = Ds(ndb_file)
    var = np.array(db.variables[var_name])
    return var

def compute_transitions(bands,in_list,fin_list):
    """
    Compute the (positive) transition energies for the bands (on a single kpoint).
    The `fast` index is associated to the bands in the fin_list list.

    Args:
        bands (list) : list with the energies
        in_list (list) : indexes of the bands used as starting points of the transitions
        fin_list (list) : indexes of the bands used as final points of the transitions

    Returns:
        :py:class:`list` : list with the transition energies for each possible couple
        of (distinct) in and out bands

    """
    transitions = []
    for v in in_list:
        for c in fin_list:
            if c > v:
                transitions.append(bands[c]-bands[v])
    return transitions

def get_evals(evals,nbands,nbands_full, set_scissor = None, set_gap = None, set_direct_gap = None, verbose = True):
    """
    Return the bands energies for each kpoint (in eV). The top of the valence band is used as the
    reference energy value. It is possible to shift the energies of the empty bands by setting an arbitrary
    value for the gap (direct or indirect) or by adding an explicit scissor. Implemented only for semiconductors.

    Args:
        evals (:py:class:`numpy.array`) : an array with the bands energy (in Hartree) for each k point
        nbands (:py:class:`int`) : number of bands
        nbands_full (:py:class:`int`) : number of occupied bands
        set_scissor (:py:class:`float`) : set the value of the scissor (in eV) that is added to the empty bands.
            If a scissor is provided the set_gap and set_direct_gap parameters are ignored
        set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system. If set_gap is provided
            the set_direct_gap parameter is ignored
        set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.

    Return:
        :py:class:`numpy.array`  : an array with the ks energies for each kpoint

    """
    evals_eV = evals*HaToeV
    homo_band = evals_eV[:,nbands_full-1]
    VBM = homo_band.max()
    evals_eV -= VBM

    if set_scissor == None and set_gap == None and set_direct_gap == None:
        return evals_eV
    elif nbands_full == nbands: # only occupied bands are present
        print('There are no empty bands. No energy shift has been applied')
        return evals_eV
    else: # shift the energy level of the empty bands if needed
        gap = get_gap(evals,nbands_full,verbose=False)
        scissor = 0.
        if set_scissor is not None:
            scissor = set_scissor
        elif set_gap is not None:
            scissor = set_gap - gap['gap']
        elif set_direct_gap is not None:
            scissor = set_direct_gap - gap['direct_gap']

        if verbose: print('Apply a scissor of',scissor,'eV')
        evals_eV[:,nbands_full:] += scissor
        return evals_eV

def get_gap(evals, nbands_full, verbose = True):
    """
    Compute the energy gap of the system (in eV). The method check if the gap is direct or
    indirect. Implemented and tested only for semiconductors.

    Args:
        evals (:py:class:`numpy.array`) : an array with the bands energy (in Hartree) for each k point
        nbands_full (:py:class:`int`) : number of occupied bands

    Return:
        :py:class:`dict` : a dictionary with the values of direct and indirect gaps and the k-point indexes
        of the VMB of the CBM and of the direct gap

    """
    val_band = evals[:,nbands_full-1]*HaToeV
    try:
        cond_band = evals[:,nbands_full]*HaToeV
    except IndexError:
        print('There are no empty states. Gap cannot be computed.')
        return None

    VBM = val_band.max()
    position_vbm = val_band.argmax()
    CBM = cond_band.min()
    position_cbm = cond_band.argmin()
    gap = (CBM-VBM)
    diff = cond_band-val_band
    position_direct_gap = np.argmin(diff)
    direct_gap = diff[position_direct_gap]

    # If there are several copies of the same point on the path it can happen that
    # the system is recognized as an indirect gap for numerical noise, so we check
    if np.allclose(gap,direct_gap,atol=1e-5,rtol=1e-5):
        gap = direct_gap
        position_cbm = position_vbm
    if verbose:
        if position_cbm == position_vbm:
            print('Direct gap system')
            print('=================')
            print('Gap :',gap,'eV')
        else :
            print('Indirect gap system')
            print('===================')
            print('Gap :',gap,'eV')
            print('Direct gap :',direct_gap,'eV')
    return {'gap':gap,'direct_gap':direct_gap,'position_cbm':position_cbm,
            'position_vbm':position_vbm,'position_direct_gap':position_direct_gap}

def get_transitions(evals, nbands, nbands_full, initial = 'full', final = 'empty',
                    set_scissor = None, set_gap = None, set_direct_gap = None):
    """
    Compute the (vertical) transitions energies. For each kpoint compute the transition energies, i.e.
    the (positive) energy difference (in eV) between the final and the initial states.

    Args:
        evals (:py:class:`numpy.array`) : an array with the bands energy (in Hartree) for each k point
        nbands (:py:class:`int`) : number of bands
        nbands_full (:py:class:`int`) : number of occupied bands
        initial (string or list) : specifies the bands from which electrons can be extracted. It can be set to `full` or
            `empty` to select the occupied or empty bands, respectively. Otherwise a list of bands can be
            provided
        final  (string or list) : specifies the final bands of the excited electrons. It can be set to `full` or
            `empty` to select the occupied or empty bands, respectively. Otherwise a list of bands can be
            provided
        set_scissor (:py:class:`float`) : set the value of the scissor (in eV) that is added to the empty bands.
            If a scissor is provided the set_gap and set_direct_gap parameters are ignored
        set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system. If set_gap is provided
            the set_direct_gap parameter is ignored
        set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.

    Return:
        :py:class:`numpy.array`  : an array with the transition energies for each kpoint

    """
    if initial == 'full':
        in_list = [ind for ind in range(nbands_full)]
    elif initial == 'empty':
        in_list = [ind for ind in range(nbands_full,nbands)]
    else:
        in_list = initial
    if final == 'full':
        fin_list = [ind for ind in range(nbands_full)]
    elif final == 'empty':
        fin_list = [ind for ind in range(nbands_full,nbands)]
    else:
        fin_list = final
    evals = get_evals(evals,nbands,nbands_full,set_scissor=set_scissor,set_gap=set_gap,set_direct_gap=set_direct_gap,verbose=False)
    transitions = []
    for bands in evals:
        transitions.append(compute_transitions(bands,in_list,fin_list))
    transitions = np.array(transitions)
    return transitions
