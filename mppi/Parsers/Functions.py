"""
Function used by the the PwParser and the NsdbsParser classes
"""
import numpy as np
from mppi.Utilities import HaToeV

def compute_transitions(bands,in_list,fin_list):
    """
    Compute the (positive) transition energies for the bands (on a single kpoint).
    The ``fast'' index is associated to the bands in the fin_list list.

    Args:
        bands (list) : list with the energies
        in_list (list) : indexes of the bands used as starting points of the transitions
        fin_list (list) : indexes of the bands used as final points of the transitions

    Returns:
        transitions (list) : list with the transition energies for each possible couple
        of (distinct) in and out bands
    """
    transitions = []
    for v in in_list:
        for c in fin_list:
            if c > v:
                transitions.append(bands[c]-bands[v])
    return transitions

def get_evals(evals,nbands,nbands_valence, set_scissor = None, set_gap = None, set_direct_gap = None, verbose = True):
    """
    Return the bands energies for each kpoint (in eV). The top of the valence band is used as the
    reference energy value. It is possible to shift the energies of the empty bands by setting an arbitrary
    value for the gap (direct or indirect) or by adding an explicit scissor. Implemented only for semiconductors.

    Args:
        evals (:py:class:`numpy.array`) : an array with the bands energy (in Hartree) for each k point
        nbands (:py:class:`int`) : number of bands
        nbands_valence (:py:class:`int`) : number of occupied bands
        set_scissor (:py:class:`float`) : set the value of the scissor (in eV) that is added to the empty bands.
            If a scissor is provided the set_gap and set_direct_gap parameters are ignored
        set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system. If set_gap is provided
            the set_direct_gap parameter is ignored
        set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.

    Return:
        :py:class:`numpy.array`  : an array with the ks energies for each kpoint

    """
    evals_eV = evals*HaToeV
    homo_band = evals_eV[:,nbands_valence-1]
    VBM = homo_band.max()
    evals_eV -= VBM

    if set_scissor == None and set_gap == None and set_direct_gap == None:
        return evals_eV
    elif nbands_valence == nbands: # only occupied bands are present
        print('There are no empty bands. No energy shift has been applied')
        return evals_eV
    else: # shift the energy level of the empty bands if needed
        gap = get_gap(evals,nbands_valence,verbose=False)
        scissor = 0.
        if set_scissor is not None:
            scissor = set_scissor
        elif set_gap is not None:
            scissor = set_gap - gap['gap']
        elif set_direct_gap is not None:
            scissor = set_direct_gap - gap['direct_gap']

        if verbose: print('Apply a scissor of',scissor,'eV')
        evals_eV[:,nbands_valence:] += scissor
        return evals_eV

def get_gap(evals, nbands_valence, verbose = True):
    """
    Compute the energy gap of the system (in eV). The method check if the gap is direct or
    indirect. Implemented and tested only for semiconductors.

    Args:
        evals (:py:class:`numpy.array`) : an array with the bands energy (in Hartree) for each k point
        nbands_valence (:py:class:`int`) : number of occupied bands

    Return:
        :py:class:`dict` : a dictionary with the values of direct and indirect gaps and the positions
        of the VMB and CBM

    """
    homo_band = evals[:,nbands_valence-1]*HaToeV
    try:
        lumo_band = evals[:,nbands_valence]*HaToeV
    except IndexError:
        print('There are no empty states. Gap cannot be computed.')
        return None

    VBM = homo_band.max()
    position_vbm = homo_band.argmax()
    CBM = lumo_band.min()
    position_cbm = lumo_band.argmin()
    gap = (CBM-VBM)
    direct_gap = (lumo_band[position_vbm]-homo_band[position_vbm])

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
    return {'gap':gap,'direct_gap':direct_gap,'position_cbm':position_cbm,'positon_vbm':position_vbm}

def get_transitions(evals, nbands, nbands_valence, initial = 'full', final = 'empty',
                    set_scissor = None, set_gap = None, set_direct_gap = None):
    """
    Compute the (vertical) transitions energies. For each kpoint compute the transition energies, i.e.
    the (positive) energy difference (in eV) between the final and the initial states.

    Args:
        evals (:py:class:`numpy.array`) : an array with the bands energy (in Hartree) for each k point
        nbands (:py:class:`int`) : number of bands
        nbands_valence (:py:class:`int`) : number of occupied bands
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
        in_list = [ind for ind in range(nbands_valence)]
    elif initial == 'empty':
        in_list = [ind for ind in range(nbands_valence,nbands)]
    else:
        in_list = initial
    if final == 'full':
        fin_list = [ind for ind in range(nbands_valence)]
    elif final == 'empty':
        fin_list = [ind for ind in range(nbands_valence,nbands)]
    else:
        fin_list = final
    evals = get_evals(evals,nbands,nbands_valence,set_scissor=set_scissor,set_gap=set_gap,set_direct_gap=set_direct_gap,verbose=False)
    transitions = []
    for bands in evals:
        transitions.append(compute_transitions(bands,in_list,fin_list))
    transitions = np.array(transitions)
    return transitions
