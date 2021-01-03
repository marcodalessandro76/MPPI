"""
This module contains some useful functions to perform computation with QuantumESPRESSO.
The module can be loaded in the notebook in one of the following way

>>> from mppi import Utilities as U

>>> U.build_kpath

or to load directly some elements

>>> from mppi.Utilities import build_kpath

>>> build_kpath

"""

def build_kpath(*kpoints,numstep=40):
    """
    Build a list of kpoints to be passed to the set_kpoints methods of the :class:`PwInput`
    for computing the band structure along a path.

    Example:
        >>> build_kpath(L,G,X,K,G,numstep=30)

    Args:
        kpoints : arguments that specify the high symmetry points along the k-path
        numstep (int): specifies the number of intermediate points used to build the path

    Returns:
        :py:class:`list` : list of kpoints as nedded by pw in the bands computation with the tpiba_b option

    """
    klist = []
    for k in kpoints[:-1]:
        klist.append(k+[numstep])
    klist.append(kpoints[-1]+[0])
    return klist
