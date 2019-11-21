"""
This file contains some low-level useful functions.
The module can be loaded in the notebook in one of the following way

import mppi.Utilities.Utils as U
from mppi.Utilities import Utils as u
from mppi.Utilities.Utils import build_kpath,...
"""

def fortran_bool(boolean):
    return {True:'.true.',False:'.false.'}[boolean]

def build_kpath(*kpoints,numstep=40):
    """
    Build a list of kpoints to be passed to the set_kpoints methods of the
    :class:`PwInput` for computing the band structure along a path.

    Args:
        klist(*list): identifies the high symmetry points of the k-path
        numstep(int): specifies the number of intermediate points used to build
                      the path
    """
    klist = []
    for k in kpoints[:-1]:
        klist.append(k+[numstep])
    klist.append(kpoints[-1]+[0])
    return klist
