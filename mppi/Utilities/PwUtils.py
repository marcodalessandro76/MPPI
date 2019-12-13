"""
This file contains some useful functions to perform computation with QuantumESPRESSO.
The module can be loaded in the notebook in one of the following way

>>> from mppi import Utilities as U

>>> U.eval_gap

or to load directly some elements

>>> from mppi.Utilities import eval_gap

>>> eval_gap

"""

def get_gap(results):
    """
    Compute the energy gap from a QuantumESPRESSO nscf computation.
    The function make usage of the PwParser of this package.
    """
    from mppi import Parsers as P
    data = P.PwParser(results,verbose=False)
    return 0.
