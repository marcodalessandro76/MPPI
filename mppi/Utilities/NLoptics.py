"""
This module defines the tools to extract the non-linear optical properties of the system from the real-time polarization.
The module can be loaded in the notebook as follows

>>> from mppi.Utilities import NLoptics as NL

>>> LR.Nonlinear_Response

"""
import numpy as np
from mppi.Utilities import FourierTransform as FT
from mppi.Utilities import Constants as C
from mppi.Utilities.Utils import damp_ft

