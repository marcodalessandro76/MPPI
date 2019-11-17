.. mppi documentation master file, created by
   sphinx-quickstart on Wed Apr 24 17:01:28 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MPPI's documentation!
================================

Multi Purpose Python Interface.

Python package for managing computations and performing post-processing in QuantumESPRESSO and Yambo.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

Module Members
--------------

QuantumESPRESSO input file representation and usage

.. toctree:: qeInput

Yambo input file representation and usage

.. toctree:: yamboInput
.. toctree:: yppIn

Calculator for QuantumESPRESSO

.. toctree:: qeCalc

Calculator for Yambo. The module contain also a parses for the output of Yambo computation

.. toctree:: yamboCalc
.. toctree:: yamboParser

Organize runs and analyze output in a dataset. A set of pre-processing functions
needed before launching the various type of datasets are given in the PreProcessings
module

.. toctree:: datasets
.. toctree:: preProcessings

MPPI Tutorial page

.. toctree::
   :maxdepth: 1

   tutorials


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
