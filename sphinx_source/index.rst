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

Class to create and modify the pw input file of QuantumESPRESSO

.. toctree:: pwInput

Class to create and modify the Yambo input file

.. toctree:: yamboInput
.. toctree:: yppIn

Calculator for QuantumESPRESSO. This Class can manage the execution of all the
programs (pw.x,ph.x,...) of the suite

.. toctree:: qeCalculator

Calculator for Yambo.

.. toctree:: yamboCalculator

Parser to manage the results of the QuantumESPRESSO and Yambo computations

.. toctree:: pwParser
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
