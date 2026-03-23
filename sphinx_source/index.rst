.. mppi documentation master file, created by
   sphinx-quickstart on Wed Apr 24 17:01:28 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MPPI's documentation!
================================

Multi Purpose Python Interface.

Python package for managing computations and performing post-processing in QuantumESPRESSO and Yambo.

Install
-------
To install the package you can clone this repository in a local folder, e.g. /home/username/Applications/MPPI.
Then move inside the folder and install the package using the pip or the python tools,

  python3 setupy install

you can add the --user option, so that the package is installed in the local folder for the python libraries.

You can also install the package using the pip installer, in this case from the same folder run

  pip3 install -e .

Note that the *editable* -e option creates a link from the location of the package to local python repository folder.
In this way you do not need to recompile the package if you make some modifications, useful for coding.

Finally, you can also install the package inside a virtual environment.


Module Members
--------------

Manage the input files with classes able to create and modify the input files of QuantumESPRESSO
and Yambo

.. toctree::
   :maxdepth: 1

   pwInput
   phInput
   yamboInput

Run calculations in QuantumESPRESSO and Yambo. The QeCalculator manage the execution of all
the programs (pw.x,ph.x,...) of the suite.

The YamboCalculator object manage the yambo and ypp executables and also their time-dependent variants

.. toctree::
   :maxdepth: 1

   runRules
   qeCalculator
   yamboCalculator
   tools

The parsers manage the results of both QuantumESPRESSO and Yambo computations.
The parsing of the results of QuantumESPRESSO makes usage of the data-file-schema.xml and is
managed by the PwParser:

.. toctree::
   :maxdepth: 1
   :caption: Parsing of QuantumESPRESSO results

   pwParser

The parsing of the Yambo results uses both the o-* files and various databases produced by Yambo.
The collective results are managed by the YamboParser class:

.. toctree::
   :maxdepth: 1
   
   yamboParser

where the parsing of the various elements are performed by the more specific classes:

.. toctree::
   :maxdepth: 1

   yamboOutputParser
   yamboDftParser
   yamboDipolesParser
   yamboQPParser
   yamboRTCarriersParser
   yamboRTGlesserParser
   yamboNLDBParser

Finally, the module ParsersUtils collects functions for the data analysis:

.. toctree::
   :maxdepth: 1

   parsersUtils

The elements of the module are used by the parsers and can be accessed also by the user to perform
various operations

Organize runs and analyze output in a dataset

.. toctree::
   :maxdepth: 1

   datasets
   postProcessing

The *Utilities* module collects some useful tools like the management of the Dos, the study of the band structure and the generation of
gaussian pulses. Moreover the module contains useful tools to ease the computations. The module is organized in several files.

.. toctree::
   :maxdepth: 1

   attributeDict.rst
   bandStructure.rst
   constants.rst
   dos.rst
   fouriertransform.rst
   utils.rst
   latticeUtils.rst
   parallel.rst
   mergeQPndb.rst
   lrOptics.rst
   nlOptics.rst

The *Models* module collects tools to perform analysis based on some (analytical or numerical) modeling of the systems.
Actually the module contains *GaussianPulse*, a tool to deal with Gaussian shaped electromagnetic pulse and *TwoLevelSystems*.

.. toctree::
   :maxdepth: 1

   gaussianPulse.rst
   twoLevelSystems.rst


MPPI notebook section
---------------------

In this page you find some tutorials and examples that explain the usage of the package.

.. toctree::
   :maxdepth: 1

   notebooks


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
