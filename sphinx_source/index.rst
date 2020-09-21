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

Note that the _editable_ -e option creates a link from the location of the package to local python repository folder.
In this way you do not need to recompile the package if you make some modifications, useful for coding.

Finally, you can also install the package inside a virtual environment.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

Module Members
--------------

Manage the input files with classes able to create and modify the input files of QuantumESPRESSO
and Yambo

.. toctree:: pwInput
.. toctree:: yamboInput

Run (many) calculations in QuantumESPRESSO and Yambo. The QeCalculator manage the execution of all
the programs (pw.x,ph.x,...) of the suite.

The YamboCalculator object manage the yambo and ypp executables and also their time-dependent variants

.. toctree:: qeCalculator
.. toctree:: yamboCalculator

Parsers to manage the results of the QuantumESPRESSO and Yambo computations

.. toctree:: pwParser
.. toctree:: yamboParser

Organize runs and analyze output in a dataset

.. toctree:: datasets
.. toctree:: postProcessing

The Utilities module collects some useful functions and classes. The module is organized in several files.

.. toctree:: constants.rst
.. toctree:: attributeDict.rst
.. toctree:: pwUtils.rst
.. toctree:: yamboUtils.rst
.. toctree:: bandStructure.rst
.. toctree:: dos.rst

MPPI Tutorials page
-------------------

In this page you find some tutorials and examples that explain the usage of the package.

.. toctree::
   :maxdepth: 1

   tutorials


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
