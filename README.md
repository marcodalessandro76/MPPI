MPPI
=======
Multi Purpose Python Interface  

Python package for managing computations and performing post-processing in QuantumESPRESSO and Yambo

Documentation
-------------
You can read the documentation in:  

https://mppi.readthedocs.io/en/devel/

Features
--------
The package is composed by several module each of which contains one ore more classes:

- InputFiles : Create and manage input files for pw.x program of the QuantumEspresso package and for Yambo.
- Calculators : prepare and run a single QuantumESPRESSO or Yambo computation.
- Datasets : organize and run several computations for both QuantumESPRESSO and Yambo.
- Parsers : Classes to extract data from the output files and database QuantumESPRESSO and Yambo.
    Other parsers, from outer packages, can be also used.
- Utilities : Collect some useful low-level functions used by the other modules of the package.

Tutorials and examples
----------------------
We provide some jupyter notebooks to show the functionality of the package.

Requirements
------------
- Quantum Espresso (tested with 6.3): http://www.quantum-espresso.org/
- yambo (tested with 4.4 devel): http://www.yambo-code.org/
- numpy: http://www.numpy.org/

Authors
------
- [Marco D'Alessandro](https://github.com/marcodalessandro76/)

The code is at an early stage of development, help us by sending bug reports and suggestions!
