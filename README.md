MPPI
=======
Multi Purpose Python Interface  

Python package for managing computations and performing post-processing in QuantumESPRESSO and Yambo

Documentation
-------------
You can read the documentation in:  

https://mppi.readthedocs.io/en/latest/

Install
-------
To install the package you can clone this repository in a local folder, e.g. /home/username/Applications/MPPI.
Then move inside the folder and install the package using the pip or the python tools,

```console
python3 setupy install
```

you can add the --user option, so that the package is installed in the local folder for the python libraries.

You can also install the package using the pip installer, in this case from the same folder run

```console
pip3 install --user -e .

```
Note that the _editable_ -e option create a link from the location of the package to local python repository folder.
In this way you do not need to recompile the package if you make some modifications, useful for coding.

Finally, you can also install the package inside a virtual environment.

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
