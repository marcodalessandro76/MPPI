# MPPI
__Multi Purpose Python Interface__  

Python package for managing computations and performing post-processing in QuantumESPRESSO and Yambo

## Documentation
You can read the documentation in:  

https://mppi.readthedocs.io/en/latest/

## Features
The package is composed by several module each of which contains one ore more classes:

- __InputFiles__ : create and manage input files for pw.x program of the QuantumEspresso package and for Yambo.
- __Calculators__ : prepare and run a single QuantumESPRESSO or Yambo computation.
- __Datasets__ : organize and run several computations for both QuantumESPRESSO and Yambo.
- __Parsers__ : classes to extract data from the output files and database QuantumESPRESSO and Yambo.
- __Utilities__ : collect some useful low-level functions used by the other modules of the package.
- __Models__ : implement some physical models. Actually the module contains one tool to deal with gaussian pulses and
               and one to analyze the time dynamics of a two-level-system.

## Installation
To install the package you can clone this repository in a local folder, e.g. `/home/username/Applications/MPPI`.
Then move inside the folder and install the package using the pip or the python tools,

    python3 setupy install

you can add the --user option, so that the package is installed in the local folder for the python libraries.

You can also install the package using the pip installer, in this case from the same folder run

    pip3 install -e .

Note that the _editable_ -e option creates a link from the location of the package to local python repository folder.
In this way you do not need to recompile the package if you make some modifications, useful for coding.

Finally, you can also install the package inside a virtual environment.

It is possible to automatically generate the documentation. To this scope run

    make html

from the folder where the package is installed. The root of the documentation is located in the file
`/package_dir/sphinx_build/html/index.html`

## Tutorials and examples
We provide many jupyter notebooks that show the functionality of each module of the package.
The tutorials are organized as `Tutorial_$name_of_the_class`. Furthermore we also include a detailed
example of Band structure calculation for various systems, this topic in included in the
`Analysis_BandStructure` notebook.

To run the notebooks you need to install the jupyter-notebook or jupyter-lab package, that can be installed as

    pip3 install jupyterlab

then from the source folder of MPPI move to `sphinx_source/tutorials` and run the jupyter environment as

    jupyter-lab &

## Requirements
- QuantumESPRESSO (tested with 6.6): http://www.quantum-espresso.org/
- yambo (tested with 5.0 devel): http://www.yambo-code.org/

## Authors
- [Marco D'Alessandro](https://github.com/marcodalessandro76/)

The code is at an early stage of development, help us by sending bug reports and suggestions!
