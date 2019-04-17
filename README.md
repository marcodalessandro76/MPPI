MPPI
=======
Multi Purpose Python Interface  

Python package for managing computations and performing post-processing in QuantumESPRESSO and Yambo

Documentation
-------------
You can read the documentation in:  

http:......

Features
--------
The package is composed by several module each of wich contains one ore more classes: 

- PwIn : Create QuantumEspresso input files for pw using python (built on a modified version of the PwIn class of YamboPy compatible with python 3)  
- TODO: YamboIn : Create Yambo input files using a transparent python script  
- QeCalculator : prepare and run a QuantumESPRESSO computation (in principle all the executable of QE compatible with the launch sintax of pw can be used)  
- TODO: yamboCalculator : prepare and run a Yambo computation  

Tutorial
--------
A jupyter notebook shows the functionality of the package

Requirements
------------
- Quantum Espresso (tested with 6.3): http://www.quantum-espresso.org/
- yambo (tested with 4.4 devel): http://www.yambo-code.org/
- numpy: http://www.numpy.org/
- matplotlib: http://matplotlib.org/
- qepppy: 

Authors
------
- [Marco D'Alessandro](https://github.com/marcodalessandro76/)


The code is at an early stage of development, help us by sending bug reports, patches and suggestions!

[comment]: # - Plot the results using matplotlib
[//]: <> (This is also a comment.)

