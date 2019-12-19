
TODO
----

- Check the compilation of the ReadTheDocs documentation.

- Implement a class that manage the Dos. The class should compute the dos starting from on object like PwParser.evals
  that contains the ks energies of all the bands for each kpoint.
  A method compute the dos for each of the k points: we define an energy range and an energy resolution and for each
  kpoint we build a list with the binning of the ks energies. Then all the binning lists are summed by including the
  weights of each kpoint. The class needs to manage the spin polarization.
  To manage the broadening study the DiracSuperposition class of PyBigDFT.

- Implement a YppBands class to manage the band structure from the output of a ypp -s b calculation.

- Complete the reference_column_names_extendOut for the real-time runlevels.

- Use the function floats_from_string in YamboParser and clean the eval_columns function.

- Update of PwParser
  1. check if it is better to to express the evals in eV and/or define a method that : convert the eval in eV (if requested),
     set the fermi energy to zero and translate the conduction bands setting the gap to a given value, and then return the
     evals

- Update of PwBands : this class can become more general purpose. To do so we have to define a structure for the ks energies
  that the class is able to manage. Then we define some methods that convert the output of Pw or Ypp calculation in this format,
  so the class can manage the BandStructure for various type of inputs. 

FUTURE DEVELOPMENT
------------------

  1. Add the possibility to run many computations in parallel (with and/or without a job manager like slurm).
     Possible strategies:
     a. (without the job manager) : add to calculator the possibility of taking a list of input objects as input
     (or define a new class with this feature). Once that this is done must run all the input in parallel
     (use multiprocessing module of python, otherwise os.system(command + '&') can work?). Then we need to modify dataset
     adding an 'is_parallel' boolean attribute, if it is true all the elements of the dataset are passed to the calculator together.
     b. (with the job manager) : defined a calculator class the write the proper script and submit it. This calculator
     has to be able to take a list of input objects as input, so the implementation of this feature is a preliminary step
     common to both the approaches. In this case dataset uses this calculator passing all its elements together.

  2.
