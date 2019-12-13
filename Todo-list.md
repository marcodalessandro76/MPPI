
TODO
----

- In a nscf computation with skip=true the computation is always skipped because the system recognize the
  data-file-schema.xls of the source scf computation. This needs to be fixed. I can delete this file after
  the copy of the scf source folder or I can use the name.xml file as a check to understand if the computation
  has been already performed.
  It seems that the data-file-schema has to be present to perform the nscf computation, I try to use the
  $name.xml file as the skipfile. __DO A CAREFUL CHECK OF THE PROCEDURE__.

- Complete the Dos method of the PwParser to add the broadening of the dos. Maybe it is better to build a class
  as for the PwBands.

- Add a function that compute the gap (both direct and indirect) from the parsing of a pw computation. I can create
  PwUtils and collect in this module function like this. Maybe the get_gap function can be inserted directly in the
  PwParser. Also a function get_valence_band can be defined (with a check that the band is the same for all the kpoints)
  I can use something like if all(v == values[0] for v in values):


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
