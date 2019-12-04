
TODO
----

- Modify the build_SAVE function of Utils so that a link of the SAVE file (and not a copy) is created.

- The QeCalculator skip the computation looking for the file self.run_options['name'].xml. It is better
  directly look for data-file-schema.xml which is the file used for parsing.

- In the QeCalculator do not set the results equal to data-file-schema.xml if the computation is not performed.
  Check if the data-file-schema.xml exists

- Complete the Dos method of the PwParser to add the broadening of the dos.    


FUTURE DEVELOPMENT
------------------

  1. Add the possibility to run many computations in parallel (with and/or without a job manager like slurm).
     Possible strategies:
     a. (without the job manager) : add to calculator the possibility of taking a list of input objects as input
     (or define a new class with this feature). Once that this is done must run all the input in parallel
     (os.system(command + '&') can work?). Then we need to modify dataset adding an 'is_parallel' boolean attribute,
     if it is true all the elements of the dataset are passed to the calculator together.
     b. (with the job manager) : defined a calculator class the write the proper script and submit it. This calculator
     has to be able to take a list of input objects as input, so the implementation of this feature is a preliminary step
     common to both the approaches. In this case dataset uses this calculator passing all its elements together.

  2.
