
TODO
----

- Check the compilation of the ReadTheDocs documentation.

- Add the spin to the Dos class.

- Complete the documentation of the Dos class.

- Add a from_Yambo method in the Dos class. The output of the YamboParser does not contain
  information on the weights of the kpoints. These information can be extracted from ns.db1
  using the YamboElectronsDB class of yambopy or from the xml file of the pw output used to
  build the SAVE folder.

- Define a get_transitions method in the PwParser.

- Complete the reference_column_names_extendOut for the real-time runlevels.

- Test the PwParser with a result for graphene.
    1. why there no self.fermi?

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
