
TODO
----

- Check the compilation of the ReadTheDocs documentation.

- Implement a class that manage the Dos. The class should compute the dos starting from on object like PwParser.evals
  that contains the ks energies of all the bands for each kpoint.
  A method compute the dos for each of the k points: we define an energy range and an energy resolution and for each
  kpoint we build a list with the binning of the ks energies. Then all the binning lists are summed by including the
  weights of each kpoint. The class needs to manage the spin polarization.
  To manage the broadening study the DiracSuperposition class of PyBigDFT.

- Complete the reference_column_names_extendOut for the real-time runlevels.

- Test the PwParser with a result for graphene.
    1. the self.fermi return the highestOccupiedLevel why not the fermi_level???
    2. carefully check the usage of VBM and fermi in the get_gap method.
    3. study the methods setfermi in electronsdb
    4. Test the get_eval method of PwParser with a result for graphene. If the number of occupied states is k dependent
          the get eval simply return the eval in eV. Check if the band structure of graphene can be plotted.

- The weights of PwParser can be converted into a numpy.array? Test it with graphene.

- remove the xml files from IO_Files folder. 

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
