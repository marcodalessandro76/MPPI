
TODO
----

- Check the compilation of the ReadTheDocs documentation.

- Complete the Dos method of the PwParser to add the broadening of the dos. Maybe it is better to build a class
  as for the PwBands.

- Set the ndb of YamboCalculator equal to the name of the folder that contain the database.  

- Add the clean_run_dir method in the YamboCalculator.

- Implement a YppBands class to manage the band structure from the output of a ypp -s b calculation. 

- Check the effect of the ExtendOut option in the Yambo input.

- The get_gap method of the Pw parser can provides wrong estimate of the Direct/Indirect nature of the system to due  
  small energy differences due to numerical noise. Fix this problem (use all close to see if the direct nature is satisfied
  inside a given numerical tolerance).

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
