
## TODO

- Check the compilation of the ReadTheDocs documentation. The are problems for the rendering of the inline math equations
  in the TLS notebook.

- Add the spin to the Dos class and add a from_Yambo method in the Dos class, this requires that the weights of the
  k points are computed by the YamboDftParser

- Complete the Analysis_Dos notebook.

- Complete the PwParser and YamboDftParser classes with add the expansion of the k points and the computation of the weigths.
  We can use the attribute weights in the PwParser for a check of the results.

- Study the getFermi method of electronsdb of YamboPy. It can be an easy addon to the fermi method of PwParser.

- Update the tutorial on the parsing of the green function.

- Add a tutorial for the GaussianPulse class.

- Improve the run_the_calculations method of Dataset using the same approach introduced in the loop function of the Parallel
  module. The loop that wait the end of the processes can be removed and the extraction of data from the Queue
  can be changed.

- Yambo uses a vector like alat parameter, whereas Pw uses a scalar alat. Maybe we need to modify the PwParser class to deal with
  vector alat (with equal components). In this way the function implemented to work with the Parser can be used by both the Yambo
  and Pw parsers. Actually only the first element of the Yambo alat is stored in the YamboDftParser alat parameter.

- The YamboCalculator does not work properly when there are several replica of the report and o- files. In this case we
  need to improve the get_output_files and get_report functions to be able to identify only the last elements that are
  associated to the new run. We can use the _number label to identify the replica for each type of o- files.

- The job_*.out file in the YamboCalculator has to be erased also when the clean_restart option is False. Otherwise
  when the calculation starts the calculator can erroneously assume that the computation is ended. Same behavior also
  for the QeCalculator?

- Check the Analysis_BandStructure tutorial and control the ypp band (since -a 2 option has been removed in p2y default).

- Define an update_from_remote function that implement the usage of rsync to fetch the results computed in a remote folder
  into a local one. The function should implement a command like:

  rsync -rLptgoDzv --exclude={'*_fragment_*','*_fragments_*'} -e ssh ismhpc:/remote_path local_path

  the options --update and --dry-run can be included


## FUTURE DEVELOPMENT

  1. ...
  2. ...
