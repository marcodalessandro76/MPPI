
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

- The command make html produce the warning
  WARNING: mathjax_config/mathjax2_config does not work for the current MathJax version, use mathjax3_config instead
  Check!

- Check the Analysis_BandStructure tutorial and control the ypp band (since -a 2 option has been removed in p2y default).

- Define an update_from_remote function that implement the usage of rsync to fetch the results computed in a remote folder
  into a local one. The function should implement a command like:

  rsync -rLptgoDzv --exclude={'*_fragment_*','*_fragments_*'} -e ssh ismhpc:/remote_path local_path

  the options --update and --dry-run can be included


## FUTURE DEVELOPMENT

  1. ...
  2. ...
