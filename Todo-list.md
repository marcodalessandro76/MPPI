
## TODO

- Check the compilation of the ReadTheDocs documentation. It seems that the are problems in the TLS doc
  and in the associated notebook.

- Add the spin to the Dos class and add a from_Yambo method in the Dos class, this requires that the weights of the
  k points are computed by the YamboDftParser

- Complete the Analysis_Dos notebook.

- Complete the PwParser and YamboDftParser classes with add the expansion of the k points and the computation of the weigths.
  We can use the attribute weiights in the PwParser for a check of the results.

- Study the getFermi method of electronsdb of YamboPy. It can be an easy addon to the fermi method of PwParser.

- Improve the skip in the YamboCalculator to manage also the skip of the ypp computation. Is it convenient to use the
  presence of the folder $name as the criterion to skip the calculation?

- Update the tutorial on the parsing of the green function.

- Add a tutorial for the GaussianPulse class.




## FUTURE DEVELOPMENT

  1. ...
  2. ...
