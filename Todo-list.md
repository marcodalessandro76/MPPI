
TODO
----

- Check the compilation of the ReadTheDocs documentation.

- Add the spin to the Dos class.

- Complete the Analysis_Dos notebook.

- Add a from_Yambo method in the Dos class. The output of the YamboParser does not contain
  information on the weights of the kpoints. These information can be extracted from ns.db1
  using the YamboElectronsDB class of yambopy or from the xml file of the pw output used to
  build the SAVE folder.

- Study the getFermi method of electronsdb of YamboPy. It can be an easy addon to the fermi method of PwParser

- Study the expand_kpoints method of latticedb and the expandEigenvalues method of electronsdb of YamboPy. This
  methods are useful to define a weigth of each k point in the NsdbsParser class.

- Change the names nbands_conduction -> nbands_empty and nbands_valence -> nbands_full.

- Improve the skip in the YamboCalculator to manage also the skip of the ypp computation

- Complete the tutorial on the YamboParser of the output files.  

- Extend and/or modify the get_transitions methods of PwParser. We need to build weight array associated to the transitions
  that is needed for computing the JDOS (if non equal weights are associated to all the kpoints) and also to build the
  absorption spectrum.

FUTURE DEVELOPMENT
------------------

  1. ...
  2. ...
