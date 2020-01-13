
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

- Update Analysis_BandStructure due to the modifications of the QeCalculator.

- Test of Dataset with the new YamboCalculator


FUTURE DEVELOPMENT
------------------

  1. ...
  2. ...
