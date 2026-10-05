
## ISSUES

- Check the compilation of the ReadTheDocs documentation. The are problems for the rendering of the inline math equations
  in the TLS notebook.


## TODO

- Complete the YamboQPParser class and the associated tutorial in the YamboParser notebook. Use the YamboQPParser to complete
  the MergeQPndb function in the Utilities module.

- Add the spin to the Dos class and add a from_Yambo method in the Dos class, this requires that the weights of the
  k points are computed by the YamboDftParser.

- Complete the PwParser and YamboDftParser classes with add the expansion of the k points and the computation of the weigths.
  We can use the attribute weights in the PwParser for a check of the results.
  For YamboDftParser restore the commented `expand_kpoints`/`expandEigenvalues` methods (they use old names: `car_kpoints`,
  `sym_car`, `rlat`, `car_red`, `vec_in_list`) with the current attributes: `self.syms` (cartesian, (nsym,3,3)),
  `get_kpoints()` and `get_reciprocal_lattice(rescale=True)` (both in units of 2pi/alat), `LatticeUtils.convert_to_crystal`.
  Store the BZ points (crystal and cartesian), the IBZ index and the symmetry index of each BZ point and the IBZ weights;
  two points are equal if their crystal coordinates differ by an integer vector (tol ~1e-5). Time reversal is not in
  `self.syms`: add an option `use_time_reversal` (default False). Add a helper that gives, for each BZ point, the index of
  the point at -k, to check the inversion symmetry of the sampling of a SAVE without inversion and time reversal (e.g. the
  fixsym/NoTr SAVEs used for a field along a given direction).
  Use case and reference implementation (explicit loop on IBZ points and symmetries): LiF project, notebook
  `NL_Chi/YamboNL_Analysis.ipynb`, section "Ground state: inversion symmetry of the k sampling" (fcc, 8 symmetries,
  8x8x8 Gamma-centered grid: 100 IBZ points -> 512 BZ points, closed under k -> -k). Add a unit test (fcc lattice with a
  few symmetries, number of BZ points and closure under k -> -k).

- Study the getFermi method of electronsdb of YamboPy. It can be an easy addon to the fermi method of PwParser.

- Update the tutorial on the parsing of the green function.

- Add a tutorial for the GaussianPulse class.

- Improve the run_the_calculations method of Dataset using the same approach introduced in the loop function of the Parallel
  module. The loop that wait the end of the processes can be removed and the extraction of data from the Queue
  can be changed.

- Define an update_from_remote function that implement the usage of rsync to fetch the results computed in a remote folder
  into a local one. The function should implement a command like:

  rsync -rLptgoDzv --exclude={'*_fragment_*','*_fragments_*'} -e ssh ismhpc:/remote_path local_path

  the options --update and --dry-run can be included
