"""
Class to perform the parsing of the ns.db1 database in the SAVE folder of a Yambo computation.
This database collects information on the lattice properties and electronic band structure of the
system.
"""

from netCDF4 import Dataset
import numpy as np
import os

from mppi.Utilities.Constants import HaToeV
from mppi.Parsers import ParsersUtils as U
from mppi.Utilities import LatticeUtils as latUtils

class YamboDftParser():
    """
    Class to read information about the lattice and electronic structure from the ``ns.db1`` database created by Yambo

    Args:
        file (:py:class:`string`) : string with the name of the file to be parsed
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        syms : the symmetries of the lattice
        lattice : array with the lattice vectors. The i-th row represents the
            i-th lattice vector in cartesian units
        alat : the lattice parameter. Yambo stores a three dimensional array in this
            field, with the length of the cell in the three dimension
        num_electrons : number of electrons
        nbands : number of bands
        nbands_full : number of occupied bands
        nbands_empty : number of empty bands
        nkpoints : number of kpoints
        kpoints : list of the kpoints expressed in cartesian coordinates in units of 2pi/alat. Note the Yambo uses
            a vector like alat parameter, so the components of the kpoints can differ from Pw ones
        evals : array of the ks energies for each kpoint (in Hartree)
        spin : number of spin components
        spin_degen : 1 if the number of spin components is 2, 2 otherwise
        time_reversal : True if the time reversal symmetry is used by Yambo (as in YamboPy, it is read from the
            tenth element of the ``DIMENSIONS`` variable)

    The method :py:meth:`expand_IBZ_kpoints` expands the k points of the irreducible Brillouin zone to the full
    Brillouin zone and adds the attributes ``kpoints_bz``, ``kpoints_bz_crystal``, ``ibz_index``, ``sym_index``
    and ``weights``.

    """

    def __init__(self,file,verbose=True):
        self.filename = file
        if verbose: print('Parse file : %s'%self.filename)
        self.readDB()

    def readDB(self):
        """
        Read the data from the ``ns.db1`` database created by Yambo. Some variables are
        extracted from the database and stored in the attributes of the object.
        """
        try:
            database = Dataset(self.filename)
        except:
            raise IOError("Error opening file %s in YamboDftParser"%self.filename)

        # lattice properties
        self.syms = np.array(database.variables['SYMMETRY'][:])
        self.lattice = np.array(database.variables['LATTICE_VECTORS'][:].T)
        self.alat = np.array(database.variables['LATTICE_PARAMETER'][:])

        # electronic structure
        self.evals  = np.array(database.variables['EIGENVALUES'][0,:])
        self.kpoints      = np.array(database.variables['K-POINTS'][:].T)
        dimensions = database.variables['DIMENSIONS'][:]
        self.nbands      = int(dimensions[5])
        self.temperature = dimensions[13]
        self.num_electrons  = int(dimensions[14])
        self.nkpoints    = int(dimensions[6])
        self.spin = int(dimensions[11])
        self.time_reversal = bool(int(dimensions[9]))
        database.close()

        #spin degeneracy if 2 components degen 1 else degen 2
        self.spin_degen = [0,2,1][int(self.spin)]

        #number of occupied bands
        self.nbands_full = int(self.num_electrons/self.spin_degen)
        self.nbands_empty = int(self.nbands-self.nbands_full)

    def get_info(self):
        """
        Provide information on the attributes of the class
        """
        print('YamboDftParser variables structure')
        print('number of k points',self.nkpoints)
        print('number of bands',self.nbands)
        print('spin degeneration',self.spin)

    def get_evals(self, set_scissor = None, set_gap = None, set_direct_gap = None, verbose = True):
        """
        Return the ks energies for each kpoint (in eV). The top of the valence band is used as the
        reference energy value. It is possible to shift the energies of the empty bands by setting an arbitrary
        value for the gap (direct or indirect) or by adding an explicit scissor. Implemented only for semiconductors.

        Args:
            set_scissor (:py:class:`float`) : set the value of the scissor (in eV) that is added to the empty bands.
                If a scissor is provided the set_gap and set_direct_gap parameters are ignored
            set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system. If set_gap is provided
                the set_direct_gap parameter is ignored
            set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.

        Return:
            :py:class:`numpy.array`  : an array with the ks energies for each kpoint

        """
        evals = U.get_evals(self.evals,self.nbands,self.nbands_full,
                set_scissor=set_scissor,set_gap=set_gap,set_direct_gap=set_direct_gap,verbose=verbose)
        return evals

    def get_transitions(self, initial = 'full', final = 'empty',set_scissor = None, set_gap = None, set_direct_gap = None):
        """
        Compute the (vertical) transitions energies. For each kpoint compute the transition energies, i.e.
        the (positive) energy difference (in eV) between the final and the initial states.

        Args:
            initial (string or list) : specifies the bands from which electrons can be extracted. It can be set to `full` or
                `empty` to select the occupied or empty bands, respectively. Otherwise a list of bands can be
                provided
            final  (string or list) : specifies the final bands of the excited electrons. It can be set to `full` or
                `empty` to select the occupied or empty bands, respectively. Otherwise a list of bands can be
                provided
            set_scissor (:py:class:`float`) : set the value of the scissor (in eV) that is added to the empty bands.
                If a scissor is provided the set_gap and set_direct_gap parameters are ignored
            set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system. If set_gap is provided
                the set_direct_gap parameter is ignored
            set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.

        Return:
            :py:class:`numpy.array`  : an array with the transition energies for each kpoint

        """
        transitions = U.get_transitions(self.evals,self.nbands,self.nbands_full,initial=initial,final=final,
                      set_scissor=set_scissor,set_gap=set_gap,set_direct_gap=set_direct_gap)
        return transitions

    def get_gap(self, verbose = True):
        """
        Compute the energy gap of the system (in eV). The method check if the gap is direct or
        indirect. Implemented and tested only for semiconductors.

        Return:
            :py:class:`dict` : a dictionary with the values of direct and indirect gaps and the positions
            of the VMB and CBM

        """
        gap = U.get_gap(self.evals,self.nbands_full,verbose=verbose)
        return gap

    def eval_lattice_volume(self, rescale = False):
        """
        Compute the volume of the direct lattice. If ``rescale`` is False the results is expressed in a.u., otherwise
            the lattice vectors are expressed in units of alat.

        Returns:
            :py:class:`float` : lattice volume

        """
        lattice = self.get_lattice(rescale=rescale)
        return latUtils.eval_lattice_volume(lattice)

    def eval_reciprocal_lattice_volume(self, rescale = False):
        """
        Compute the volume of the reciprocal lattice. If ``rescale`` is True the reciprocal lattice vectors are expressed
        in units of 2*np.pi/alat.

        Returns:
            :py:class:`float` : reciprocal lattice volume

        """
        lattice = self.get_reciprocal_lattice(rescale=rescale)
        return latUtils.eval_lattice_volume(lattice)

    def get_lattice(self, rescale = False):
        """
        Compute the lattice vectors. If rescale = True the vectors are expressed in units
        of the lattice constant. We use the first component of the (vector) lattice constant of Yambo. Note that
        it can differ from the `alat` of the :class:`PwParser` (for a fcc cell it is alat/2), in this case
        the rescaled vectors of the two classes differ by a factor, while the vectors in a.u. are equal

        Args:
            rescale (:py:class:`bool`)  : if True express the lattice vectors in units alat

        Returns:
            :py:class:`array` : array with the lattice vectors a_i as rows

        """
        alat = self.alat[0]
        return latUtils.get_lattice(self.lattice,alat,rescale=rescale)

    def get_reciprocal_lattice(self, rescale = False):
        """
        Compute the reciprocal lattice vectors. If rescale = False the vectors are normalized
        so that np.dot(a_i,b_j) = 2*np.pi*delta_ij, where a_i is a basis vector of the direct
        lattice. If rescale = True the reciprocal lattice vectors are expressed in units of
        2*np.pi/alat. We use the first component of the (vector) lattice constant of Yambo. Note that
        it can differ from the `alat` of the :class:`PwParser` (for a fcc cell it is alat/2), in this case
        the rescaled vectors of the two classes differ by a factor, while the vectors in a.u. are equal

        Args:
            rescale (:py:class:`bool`)  : if True express the reciprocal vectors in units of 2*np.pi/alat

        Returns:
            :py:class:`array` : array with the reciprocal lattice vectors b_i as rows

        """
        alat = self.alat[0]
        return latUtils.get_reciprocal_lattice(self.lattice,alat,rescale=rescale)

    def get_kpoints(self, use_scalar_alat = True):
        """
        Get the kpoints using cartesian coordinates in units of 2*np.pi/alat (with a vector alat).

        Args:
            use_scalar_alat (:py:class:`bool`)  : if True express the kpoints in units of 2*np.pi/alat[0]

        Returns:
            :py:class:`array` : array with the kpoints

        """
        return latUtils.get_yambo_kpoints(self.kpoints,self.alat,use_scalar_alat=use_scalar_alat)


    def expand_IBZ_kpoints(self, use_time_reversal = None, atol = 1e-4, verbose = True):
        r"""
        Expand the k points of the irreducible Brillouin zone (IBZ) to the full Brillouin zone (BZ). Each IBZ point
        :math:`\mathbf{k}` is rotated by all the symmetries :math:`S` of the system (and by :math:`-S` if the time
        reversal is used) and the point :math:`S\mathbf{k}` is added to the BZ if it is not equivalent to a point
        already found, i.e. if their crystal coordinates do not differ by an integer vector within the tolerance
        ``atol``. Note that Yambo stores the k points in single precision, so the tolerance cannot be much smaller
        than the default value. The method sets the attributes:

        * ``kpoints_bz``: the BZ points :math:`S\mathbf{k}` in cartesian coordinates, in units of 2*np.pi/alat[0]
          (as :py:meth:`get_kpoints`)
        * ``kpoints_bz_crystal``: the BZ points in crystal coordinates, folded in the interval [0,1)
        * ``ibz_index``: for each BZ point the index of the IBZ point it comes from. For instance, the energies in the
          full BZ are given by ``evals[ibz_index]``
        * ``sym_index``: for each BZ point the index of the symmetry used. If the time reversal is used the
          indexes larger or equal to the number of symmetries refer to :math:`-S`, with S = syms[index-len(syms)]
        * ``weights``: the weights of the IBZ points (normalized to one)

        Args:
            use_time_reversal (:py:class:`bool`) : if True the symmetries :math:`-S` are also used. If None (default)
                the attribute ``time_reversal`` read from the database is used
            atol (:py:class:`float`) : tolerance on the crystal coordinates used to identify equivalent points
            verbose (:py:class:`bool`) : define the amount of information provided on terminal

        Returns:
            :py:class:`numpy.array` : array with the BZ points in crystal coordinates (folded in the interval [0,1))

        """
        if use_time_reversal is None: use_time_reversal = self.time_reversal
        syms = list(self.syms)
        if use_time_reversal: syms += [-s for s in self.syms]
        blat = self.get_reciprocal_lattice(rescale=True)

        self._bz_atol = atol
        self._bz_cells = {}
        kpoints_bz, kpoints_bz_crystal, ibz_index, sym_index = [], [], [], []
        for ik, k in enumerate(self.get_kpoints()):
            for isym, sym in enumerate(syms):
                k_rot = np.dot(sym,k)
                k_crys = self._fold_crystal(latUtils.convert_to_crystal(blat,k_rot))
                if self._find_bz_point(k_crys,kpoints_bz_crystal) >= 0: continue
                self._bz_cells.setdefault(self._bz_cell(k_crys),[]).append(len(kpoints_bz))
                kpoints_bz.append(k_rot)
                kpoints_bz_crystal.append(k_crys)
                ibz_index.append(ik)
                sym_index.append(isym)

        self.kpoints_bz = np.array(kpoints_bz)
        self.kpoints_bz_crystal = np.array(kpoints_bz_crystal)
        self.ibz_index = np.array(ibz_index)
        self.sym_index = np.array(sym_index)
        self.weights = np.bincount(self.ibz_index,minlength=self.nkpoints)/len(self.ibz_index)
        if verbose:
            print('Number of symmetries used: %s (time reversal: %s)'%(len(syms),use_time_reversal))
            print('%s IBZ k points expanded to %s BZ k points'%(self.nkpoints,len(self.ibz_index)))
        return self.kpoints_bz_crystal

    def get_minus_k_indexes(self):
        r"""
        For each point :math:`\mathbf{k}` of the BZ find the index of the BZ point equivalent to :math:`-\mathbf{k}`.
        It can be used to check if the k sampling is closed under the inversion, for instance when the database is
        computed without inversion and time reversal symmetries (as for a field along a given direction). The method
        :py:meth:`expand_IBZ_kpoints` has to be called first.

        Returns:
            :py:class:`numpy.array` : for each BZ point the index of the BZ point at :math:`-\mathbf{k}`, or -1 if
            :math:`-\mathbf{k}` is not in the BZ grid

        """
        if not hasattr(self,'ibz_index'):
            raise AttributeError('The k points are not expanded. Call expand_IBZ_kpoints first')
        return np.array([self._find_bz_point(self._fold_crystal(-k),self.kpoints_bz_crystal)
                         for k in self.kpoints_bz_crystal])

    def _fold_crystal(self, k_crys):
        """
        Fold the crystal coordinates in the interval [0,1). Values closer to one than the tolerance are mapped to zero.
        """
        k_crys = np.mod(k_crys,1.0)
        k_crys[k_crys > 1.0-self._bz_atol] = 0.0
        return k_crys

    def _bz_cell(self, k_crys):
        """
        Index of the cell of size atol that contains the (folded) crystal coordinates k_crys.
        """
        ncells = int(np.ceil(1.0/self._bz_atol))
        return tuple(np.floor(k_crys/self._bz_atol).astype(int) % ncells)

    def _find_bz_point(self, k_crys, kpoints_crystal):
        """
        Index of the BZ point equivalent to k_crys (crystal coordinates within the tolerance, modulo integer vectors),
        or -1 if it is not found. Only the points in the cell of k_crys and in the neighboring ones are compared.
        """
        ncells = int(np.ceil(1.0/self._bz_atol))
        cell = np.array(self._bz_cell(k_crys))
        for shift in np.ndindex(3,3,3):
            for j in self._bz_cells.get(tuple((cell+np.array(shift)-1) % ncells),[]):
                diff = kpoints_crystal[j]-k_crys
                diff -= np.round(diff)
                if np.all(np.abs(diff) < self._bz_atol): return j
        return -1
