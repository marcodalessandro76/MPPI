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
            field but here only the first element is parsed to have to same structure
            of the :class:`PwParser` class
        num_electrons : number of electrons
        nbands : number of bands
        nbands_full : number of occupied bands
        nbands_empty : number of empty bands
        nkpoints : numer of kpoints
        kpoints : list of the kpoints expressed in cartesian coordinates in units of 2pi/alat. Note the Yambo uses
            a vector like alat parameter, so the components of the kpoints can differ from Pw ones
        evals : array of the ks energies for each kpoint (in Hartree)
        spin : number of spin components
        spin_degen : 1 if the number of spin components is 2, 2 otherwise

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
        self.alat = database.variables['LATTICE_PARAMETER'][:][0]

        # electronic structure
        self.evals  = np.array(database.variables['EIGENVALUES'][0,:])
        self.kpoints      = np.array(database.variables['K-POINTS'][:].T)
        dimensions = database.variables['DIMENSIONS'][:]
        self.nbands      = int(dimensions[5])
        self.temperature = dimensions[13]
        self.num_electrons  = int(dimensions[14])
        self.nkpoints    = int(dimensions[6])
        self.spin = int(dimensions[11])
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
        of the lattice constant.

        Args:
            rescale (:py:class:`bool`)  : if True express the lattice vectors in units alat

        Returns:
            :py:class:`array` : array with the lattice vectors a_i as rows

        """
        return latUtils.get_lattice(self.lattice,self.alat,rescale=rescale)

    def get_reciprocal_lattice(self, rescale = False):
        """
        Compute the reciprocal lattice vectors. If rescale = False the vectors are normalized
        so that np.dot(a_i,b_j) = 2*np.pi*delta_ij, where a_i is a basis vector of the direct
        lattice. If rescale = True the reciprocal lattice vectors are expressed in units of
        2*np.pi/alat.

        Args:
            rescale (:py:class:`bool`)  : if True express the reciprocal vectors in units of 2*np.pi/alat

        Returns:
            :py:class:`array` : array with the reciprocal lattice vectors b_i as rows

        """
        return latUtils.get_reciprocal_lattice(self.lattice,self.alat,rescale=rescale)


    #####################################################################################
    #
    # def expand_kpoints(self,atol=1e-6,verbose=0):
    #     """
    #     Take a list of qpoints and symmetry operations and return the full brillouin zone
    #     with the corresponding index in the irreducible brillouin zone
    #     """
    #
    #     #check if the kpoints were already exapnded
    #     kpoints_indexes  = []
    #     kpoints_full     = []
    #     symmetry_indexes = []
    #
    #     #kpoints in the full brillouin zone organized per index
    #     kpoints_full_i = {}
    #
    #     #expand using symmetries
    #     for nk,k in enumerate(self.car_kpoints):
    #         #if the index in not in the dicitonary add a list
    #         if nk not in kpoints_full_i:
    #             kpoints_full_i[nk] = []
    #
    #         for ns,sym in enumerate(self.sym_car):
    #
    #             new_k = np.dot(sym,k)
    #
    #             #check if the point is inside the bounds
    #             k_red = car_red([new_k],self.rlat)[0]
    #             k_bz = (k_red+atol)%1
    #
    #             #if the vector is not in the list of this index add it
    #             if not vec_in_list(k_bz,kpoints_full_i[nk]):
    #                 kpoints_full_i[nk].append(k_bz)
    #                 kpoints_full.append(new_k)
    #                 kpoints_indexes.append(nk)
    #                 symmetry_indexes.append(ns)
    #                 continue
    #
    #     #calculate the weights of each of the kpoints in the irreducible brillouin zone
    #     nkpoints_full = len(kpoints_full)
    #     weights = np.zeros([nkpoints_full])
    #     for nk in kpoints_full_i:
    #         weights[nk] = float(len(kpoints_full_i[nk]))/nkpoints_full
    #
    #     if verbose: print("%d kpoints expanded to %d"%(len(self.car_kpoints),len(kpoints_full)))
    #
    #     #set the variables
    #     self.weights_ibz      = np.array(weights)
    #     self.kpoints_indexes  = np.array(kpoints_indexes)
    #     self.symmetry_indexes = np.array(symmetry_indexes)
    #     self.iku_kpoints      = [k*self.alat for k in kpoints_full]
    #
    # def expandEigenvalues(self):
    #     """
    #     Expand eigenvalues to the full brillouin zone
    #     """
    #
    #     self.eigenvalues = self.eigenvalues_ibz[self.lattice.kpoints_indexes]
    #
    #     self.nkpoints_ibz = len(self.eigenvalues_ibz)
    #     self.weights_ibz = np.zeros([self.nkpoints_ibz],dtype=np.float32)
    #     self.nkpoints = len(self.eigenvalues)
    #
    #     #counter counts the number of occurences of element in a list
    #     for nk_ibz,inv_weight in list(collections.Counter(self.lattice.kpoints_indexes).items()):
    #         self.weights_ibz[nk_ibz] = float(inv_weight)/self.nkpoints
    #
    #     #kpoints weights
    #     self.weights = np.full((self.nkpoints), 1.0/self.nkpoints,dtype=np.float32)
