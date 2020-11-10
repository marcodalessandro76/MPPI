"""
Class to perform the parsing of the ns.db1 database in the SAVE folder of a Yambo computation.
This database collects information on the lattice properties and electronic band structure of the
system.
"""

from netCDF4 import Dataset
import numpy as np
import os

from mppi.Utilities import HaToeV
from mppi.Parsers import Functions as F


class NsdbsParser():
    """
    Class to read information about the lattice and electronic structure from the ``ns.db1`` database created by Yambo

    Args:
        save(:py:class:`string`) : Yambo SAVE folder that contains the database
        nbs(:py:class:`string`) : name of the electrons database

    Attributes:
        sym : the symmetries of the lattice
        lat : the lattice vectors in cart. coord. in a.u.
        alat : the lattice parameter(s) in a.u.

        num_electrons : number of electrons
        nbands : number of bands
        nbands_valence : number of occupied bands
        nbands_conduction : number of empty bands
        nkpoints : numer of kpoints
        kpoints : list of the kpoints
        evals : array of the ks energies for each kpoint (in Hartree)
        spin : number of spin components
        spin_degen : 1 if the number of spin components is 2, 2 otherwise
        time_rev : ?
    """

    def __init__(self,save='SAVE',dbs='ns.db1',verbose=True):
        self.filename = os.path.join(save,dbs)
        if verbose: print('Parse file : %s'%self.filename)
        self.readDB()

    def readDB(self):
        """
        Read the data from the NSCF database created by Yambo. Some variables are
        extracted from the database and stored in the attributes of the object.
        """
        try:
            database = Dataset(self.filename)
        except:
            raise IOError("Error opening file %s in NsdbParser"%self.filename)

        # lattice properties
        self.sym = np.array(database.variables['SYMMETRY'][:])
        self.lat = np.array(database.variables['LATTICE_VECTORS'][:].T)
        self.alat = np.array(database.variables['LATTICE_PARAMETER'][:].T)

        # electronic structure
        self.evals  = np.array(database.variables['EIGENVALUES'][0,:])
        self.kpoints      = np.array(database.variables['K-POINTS'][:].T)
        dimensions = database.variables['DIMENSIONS'][:]
        self.nbands      = int(dimensions[5])
        self.temperature = dimensions[13]
        self.num_electrons  = int(dimensions[14])
        self.nkpoints    = int(dimensions[6])
        self.spin = int(dimensions[11])
        self.time_rev = dimensions[9]
        database.close()

        #spin degeneracy if 2 components degen 1 else degen 2
        self.spin_degen = [0,2,1][int(self.spin)]

        #number of occupied bands
        self.nbands_valence = int(self.num_electrons/self.spin_degen)
        self.nbands_conduction = int(self.nbands-self.nbands_valence)

    def get_evals(self, set_gap = None, set_direct_gap = None):
        """
        Return the ks energies for each kpoint (in eV). The top of the valence band is used as the
        reference energy value. It is possible to shift the energies of the empty bands by setting an arbitrary
        value for the gap (direct or indirect).
        Implemented only for semiconductors, the energy shift of empty bands does not update their occupation levels.

        Args:
            set_gap (float) : set the value of the gap (in eV) of the system
            set_direct_gap (float) : set the value of the direct gap (in eV) of the system. If set_gap
                            is provided this parameter is ignored

        Return:
            :py:class:`numpy.array`  : an array with the ks energies for each kpoint

        """
        evals = F.get_evals(self.evals,self.nbands,self.nbands_valence,set_gap=set_gap,set_direct_gap=set_direct_gap)
        return evals

    def get_transitions(self, initial = 'full', final = 'empty',set_gap = None, set_direct_gap = None):
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
            set_gap (float) : set the value of the gap (in eV) of the system
            set_direct_gap (float) : set the value of the direct gap (in eV) of the system. If set_gap
                            is provided this parameter is ignored

        Return:
            :py:class:`numpy.array`  : an array with the transition energies for each kpoint

        """
        transitions = F.get_transitions(self.evals,self.nbands,self.nbands_valence,initial=initial,final=final,set_gap=set_gap,set_direct_gap=set_direct_gap)
        return transitions

    def get_gap(self, verbose = True):
        """
        Compute the energy gap of the system (in eV). The method check if the gap is direct or
        indirect. Implemented and tested only for semiconductors.

        Return:
            :py:class:`dict` : a dictionary with the values of direct and indirect gaps and the positions
            of the VMB and CBM

        """
        gap = F.get_gap(self.evals,self.nbands_valence,verbose=verbose)
        return gap
