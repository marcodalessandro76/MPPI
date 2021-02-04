"""
Module that manages the parsing of a Yambo ``dipoles`` database.
"""
from netCDF4 import Dataset
import numpy as np

class YamboDipolesParser():
    """
    Class to read information about the dipoles from the ``dipoles`` database
    created by Yambo.

    Args:
        file (:py:class:`string`) : string with the name of the database to be parsed
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        dip_ir (:py:class:`np.array`): Array with the variable ``DIP_iR`` that contains the position
            dipoles times the imaginary unit. The structure of the array is
            [kpoints][full_bands][empty_bands][cartesian_component][real and imaginary part]
        dip_p (:py:class:`np.array`): Array with the variable ``DIP_p`` that contains the matrix
            elements of the momentum operator. The structure of the array
            is [kpoints][full_bands][empty_bands][cartesian_component][real and imaginary part]
        dip_v (:py:class:`np.array`): Array with the variable ``DIP_v`` that contains the matrix
            elements of the velocity operator (commutator of x operator with the Hamiltonian).
            The structure of the array is
            [kpoints][full_bands][empty_bands][cartesian_component][real and imaginary part]
        dip_spin (:py:class:`np.array`): Array with the variable ``DIP_spin`` that contains the
            matrix elements of the spin operator. The structure of the array
            is [kpoints][full_bands][empty_bands][cartesian_component][real and imaginary part]

    Note:
        If the option ``DipBandsALL`` is enabled the dipoles are computed for all the bands range,
        not only for the valence-conduction pairs.

    """

    def __init__(self,file,verbose=True):
        self.filename = file
        if verbose: print('Parse file : %s'%self.filename)
        self.readDB(verbose)

    def readDB(self,verbose):
        """
        Read the data from the ``dipoles`` database created by Yambo. The variables
        are extracted from the database and stored in the attributes of the object.

        Args:
            verbose (:py:class:`boolean`) : define the amount of information provided on terminal

        """
        try:
            database = Dataset(self.filename)
        except:
            raise IOError("Error opening file %s in YamboDipolesParser"%self.filename)

        self.dip_ir = np.array(database.variables['DIP_iR'][0])
        self.dip_p = np.array(database.variables['DIP_P'][0])
        self.dip_v = np.array(database.variables['DIP_v'][0])
        try:
            self.dip_spin = np.array(database.variables['DIP_spin'][0])
        except KeyError:
            if verbose: print('Spin dipoles not found in the ndb.dipoles')
            self.dip_spin = np.array([0])

    def get_info(self):
        """
        Provide information on the attributes of the class
        """
        print('YamboDipolesParser variables structure')
        print('dip_ir shape',self.dip_ir.shape)
        print('dip_p shape',self.dip_v.shape)
        print('dip_v shape',self.dip_v.shape)
        print('dip_spin shape',self.dip_spin.shape)

    def r_dipole(self,component,first_k = None, last_k = None,
        first_band = None, last_band = None):
        """
        Compute the array with the matrix elements of the position operator. The `i` factor in the
        definition of the dip_ir attribute produces a flip of the real and imaginary
        components (with a change of sign).

        Args:
            component (py:class:`int`)   : cartesian component of the dipole operator
            first_k (:py:class:`int`) : index of the first k variable. If None the
                first value of the dip_ir variable is used
            last_k (:py:class:`int`) : index of the last k variable. If None the
                last value of the dip_ir variable is used
            first_band (:py:class:`int`) : index of the first band. If None the
                first value of the dip_ir variable is used
            last_band (:py:class:`int`) : index of the last band. If None the
                last value of the dip_ir variable is used

        Returns:
            :py:class:`array` : complex array with the dipoles matrix element
            with shape(kpoint,band1,band2)

        """
        r_real = self.dip_ir[slice(first_k,last_k),slice(first_band,last_band),slice(first_band,last_band),component,1]
        r_imag = -self.dip_ir[slice(first_k,last_k),slice(first_band,last_band),slice(first_band,last_band),component,0]
        return r_real+1j*r_imag

    def spin_dipole(self,component,first_k = None, last_k = None,
        first_band = None, last_band = None):
        """
        Compute the array with the matrix elements of the spin operator.

        Args:
            component (py:class:`int`)   : cartesian component of the spin operator
            first_k (:py:class:`int`) : index of the first k variable. If None the
                first value of the dip_spin variable is used
            last_k (:py:class:`int`) : index of the last k variable. If None the
                last value of the dip_spin variable is used
            first_band (:py:class:`int`) : index of the first band. If None the
                first value of the dip_spin variable is used
            last_band (:py:class:`int`) : index of the last band. If None the
                last value of the dip_spin variable is used

        Returns:
            :py:class:`array` : complex array with the spin dipole matrix element
            with shape(kpoint,band1,band2)

        """
        spin_real = self.dip_spin[slice(first_k,last_k),slice(first_band,last_band),slice(first_band,last_band),component,0]
        spin_imag = self.dip_spin[slice(first_k,last_k),slice(first_band,last_band),slice(first_band,last_band),component,1]
        return spin_real+1j*spin_imag
