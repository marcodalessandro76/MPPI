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
        dip_v (:py:class:`np.array`): Array with the variable ``DIP_v`` that contains the matrix
            elements of the velocity operator. The structure of the array
            is [kpoints][full_bands][empty_bands][cartesian_component][real and imaginary part]
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
        print('dip_v shape',self.dip_v.shape)
        print('dip_spin shape',self.dip_spin.shape)

    def get_r_dipole(self,kpoint,band_l,band_r,cart):
        """
        Get the matrix element <band_l|r[cart]|band_r> of the position operator.
        The `i` factor in the definition of the dip_ir attribute produces a flip
        of the components (with a change of sign).

        Args:
            kpoint (py:class:`int`) : kpoint
            band1_l (py:class:`int`)  : left band of the bracket
            band_r (py:class:`int`)  : right band of the bracket
            cart (py:class:`int`)   : cartesian component of the dipole operator

        Returns:
            :py:class:`array` : array with the real and imaginary part of the matrix element

        """
        ir = self.dip_ir[kpoint][band_l][band_r][cart]
        return np.array([ir[1],-ir[0]])
