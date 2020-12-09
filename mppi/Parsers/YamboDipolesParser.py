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
        verbose (:py:class:`boolean`) : Determine the amount of information provided on terminal

    Attributes:
        dip_r (:py:class:`np.array`): Array with the position dipoles. The structure of the array
            is [kpoints][full_bands][empty_bands][cartesian_component][real and imaginary part]
        dip_v (:py:class:`np.array`): Array with the velocity dipoles. The structure of the array
            is [kpoints][full_bands][empty_bands][cartesian_component][real and imaginary part]
        dip_spin (:py:class:`np.array`): Array with the spin dipoles. The structure of the array
            is [kpoints][full_bands][empty_bands][cartesian_component][real and imaginary part]

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
            verbose (:py:class:`boolean`) : Determine the amount of information provided on terminal

        """
        try:
            database = Dataset(self.filename)
        except:
            raise IOError("Error opening file %s in YamboDipolesParser"%self.filename)

        self.dip_r = np.array(database.variables['DIP_iR'][0])
        self.dip_v = np.array(database.variables['DIP_v'][0])
        try:
            self.dip_spin = np.array(database.variables['DIP_spin'][0])
        except KeyError:
            if verbose: print('Spin dipoles not found in the ndb.dipoles')
            self.dip_spin = np.array([0])

    def get_info(self):
        """
        Provide information on the atributes of the class
        """
        print('YamboDipolesParser variables structure')
        print('dip_r shape',self.dip_r.shape)
        print('dip_v shape',self.dip_v.shape)
        print('dip_spin shape',self.dip_spin.shape)
