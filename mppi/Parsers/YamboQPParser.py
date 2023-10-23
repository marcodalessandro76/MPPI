"""
Module that manages the parsing of the ``ndb.RT_carriers`` database created by `yambo_rt`.
"""

from netCDF4 import Dataset

import numpy as np

class YamboQPParser():
    """
    Class that manages the quasi-particle correction written in the ``ndb.QP`` database
    created by `yambo`.

    Args:
        file (:py:class:`string`) : string with the name of the database to be parsed
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        QP_table (:py:class:`np.array`) : Array ...
        QP_kpts (:py:class:`np.array`) : Array ...
        QP_E (:py:class:`np.array`) : Array with the quasi-particle energies (in Hartree)
        QP_Eo (:py:class:`np.array`) : Array with the dft energies (in Hartree)
        QP_Z (:py:class:`np.array`) : Array ...
    """

    def __init__(self,file,verbose=True):
        self.filename = file
        if verbose: print('Parse file : %s'%self.filename)
        self.readDB(verbose)

    def readDB(self,verbose):
        """
        Read the data from the ``ndb.RT_carriers`` database created by `yambo_rt`. The variables
        are extracted from the database and stored in the attributes of the object.

        Args:
            verbose (:py:class:`boolean`) : define the amount of information provided on terminal

        """
        try:
            database = Dataset(self.filename)
        except:
            raise IOError("Error opening file %s in YamboQPParser"%self.filename)
        try:
             qp_test = database['QP_E']
        except IndexError:
             print('Old version of database detected')
             raise IndexError('Problem with the ndb.QP database')
        self.QP_table = np.array(database.variables['QP_table'][:].T)
        self.QP_kpts = np.array(database.variables['QP_kpts'][:].T)
        self.QP_E = np.array(database.variables['QP_E'][:])
        self.QP_Eo = np.array(database.variables['QP_Eo'][:])
        self.QP_Z = np.array(database.variables['QP_Z'][:])
        #print(database['PARS'])
        print(list(map(int,database['PARS'][:])))

    def get_info(self):
        """
        Provide information on the attributes of the class
        """
        print('YamboQPParser variables structure')
        print('QP_table shape',self.QP_table.shape)
        print('QP_kpts shape',self.QP_kpts.shape)
        print('QP_E shape',self.QP_E.shape)
        print('QP_Eo shape',self.QP_Eo.shape)
        print('QP_Z shape',self.QP_Z.shape)
