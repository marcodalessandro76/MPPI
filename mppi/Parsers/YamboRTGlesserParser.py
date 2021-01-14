"""
Module that manages the parsing of the ``ndb.RT_G_PAR`` database created by `yambo_rt`.
"""
from netCDF4 import Dataset
import numpy as np

class YamboRTGlesserParser():
    """
    Class to manage information about the :math:`G^<` from the ``ndb.RT_G_PAR`` database
    created by `yambo_rt`.

    Args:
        file (:py:class:`string`) : string with the name of the database to be parsed
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        Gless (:py:class:`np.array`): Array with the variable ``dG_lesser`` that contains the
            lesser component of the time-dependent Green function expressed in the Kohn-Sham basis.
            The structure of the array is [time,kpoint,band1,band2,real and imaginary part]

    """

    def __init__(self,file,verbose=True):
        self.filename = file
        if verbose: print('Parse file : %s'%self.filename)
        self.readDB(verbose)

    def readDB(self,verbose):
        """
        Read the data from the ``ndb.RT_G_PAR`` database created by `yambo_rt`. The variables
        are extracted from the database and stored in the attributes of the object.

        Args:
            verbose (:py:class:`boolean`) : define the amount of information provided on terminal

        """
        try:
            database = Dataset(self.filename)
        except:
            raise IOError("Error opening file %s in YamboRTGlesserParser"%self.filename)
        self.Gless = np.array(database.variables['dG_lesser'])

    def get_info(self):
        """
        Provide information on the attributes of the class
        """
        print('YamboRTGlesserParser variables structure')
        print('Gless shape',self.Gless.shape)

    def getDensityMatrix(self, first_time = None, last_time = None, first_k = None, last_k = None,
        first_band = None, last_band = None):
        """
        Compute the density matrix from the Green function according to the
        relation

        .. math::
            \\rho_{b1,b2,k}(t) = -iG^<_{b1,b2,k}(t)

        The real and complex parts of the ``dG`` array are recasted to produce the
        complex structure of the density matrix.

        Args:
            first_time (:py:class:`int`) : index of the first time variable. If None the
                first value of the Green function is used
            last_time (:py:class:`int`) : index of the last time variable. If None the
                last value of the Green function is used
            first_k (:py:class:`int`) : index of the first k variable. If None the
                first value of the Green function is used
            last_k (:py:class:`int`) : index of the last k variable. If None the
                last value of the Green function is used
            first_band (:py:class:`int`) : index of the first band. If None the
                first value of the Green function is used
            last_band (:py:class:`int`) : index of the last band. If None the
                last value of the Green function is used

        Returns:
            (:py:class:`array`) : the `complex` array with the density matrix.
                The structure of the array is [time,kpoint,band1,band2]

        """
        # define the density matrix in the conduction sector as -i*G_<
        dm = - 1j*self.Gless[slice(first_time,last_time),slice(first_k,last_k),slice(first_band,last_band),slice(first_band,last_band),0] + \
            self.Gless[slice(first_time,last_time),slice(first_k,last_k),slice(first_band,last_band),slice(first_band,last_band),1]
        return dm

    def evalDensityMatrixTrace(self, first_time = None, last_time = None, first_k = None, last_k = None,
        first_band = None, last_band = None):
        """
        Compute the trace (for each time step and for each k point) of the density matrix. Useful for the computation of
        the expectation values and for the construction of the Bloch vectors.

        Args:
            first_time (:py:class:`int`) : index of the first time variable. If None the
                first value of the Green function is used
            last_time (:py:class:`int`) : index of the last time variable. If None the
                last value of the Green function is used
            first_k (:py:class:`int`) : index of the first k variable. If None the
                first value of the Green function is used
            last_k (:py:class:`int`) : index of the last k variable. If None the
                last value of the Green function is used
            first_band (:py:class:`int`) : index of the first band. If None the
                first value of the Green function is used
            last_band (:py:class:`int`) : index of the last band. If None the
                last value of the Green function is used

        Returns:
            (:py:class:`array`) : the `real` array with the traces of the (last_band-first_band)x(last_band-first_band)
                bloch of the density matrix. The structure of the array is [time,kpoint]

        """
        dm = self.getDensityMatrix(first_time,last_time,first_k,last_k,first_band,last_band)
        ntime,nk = dm.shape[0],dm.shape[1]
        trace = np.zeros([ntime,nk])
        for t in range(ntime):
                for k in range(nk):
                    trace[t][k] = np.trace(dm[t][k]).real
        return trace

    def buildBlochVectors(self, first_time = None, last_time = None, first_k = None, last_k = None,
        first_band = 0):
        """
        Build the Bloch vectors :math:`u_i(t,k)` (for each time step and for each k point) associated to a 2x2 bloch of the density
        matrix specified by the ``first_band`` parameter. Bloch vectors are computed according to the relation

        .. math::
            u_1(t,k) = 2Re(\\rho_{b1,b1+1,k}(t)) \, , \\
            u_2(t,k) = -2Im(\\rho_{b1,b1+1,k}(t)) \, , \\
            u_3(t,k) = \\rho_{b1,b1,k}(t) - \\rho_{b1+1,b1+1,k}(t)

        The Bloch vectors are rescaled by the value of the trace of the density matrix (for each time and k point). If one trace
        vanishes the vector associated to the same time and k is set to zero

        Args:
            first_time (:py:class:`int`) : index of the first time variable. If None the
                first value of the Green function is used
            last_time (:py:class:`int`) : index of the last time variable. If None the
                last value of the Green function is used
            first_k (:py:class:`int`) : index of the first k variable. If None the
                first value of the Green function is used
            last_k (:py:class:`int`) : index of the last k variable. If None the
                last value of the Green function is used
            first_band (:py:class:`int`) : index of the first band. The 2x2 bloch of the
                dm starts at this index value

        Returns:
            (:py:class:`array`) : the `real` array with the Bloch vectors for all the times and k points.
                The structure of the array is [component,time,kpoint]

        """
        dm = self.getDensityMatrix(first_time,last_time,first_k,last_k,first_band,first_band+2)
        trace = self.evalDensityMatrixTrace(first_time,last_time,first_k,last_k,first_band,first_band+2)
        ntime,nk = dm.shape[0],dm.shape[1]
        u1 = 2*dm[:,:,0,1].real
        u2 = -2*dm[:,:,0,1].imag
        u3 = dm[:,:,0,0].real-dm[:,:,1,1].real
        Bloch = np.array([u1,u2,u3])
        for t in range(ntime):
            for k in range(nk):
                tr = trace[t,k]
                if tr == 0:
                    Bloch[:,t,k] = 0.
                else:
                    Bloch[:,t,k] = Bloch[:,t,k]/tr
        return Bloch
