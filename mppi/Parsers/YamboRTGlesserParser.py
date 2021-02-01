"""
Module that manages the parsing of the ``ndb.RT_G_PAR`` database created by `yambo_rt`.
"""
from netCDF4 import Dataset
import numpy as np

def buildBlochVectors(dm):
    """
    Build the Bloch vectors :math:`u_i(t,k)` (for each time step and for each k point) associated to the
    2x2 density matrix provided as input. Bloch vectors are computed according to the relation

    .. math::
        u_1(t,k) = 2Re(\\rho_{0,1,k}(t)) \, , \\
        u_2(t,k) = -2Im(\\rho_{0,1,k}(t)) \, , \\
        u_3(t,k) = \\rho_{1,1,k}(t) - \\rho_{0,0,k}(t)

    Note that the definition of :math:`u_3` respects the condition `excited state occupation` minus
    `ground state occupation` (the index is inverted with respect to the TLS where 1 is excited state
    and 2 is the ground state). The Bloch vectors are rescaled by the value of the trace of the density
    matrix (for each time and k point). If one trace vanishes the vector associated to the same time
    and k is set to zero

    Args:
        dm (:py:class:`array`) : time and k dependent array(s) with the 2x2 density matrix. The
            function assumes that the shape of dm is (time,kpoint,2,2)

    Returns:
        (:py:class:`array`) : the `real` array with the Bloch vectors for all the times and k points.
            The shape of the array is (component,time,kpoint)

    """
    trace = np.trace(dm,axis1=2,axis2=3).real
    ntime,nk = dm.shape[0],dm.shape[1]
    u1 = 2*dm[:,:,0,1].real
    u2 = -2*dm[:,:,0,1].imag
    u3 = dm[:,:,1,1].real-dm[:,:,0,0].real
    Bloch = np.array([u1,u2,u3])
    for t in range(ntime):
        for k in range(nk):
            if trace[t,k] == 0:
                Bloch[:,t,k] = 0.
            else:
                Bloch[:,t,k] = Bloch[:,t,k]/trace[t,k]
    return Bloch

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

    def buildEqDensityMatrix(self, num_occupied_bands = 1, band_occupation = 2):
        """
        Compute the equilibrium density matrix, that is a (time and k independent)
        diagonal array filled with the occupations of the bands

        Args:
            num_occupied_bands (:py:class:`int`) : number of occupied bands
            occupation (:py:class:`int`) : set the number of electrons in the
                occupied bands

        Returns:
            :py:class:`array` : diagonal array (of the dimension of the Green function)
                with the equilibrium occupations
        """
        dimension = self.Gless.shape[2] # number of bands
        occupations = [band_occupation if band < num_occupied_bands else 0 for band in range(dimension)]
        dm0 = np.diag(occupations)
        return dm0

    def buildDensityMatrix(self, first_time = None, last_time = None, first_k = None, last_k = None,
        first_band = None, last_band = None, equilibrium_dm = None):
        """
        Compute the density matrix from the Green function according to the
        relation

        .. math::
            \\rho_{b1,b2,k}(t) = -iG^<_{b1,b2,k}(t)

        The real and complex parts of the ``dG`` array are recasted to produce the
        complex structure of the density matrix.
        If the ``equilibrium_dm`` argument is not None add the equilibrium occupations
        provided in the argument to the density matrix.

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
            equilibrium_dm (:py:class:`array`) : if not None include the equilibrium
                occupations. The equilibrium_dm has to be an array with shape equal
                to (band1,band2) and the slice corresponding to the block of the Green
                function is extracted

        Returns:
            (:py:class:`array`) : the `complex` array with the density matrix.
                The structure of the array is [time,kpoint,band1,band2]

        """
        dm = - 1j*self.Gless[slice(first_time,last_time),slice(first_k,last_k),slice(first_band,last_band),slice(first_band,last_band),0] + \
            self.Gless[slice(first_time,last_time),slice(first_k,last_k),slice(first_band,last_band),slice(first_band,last_band),1]
        if equilibrium_dm is not None:
            dm[:,:] = dm[:,:]+equilibrium_dm[slice(first_band,last_band),slice(first_band,last_band)]
        return dm
