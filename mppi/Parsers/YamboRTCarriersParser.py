"""
Module that manages the parsing of the ``ndb.RT_carriers`` database created by `yambo_rt`.
"""
from netCDF4 import Dataset
from mppi.Utilities.Constants import HaToeV
from mppi.Utilities import Dos
from mppi.Utilities import LatticeUtils as latUtils

import numpy as np

def distance_point_to_segment(point,start,end):
    """
    Compute the distance between a point and a segment that
    has the points `start` and `end` as edges.
    The function checks if the orthogonal projection of the point to
    the segment line belongs to the segment itself. If not, the distance
    is computed as the minum value between the distances of the `point`
    w.r.t. the start and end points of the segment.
    Note that the function uses vectors analysis so it can be used in both
    2D and 3D.

    Args:
        point (:py:class:`np.array`) : cartesian coordinates of the point
        start (:py:class:`np.array`) : cartesian coordinates of the start point of the segment
        end (:py:class:`np.array`) : cartesian coordinates of the end point of the segment

    Returns :
        :py:class:`tuple` : the tuple (distance,ascissa_path) with the value of the distance
        and the value of the curvilinear ascissa on the line associated to the line element closest
        to the `point`

    """
    line_vec = end-start
    point_vec = point-start # coordinates of vec w.r.t. the start of the line
    line_unit = line_vec/np.linalg.norm(line_vec)
    dot_prod = np.dot(point_vec,line_unit)
    proj_vec = dot_prod*line_unit # projection of the point_vec in the direction of the segment
    if dot_prod < 0:
        distance = np.linalg.norm(point-start)
        ascissa_path = 0.
        return distance,ascissa_path
    if np.linalg.norm(proj_vec) > np.linalg.norm(line_vec):
        distance = np.linalg.norm(point-end)
        ascissa_path = np.linalg.norm(line_vec)
        return distance,ascissa_path
    else:
        distance = np.linalg.norm(proj_vec-point_vec)
        ascissa_path = np.linalg.norm(proj_vec)
        return distance,ascissa_path

class YamboRTCarriersParser():
    """
    Class to manage information about the real time distrubtion of carriers from the
    ``ndb.RT_carriers`` database created by `yambo_rt`.

    Args:
        file (:py:class:`string`) : string with the name of the database to be parsed
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        E_bare (:py:class:`np.array`) : Array that contains the bare bands energies.
            The structure of the array is [numk*numbands]. Energies are expressed in eV
            The structure is k1*b1,k1*b2,..,k1*bn,k2*b1,...
        f_bare (:py:class:`np.array`) : Array that contains the bare bands occupations.
            The structure of the array is [numk*numbands]
        kpoints (:py:class:`np.array`) : Array that contains the kpoints in cartesian  coordinates
            in units of 2*np.pi/alat (with a vector alat). The structure of the array is [numk,3],
            the first index runs over the kpoints and the second one gives the component
        bands_kpts (:py:class:`np.array`) : Array with the [first_band,last_band,numk]
        k_weight (:py:class:`np.array`) : Array that contains the weights of the kpoints.
            The structure of the array is [numk]
        delta_E (:py:class:`np.array`) : Array that contains the time-dependent correction
            to the bands energies. The structure of the array is [time,numk*numbands].
            Energies are expressed in eV
        delta_f (:py:class:`np.array`) : Array that contains the time-dependent corrections to
            the bands occupations. The structure of the array is [time,numk*numbands]

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
            raise IOError("Error opening file %s in YamboRTCarriersParser"%self.filename)
        self.E_bare = HaToeV*np.array(database.variables['RT_carriers_E_bare'])
        self.f_bare = np.array(database.variables['RT_carriers_f_bare'])
        self.kpoints = np.array(database.variables['RT_kpt'][:].T)
        self.bands_kpts = np.array(database.variables['RT_bands_kpts'])
        self.k_weight = np.array(database.variables['RT_k_weight'])
        self.delta_E = HaToeV*np.array(database.variables['RT_carriers_delta_E'])
        self.delta_f = np.array(database.variables['RT_carriers_delta_f'])

    def get_info(self):
        """
        Provide information on the attributes of the class
        """
        print('YamboRTCarriersParser variables structure')
        print('Bands used and number of k-points',self.bands_kpts)
        print('kpoints shape',self.kpoints.shape)
        print('E_bare shape',self.E_bare.shape)
        print('f_bare shape',self.f_bare.shape)
        print('delta_E shape',self.delta_E.shape)
        print('delta_f',self.delta_f.shape)

    def build_f_bare_dos(self, dE = 0.1, eta = 0.05, broad_kind = 'lorentzian'):
        """
        For each kpoint build a dos which expresses the bare occupation level in terms
        of the energy. The energy ranges from the minum to the maximum of the
        E_bare variable.

        Args:
            dE (:py:class:`float`) : energy step in eV
            eta (:py:class:`float`) : magnitude of the broading parameter (in the same units used for the values array)
            broad_kind (:py:class:`string`) : type of broading function used (lorentzian, gaussian)

        Returns:
            :py:class:`Dos` : Instance of the ``Dos`` class. The object is an array of dos, one for
                each kpoint

        """
        numbnds = self.bands_kpts[1]-self.bands_kpts[0]
        numkp = self.bands_kpts[2]
        Emin,Emax = min(self.E_bare)-10*eta,max(self.E_bare)+10*eta

        dos = Dos()
        for kind in range(numkp):
            E_bare_k = self.E_bare[kind*numbnds:(kind+1)*numbnds]
            f_bare_k = self.f_bare[kind*numbnds:(kind+1)*numbnds]
            dos.append(E_bare_k,weights=f_bare_k,minVal=Emin,maxVal=Emax,
                step=dE,eta=eta,broad_kind=broad_kind)

        return dos
