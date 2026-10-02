"""
Module that manages the parsing of the ``ndb.Nonlinear`` database created by `yambo_nl`.
The module is an (almost) identical clone of the ``nldb.py`` of yambopy.
"""

from netCDF4 import Dataset
from mppi.Utilities.Constants import HaToeV, Light_speed_au, FsToAu
from mppi.Utilities.Utils import Plot_3dArray
import numpy as np

def read_string(database, name):
    """
    Read a string variable of a netCDF database.
    """
    return database.variables[name][...].tobytes().decode().strip()

def read_value(database, name, default=None, dtype=np.double):
    """
    Read the first element of a variable of a netCDF database. If the variable is not present
    and a default is given, return the default.
    """
    try:
        return database.variables[name][0].astype(dtype)
    except (KeyError, IndexError):
        if default is None: raise
        return default

def read_array(database, name, dtype=np.double):
    """
    Read a variable of a netCDF database as a (not masked) numpy array.
    """
    return np.array(database.variables[name][:],dtype=dtype)

class YamboNLDBParser(object):
    """
    Read the ``ndb.Nonlinear`` database written by yambo_nl and its fragments ``ndb.Nonlinear_fragment_n``
    (one for each run, i.e. for each frequency or angle of the external field). The time series is stored in the
    IO_TIME_points attribute, while the external fields, the current and the polarization of the different runs
    are stored in the Efield, Polarization and Current lists.

    Args:
        file (:py:class:`string`) : name of the ``ndb.Nonlinear`` database (including the path)
        nl_db (:py:class:`string`) : not used, kept for backward compatibility
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        Gauge (:py:class:`string`) : gauge used in the computation
        NE_steps (:py:class:`int`) : number of time steps of the real-time propagation
        RT_step (:py:class:`float`) : time step of the real-time propagation (in atomic units)
        n_frequencies (:py:class:`int`) : number of frequencies of the external field
        n_angles (:py:class:`int`) : number of angles of the external field
        n_runs (:py:class:`int`) : number of runs (frequencies or angles) of the database
        IO_TIME_points (:py:class:`np.array`) : time points of the real-time propagation (in atomic units)
        Efield_general (:py:class:`list`) : dictionaries with the parameters of the external fields (1, 2, 3) found
            in the ``ndb.Nonlinear`` database
        N_ext_fields (:py:class:`int`) : number of external fields found in the database
        Polarization (:py:class:`list`) : for each run, array with shape (3,time) with the polarization
        Current (:py:class:`list`) : for each run, array with shape (3,time) with the current
        E_ext (:py:class:`list`) : for each run, complex array with shape (3,time) with the external field
        E_tot (:py:class:`list`) : for each run, complex array with shape (3,time) with the total field
        E_ks (:py:class:`list`) : for each run, complex array with shape (3,time) with the Kohn-Sham field
        Efield (:py:class:`list`) : for each run, dictionary with the parameters of the first external field
        Efield2 (:py:class:`list`) : for each run, dictionary with the parameters of the second external field (the
            first one if a second field is not present)

    """

    def __init__(self,file,nl_db='ndb.Nonlinear',verbose=True):
        self.nl_path = file
        if verbose: print('Parse file : %s'%self.nl_path)
        try:
            data = Dataset(self.nl_path)
        except OSError:
            raise ValueError("Error reading NONLINEAR database at %s"%self.nl_path)
        self.read_observables(data)
        data.close()

    def read_Efield(self,database,RT_step,n):
        """
        Read the parameters of the n-th external field. Raise a KeyError if the field is not present.

        Args:
            database (:py:class:`netCDF4.Dataset`) : the database
            RT_step (:py:class:`float`) : time step of the real-time propagation (in atomic units)
            n (:py:class:`int`) : index of the field

        Returns:
            :py:class:`dict` : the parameters of the field (in atomic units)

        """
        n = str(n)
        efield = {}
        efield["name"]       = read_string(database,'Field_Name_'+n)
        efield["versor"]     = read_array(database,'Field_Versor_'+n)
        efield["intensity"]  = read_value(database,'Field_Intensity_'+n)
        efield["damping"]    = read_value(database,'Field_FWHM_'+n,
                                          default=read_value(database,'Field_Damping_'+n,default=0.))
        try:
            efield["freq_range"] = read_array(database,'Field_Freq_range_'+n)
        except KeyError:
            efield["freq_range"] = read_array(database,'Field_Freq_'+n)
        try:
            efield["freq_steps"] = read_array(database,'Field_Freq_steps_'+n)
        except KeyError:
            efield["freq_steps"] = 1.0
        efield["freq_step"]    = read_value(database,'Field_Freq_step_'+n,default=0.)
        efield["initial_time"] = read_value(database,'Field_Initial_time_'+n)
        efield["peak"]         = read_value(database,'Field_peak_'+n,default=10.)

        # set t_initial according to Yambo
        efield["initial_indx"] = max(round(efield["initial_time"]/RT_step)+1,2)
        efield["initial_time"] = (efield["initial_indx"]-1)*RT_step

        # define the field amplitude
        efield["amplitude"] = np.sqrt(efield["intensity"]*4.0*np.pi/Light_speed_au)

        return efield

    def read_observables(self,database):
        """
        Read all data from the database and from its fragments
        """
        self.Gauge          = read_string(database,'GAUGE')
        self.NE_steps       = read_value(database,'NE_steps',dtype='int')
        self.RT_step        = read_value(database,'RT_step')
        self.n_frequencies  = read_value(database,'n_frequencies',dtype='int')
        self.n_angles       = read_value(database,'n_angles',default=0,dtype='int')
        try:
            self.NL_initial_versor = read_array(database,'NL_initial_versor')
        except KeyError:
            self.NL_initial_versor = np.zeros(3)
        self.NL_damping     = read_value(database,'NL_damping')
        self.RT_bands       = read_array(database,'RT_bands',dtype='int')
        self.NL_er          = read_array(database,'NL_er')
        self.l_force_SndOrd = read_value(database,'l_force_SndOrd',dtype='bool')
        self.l_use_DIPOLES  = read_value(database,'l_use_DIPOLES',dtype='bool')
        self.l_eval_CURRENT = read_value(database,'l_eval_CURRENT',default=False,dtype='bool')
        self.QP_ng_SH       = read_value(database,'QP_ng_SH',dtype='int')
        self.QP_ng_Sx       = read_value(database,'QP_ng_Sx',dtype='int')
        self.RAD_LifeTime   = read_value(database,'RAD_LifeTime')
        self.Integrator     = read_string(database,'Integrator')
        self.Correlation    = read_string(database,'Correlation')

        # Time variables
        self.IO_TIME_N_points  = read_value(database,'IO_TIME_N_points',dtype='int')
        self.IO_TIME_LAST_POINT= read_value(database,'IO_TIME_LAST_POINT',dtype='int')
        self.IO_TIME_points    = read_array(database,'IO_TIME_points')

        # External fields (the fields not present in the database are skipped)
        self.Efield_general = []
        for n in range(1,4):
            try:
                self.Efield_general.append(self.read_Efield(database,self.RT_step,n))
            except KeyError:
                pass
        self.N_ext_fields = len(self.Efield_general)

        # Number of runs
        if self.n_angles != 0 and self.n_frequencies != 0:
            raise ValueError('Both n_angles and n_frequencies are different from zero in %s'%self.nl_path)
        self.n_runs = max(self.n_angles,self.n_frequencies,1)

        # Polarization, current and fields of each run
        self.Polarization = []
        self.Current      = []
        self.E_ext        = []
        self.E_tot        = []
        self.E_ks         = []
        self.Efield       = [] # the first external field of each run
        self.Efield2      = [] # the second external field of each run
        for f in range(1,self.n_runs+1):
            fragment = self.nl_path+"_fragment_"+str(f)
            try:
                data_p_and_j = Dataset(fragment)
            except OSError:
                print("Error reading database: %s"%fragment)
                continue
            freq = str(f).zfill(4)
            self.Polarization.append(read_array(data_p_and_j,'NL_P_freq_'+freq))
            self.Current.append(read_array(data_p_and_j,'NL_J_freq_'+freq))
            for name,store in (('E_ext',self.E_ext),('E_tot',self.E_tot),('E_ks',self.E_ks)):
                field = read_array(data_p_and_j,name+'_freq_'+freq)
                store.append(field[:,:,0]+1j*field[:,:,1])
            efield = self.read_Efield(data_p_and_j,self.RT_step,1)
            try:
                efield2 = self.read_Efield(data_p_and_j,self.RT_step,2)
            except KeyError:
                efield2 = efield
            self.Efield.append(efield)
            self.Efield2.append(efield2)
            data_p_and_j.close()

    def get_info(self):
        """
        Provide information on the attributes of the class
        """
        print('YamboNLDBParser variables structure')
        s="\n * * * ndb.Nonlinear db data * * * \n\n"
        s+="Gauge             : "+str(self.Gauge)+"\n"
        s+="NE_steps          : "+str(self.NE_steps)+"\n"
        s+="RT_step           : "+str(self.RT_step/FsToAu)+" [fs] \n"
        s+="n_frequencies     : "+str(self.n_frequencies)+"\n"
        s+="n_angles          : "+str(self.n_angles)+"\n"
        s+="NL_initial_versor : "+str(self.NL_initial_versor)+"\n"
        s+="NL_damping        : "+str(self.NL_damping*HaToeV)+" [eV] \n"
        s+="RT_bands          : "+str(self.RT_bands)+"\n"
        s+="NL_er             : "+str(self.NL_er*HaToeV)+" [eV] \n"
        s+="Second Order      : "+str(self.l_force_SndOrd)+"\n"
        s+="Use Dipoles       : "+str(self.l_use_DIPOLES)+"\n"
        s+="QP_ng_SH          : "+str(self.QP_ng_SH)+"\n"
        s+="QP_ng_Sx          : "+str(self.QP_ng_Sx)+"\n"
        s+="RAD_LifeTime      : "+str(self.RAD_LifeTime/FsToAu)+" [fs] \n"
        s+="Integrator        : "+str(self.Integrator)+"\n"
        s+="Correlation       : "+str(self.Correlation)+"\n"
        s+="External fields   : "+str(self.N_ext_fields)+"  runs : "+str(self.n_runs)+"\n"
        for efield in self.Efield:
            if efield["name"] == "none":
                continue
            s+="\nEfield name         : "+str(efield["name"])+"\n"
            s+="Efield versor       : "+str(efield["versor"])+"\n"
            s+="Efield Intensity    : "+str(efield["intensity"])+"\n"
            s+="Efield Damping      : "+str(efield["damping"])+"\n"
            s+="Efield Freq range   : "+str(efield["freq_range"]*HaToeV)+" [eV] \n"
            s+="Efield Initial time : "+str(efield["initial_time"]/FsToAu)+" [fs] \n"
        print(s)

    def get_time(self,convert_to_fs=True):
        """
        Get the time points of the real-time propagation.

        Args:
            convert_to_fs (:py:class:`boolean`) : if True, convert the time points from atomic units to femto-seconds

        Returns:
            :py:class:`array` : array with the time points of the real-time propagation
        """
        if convert_to_fs:
            return self.IO_TIME_points/FsToAu
        return self.IO_TIME_points

    def plot_polarization(self,convert_to_fs=True,xlim=None,run_index=0):
        """
        Plot the polarization of the system as a function of time.

        Args:
            convert_to_fs (:py:class:`boolean`) : if True, convert the time points from atomic units to femto-seconds
            xlim (:py:class:`tuple`) : tuple with the limits of the x-axis
            run_index (:py:class:`int`) : index of the run for which to plot the polarization
        """
        Plot_3dArray(self.get_time(convert_to_fs),self.Polarization[run_index],xlim=xlim,label='Polarization')

    def plot_current(self,convert_to_fs=True,xlim=None,run_index=0):
        """
        Plot the current of the system as a function of time.

        Args:
            convert_to_fs (:py:class:`boolean`) : if True, convert the time points from atomic units to femto-seconds
            xlim (:py:class:`tuple`) : tuple with the limits of the x-axis
            run_index (:py:class:`int`) : index of the run for which to plot the current
        """
        Plot_3dArray(self.get_time(convert_to_fs),self.Current[run_index],xlim=xlim,label='J')
