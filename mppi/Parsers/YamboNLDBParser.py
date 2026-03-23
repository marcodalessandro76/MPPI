"""
Module that manages the parsing of the ``ndb.Nonlinear`` database created by `yambo_nl`.
The module is an (almost) identical clone of the ``nldb.py`` of yambopy.  
"""

from netCDF4 import Dataset
from mppi.Utilities.Constants import HaToeV, Light_speed_au, FsToAu
from mppi.Utilities.Utils import Plot_3dArray  
import numpy as np
import sys, os

class YamboNLDBParser(object):
    """
    Open the NL databases and store it in a NLDB class.
    This class reads all data from the ndb.Nonlinear database
    and its fragment ndb.Nonlinear_fragment_xxx
    Time series is stored in IO_TIME_points variable
    while external fields, current and polarization of the different runs 
    are stored in Efield[x],Polarization[x],Current[x].
    Args:
        file (:py:class:`string`) : string with the name of the database to be parsed
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        Gauge (:py:class:`string`) : string with the gauge used in the calculation
        NE_steps (:py:class:`int`) : number of time steps in the real-time propagation
        RT_step (:py:class:`float`) : time step in the real-time propagation (in atomic units)
        n_frequencies (:py:class:`int`) : number of frequencies in the frequency-dependent calculation
        n_angles (:py:class:`int`) : number of angles in the angle-dependent calculation
        polarization (:py:class:`np.array`) : array with the polarization of the system at each time step
        current (:py:class:`np.array`) : array with the current of the system at each time step
        E_ext (:py:class:`np.array`) : array with the external field of the system at each time step
        E_tot (:py:class:`np.array`) : array with the total field of the system at each time step
        E_ks (:py:class:`np.array`) : array with the Kohn-Sham field of the system at each time step
        Efield (:py:class:`list`) : list of dictionaries with the information on the external fields used in the calculation           
    """

    def __init__(self,file,nl_db='ndb.Nonlinear',verbose=True):
        self.nl_path = file
        if verbose: print('Parse file : %s'%self.nl_path)
        try:
            data= Dataset(self.nl_path)
        except:
            raise ValueError("Error reading NONLINEAR database at %s"%self.nl_path)

        self.read_observables(data)
        data.close()

    def read_Efield(self,database,RT_step,n):
         efield={}
         efield["name"]       =database.variables['Field_Name_'+str(n)][...].tobytes().decode().strip()
         efield["versor"]     =database.variables['Field_Versor_'+str(n)][:].astype(np.double)
         efield["intensity"]  =database.variables['Field_Intensity_'+str(n)][0].astype(np.double)
         try:
            efield["damping"]    =database.variables['Field_FWHM_'+str(n)][0].astype(np.double)
         except:      
            efield["damping"]    =database.variables['Field_Damping_'+str(n)][0].astype(np.double)
         try:
             efield["freq_range"] =database.variables['Field_Freq_range_'+str(n)][:].astype(np.double)
         except:
             efield["freq_range"] =database.variables['Field_Freq_'+str(n)][:].astype(np.double)
         try:
             efield["freq_steps"] =database.variables['Field_Freq_steps_'+str(n)][:].astype(np.double)
         except:
             efield["freq_steps"] =1.0
         try:
             efield["freq_step"]  =database.variables['Field_Freq_step_'+str(n)][0].astype(np.double)
         except:
             efield["freq_step"]  =0.0
         efield["initial_time"]  =database.variables['Field_Initial_time_'+str(n)][0].astype(np.double)
         try:
             efield["peak"]  =database.variables['Field_peak_'+str(n)][0].astype(np.double)
         except:
             efield["peak"]  =10.0
         #
         # set t_initial according to Yambo 
         #
         efield["initial_indx"] =max(round(efield["initial_time"]/RT_step)+1,2)
         efield["initial_time"] =(efield["initial_indx"]-1)*RT_step
         #
         # define the field amplitude
         #
         efield["amplitude"]    =np.sqrt(efield["intensity"]*4.0*np.pi/Light_speed_au)

         return efield

    def read_observables(self,database):
        """
        Read all data from the database
        """
        self.Gauge          = database.variables['GAUGE'][...].tobytes().decode().strip()
        self.NE_steps       = database.variables['NE_steps'][0].astype('int')
        self.RT_step        = database.variables['RT_step'][0].astype(np.double)
        self.n_frequencies  = database.variables['n_frequencies'][0].astype('int')

        try:
            self.n_angles       = database.variables['n_angles'][0].astype('int')
        except:
            self.n_angles   = 0
        try:
            self.NL_initial_versor = database.variables['NL_initial_versor'][:].astype(np.double)
        except:
            self.NL_initial_versor = [0.0, 0.0, 0.0] 

        self.NL_damping     = database.variables['NL_damping'][0].astype(np.double)
        self.RT_bands       = database.variables['RT_bands'][:].astype('int')
        self.NL_er          = database.variables['NL_er'][:].astype(np.double)
        self.l_force_SndOrd = database.variables['l_force_SndOrd'][0].astype('bool')
        self.l_use_DIPOLES  = database.variables['l_use_DIPOLES'][0].astype('bool')
        try:
            self.l_eval_CURRENT = database.variables['l_eval_CURRENT'][0].astype('bool')
        except:
            self.l_eval_CURRENT = False
        self.QP_ng_SH       = database.variables['QP_ng_SH'][0].astype('int')
        self.QP_ng_Sx       = database.variables['QP_ng_Sx'][0].astype('int')
        self.RAD_LifeTime   = database.variables['RAD_LifeTime'][0].astype(np.double)
        self.Integrator     = database.variables['Integrator'][...].tobytes().decode().strip()
        self.Correlation    = database.variables['Correlation'][...].tobytes().decode().strip()
        #
        # Time variables
        #
        self.IO_TIME_N_points  = database.variables['IO_TIME_N_points'][0].astype('int')
        self.IO_TIME_LAST_POINT= database.variables['IO_TIME_LAST_POINT'][0].astype('int')
        self.IO_TIME_points    = database.variables['IO_TIME_points'][:].astype(np.double)
        #
        # External fields
        # 
        self.Efield_general=[]
        self.N_ext_fields=0
        for n in range(1,4):
            try:
                efield=self.read_Efield(database,self.RT_step,n)
            except:
                print("Field %d not found" % n)
            else:
                self.Efield_general.append(efield.copy())
                self.N_ext_fields=self.N_ext_fields+1

        #
        # Read polarization and currect files 
        #
        self.Polarization=[]
        self.Current     =[]
        self.E_ext       =[]
        self.E_tot       =[]
        self.E_ks        =[]
        self.Efield      =[] # Store the first external field for each run at different frequencies
        self.Efield2     =[]
        #
        if self.n_angles!=0:
            self.n_runs=self.n_angles
        if self.n_frequencies!=0:
            self.n_runs=self.n_frequencies
        if (self.n_angles!=0 and self.n_frequencies!=0):
            print("Error both n_angles and n_frequencies !=0 ")
            sys.exit(0)
        if (self.n_angles==0 and self.n_frequencies==0):
            self.n_runs=1

        #
        for f in range(self.n_runs):
            try:
                data_p_and_j= Dataset(self.nl_path+"_fragment_"+str(f+1))
            except:
                print("Error reading database: %s" % self.nl_path+"_fragment_"+str(f+1))
                continue
            pol  = data_p_and_j.variables['NL_P_freq_'+str(f+1).zfill(4)][:,:].astype(np.double)
            curr = data_p_and_j.variables['NL_J_freq_'+str(f+1).zfill(4)][:,:].astype(np.double)
            e_ext= data_p_and_j.variables['E_ext_freq_'+str(f+1).zfill(4)][:,:,:].astype(np.double)
            e_ext_c=e_ext[:,:,0]+1j*e_ext[:,:,1]
            e_tot= data_p_and_j.variables['E_tot_freq_'+str(f+1).zfill(4)][:,:,:].astype(np.double)
            e_tot_c=e_tot[:,:,0]+1j*e_tot[:,:,1]
            e_ks = data_p_and_j.variables['E_ks_freq_'+str(f+1).zfill(4)][:,:,:].astype(np.double)
            e_ks_c=e_ks[:,:,0]+1j*e_ks[:,:,1]

            self.Polarization.append(pol.copy())
            self.Current.append(curr.copy())
            self.E_ext.append(e_ext_c.copy())
            self.E_tot.append(e_tot_c.copy())
            self.E_ks.append(e_ks_c.copy())

            # Read only the first field for SHG
            efield=self.read_Efield(data_p_and_j,self.RT_step,1)
            try:
                efield2=self.read_Efield(data_p_and_j,self.RT_step,2)
            except:
                efield2=efield
            self.Efield.append(efield.copy())
            self.Efield2.append(efield2.copy())

    def get_info(self):
        """
        Provide information on the attributes of the class
        """
        print('YamboNLDBParser variables structure')
        s="\n * * * ndb.Nonlinear db data * * * \n\n"
        s+="Gauge         : "+str(self.Gauge)+"\n"
        s+="NE_steps      : "+str(self.NE_steps)+"\n"
        s+="RT_step       : "+str(self.RT_step/FsToAu)+" [fs] \n"
        s+="n_frequencies : "+str(self.n_frequencies)+"\n"   
        s+="n_angles      : "+str(self.n_angles)+"\n"   
        s+="NL_initial_versor   : "+str(self.NL_initial_versor)+"\n"   
        s+="NL_damping    : "+str(self.NL_damping*HaToeV)+"\n"
        s+="RT_bands      : "+str(self.RT_bands)+"\n"
        s+="NL_er         : "+str(self.NL_er*HaToeV)+" [eV] \n"
        s+="Sencond Order : "+str(self.l_force_SndOrd)+"\n" 
        s+="Use Dipoles   : "+str(self.l_use_DIPOLES)+"\n"
        s+="QP_ng_SH      : "+str(self.QP_ng_SH)+"\n"
        s+="QP_ng_Sx      : "+str(self.QP_ng_Sx)+"\n"  
        s+="RAD_LifeTime  : "+str(self.RAD_LifeTime/FsToAu)+" [fs] \n" 
        s+="Integrator    : "+str(self.Integrator)+"\n"
        s+="Correlation   : "+str(self.Correlation)+"\n"
        for efield in self.Efield:
            if efield["name"] == "none":
                continue
            s+="\nEfield name         : "+str(efield["name"])+"\n"
            s+="Efield versor       : "+str(efield["versor"])+"\n"  
            s+="Efield Intesity     : "+str(efield["intensity"])+"\n"
            s+="Efield Damping      : "+str(efield["damping"])+"\n"
            s+="Efield Freq range   : "+str(efield["freq_range"]*HaToeV)+" [ev] \n"
            s+="Efield Initial time : "+str(efield["initial_time"]/FsToAu)+" [fs] \n"
        print(s)

    def get_time(self,convert_to_fs=True):
        """
        Get the time points of the real-time propagation.
        
        Args:
            convert_to_fs (:py:class:`boolean`) : if True, convert the time points from atomic units to femto-seconds
        
        Returns:
            :py:class:`array` : array with the time points of the real-time propagation (in fs)
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
        time = self.get_time(convert_to_fs)
        pol = self.Polarization[run_index]
        Plot_3dArray(time,pol,xlim=xlim,label='Polarization')
    
    def plot_current(self,convert_to_fs=True,xlim=None,run_index=0):
        """
        Plot the current of the system as a function of time.

        Args:
            convert_to_fs (:py:class:`boolean`) : if True, convert the time points from atomic units to femto-seconds
            xlim (:py:class:`tuple`) : tuple with the limits of the x-axis        
            run_index (:py:class:`int`) : index of the run for which to plot the current  
        """
        time = self.get_time(convert_to_fs)
        curr = self.Current[run_index]
        Plot_3dArray(time,curr,xlim=xlim,label='J')

