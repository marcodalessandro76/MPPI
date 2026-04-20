"""
This module contains a (almost exact and partially simplified )vcopy of the nl_analysis.py module of the yambopy project. 
This material is kept here (at least for the moment) for comparison with the actual implementation of the non-linear optics analysis 
in the mppi project.
The module can be loaded in the notebook as follows

>>> from mppi.Utilities import NLanalysisYamboPy as NLYP

>>> NLYP.Xn_from_sine

"""

# Copyright (c) 2023-2026, Claudio Attaccalite,
#                          Myrta Gruening
# All rights reserved.
#
# This file is part of the yambopy project
# Calculate linear response from real-time calculations (yambo_nl)
# Modified to allow generalisation
#


import numpy as np
from mppi.Utilities import Constants as C
from tqdm import tqdm
from scipy.optimize import least_squares
import os, sys
import itertools
from abc import ABC,abstractmethod

#####################################################################################
SVCMm12VMm1 =29.98*np.power(10,3,dtype=np.double) #  from [statV/cm] to [V/m]
AU2VMm1     =5.14220632*np.power(10,11,dtype=np.double)
######################################################################################

def Divide_by_the_Field(efield,order):
    
    if efield['name']=='SIN' or efield['name']=='SOFTSIN':
        if order !=0:
            divide_by_field=np.power(-2.0*1.0j/efield['amplitude'],order,dtype=np.cdouble)
        elif order==0:
            divide_by_field=-4.0/np.power(efield['amplitude'],2.0,dtype=np.cdouble)
    else:
        raise ValueError("Electric field not implemented in Divide_by_the_Field!")

    return divide_by_field

#
#
# Template class
#
class Xn_from_signal(ABC):
    
    def __init__(self,nldb,X_order=4,T_range=[-1, -1],nsamp=-1,solver='full',tol=1e-10,debug_mode=False): 
        self.time = nldb.IO_TIME_points # Time series
        self.T_step =self.time[1] - self.time[0]     # Time step of the simulation
        self.T_deph =12.0/nldb.NL_damping # NL_damping is the energy damping in au. Why 12?
        self.efield = nldb.Efield[0] # External field of the first run
        self.pumps = nldb.Efield_general[1:] 
        self.efields = nldb.Efield
        self.n_runs = len(nldb.Polarization)     # Number of external laser frequencies
        self.polarization = nldb.Polarization     # Array of polarizations for each laser frequency
        self.X_order = X_order
        SOLVE_MODES = ['', 'full', 'lstsq', 'lstsq_opt', 'svd']
        if solver not in SOLVE_MODES: 
            raise ValueError("Invalid solver mode. Expected one of: %s" % SOLVE_MODES)
        self.solver = solver
        self.nsamp = nsamp
        self.T_urange = T_range 
        self.freqs = np.array([efield["freq_range"][0] for efield in nldb.Efield], dtype=np.double)     # Energies in au
        self.out_dim = 0
        self.tol = tol
        self.debug_mode = debug_mode
        super().__init__()
        
    def __str__(self): 
        """
        Print info of the class
        """
        s="\n * * *  Xn from signal class  * * * \n\n"
        s+="Max time: "+str(self.time[-1])+"\n"
        s+="Time step : "+str(self.T_step)+"\n"
        s+="Type Efield    : "+str(self.efield["name"])+"\n"
        s+="Number of runs   : "+str(self.n_runs)+"\n"
        s+="Max harmonic order   : "+str(self.X_order)+"\n"
        s+="Solver           : "+str(self.solver)+"\n"
        if self.nsamp > 0:
            s+="Sampling points  : "+str(self.nsamp)+"\n"
        if self.T_urange!=[-1, -1]:
            s+="User time range      : "+str(self.T_urange)+" [au] \n"
        s+="Frequency range: ["+str(self.freqs[0])+","+str(self.freqs[-1])+"] [au] \n"
        return s
    
    @abstractmethod
    def set_defaults(self):
        pass

    @abstractmethod
    def set_sampling(self,ifrq): 
       pass
   
    @abstractmethod
    def define_matrix(self,samp,ifrq):
        pass

    @abstractmethod
    def update_time_range(self):
        pass

    @abstractmethod
    def get_Unit_of_Measure(self,i_order):
        pass

    def get_sampling(self,T_range,idir,ifrq): 
        i_t_start = int(np.round( T_range[0] / self.T_step)) 
        i_deltaT  = int(np.round((T_range[1]-T_range[0])/self.T_step)/self.nsamp)
        i_t = i_t_start + i_deltaT * np.arange(self.nsamp)
        T_i = i_t*self.T_step - self.efield["initial_time"] 
        S_i = np.array([self.polarization[ifrq][idir,i] for i in i_t]) 
        return T_i,S_i
    
    def solve_lin_system(self,mat,samp,init=None):
        mat_dim = int(mat.shape[1])
        out=np.zeros(mat_dim,dtype=np.cdouble)
        if self.solver=="full" and ((not self.IsSquare(mat)) or (not self.IsWellConditioned(mat))):
            print(f'WARNING: solver changed to least square since square:{self.IsSquare(mat)} and well-conditioned:{self.IsWellConditioned(mat)}')
            self.solver = "lstsq"
        if self.solver=="full":
            out = np.linalg.solve(mat,samp)
        if self.solver=="lstsq":
            out = np.linalg.lstsq(mat,samp,rcond=self.tol)[0]
        if self.solver=="svd":
            inv = np.linalg.pinv(mat,rcond=self.tol)
            for i_n in range(mat_dim):
                out[i_n]=out[i_n]+np.sum(inv[i_n,:]*samp[:])
        if self.solver=="lstsq_opt":
            if(init is None):
                x0_cmplx = np.linalg.lstsq(mat, samp, rcond=tol)[0]
            else:
                x0_cmplx = init
            x0 = np.concatenate((x0_cmplx.real, x0_cmplx.imag))
            res = least_squares(residuals_func, x0, ftol=1e-11,gtol=1e-11,xtol=1e-11,verbose=0,x_scale='jac',args=(mat,samp))
            out = res.x[0:int(res.x.size/2)] + 1j * res.x[int(res.x.size/2):res.x.size]
        return out
    
    def perform_analysis(self):
        _ = self.set_defaults()
        out = np.zeros((self.out_dim, self.n_runs, 3), dtype=np.cdouble) 
        for i_f in tqdm(range(self.n_runs)):
            T_r = self.set_sampling(i_f)
            for i_d in range(3):
                samp_time, samp_sig= self.get_sampling(T_r,i_d,i_f)
                matrix = self.define_matrix(samp_time,i_f)
                raw = self.solve_lin_system(matrix,samp_sig)
                out[:, i_f, i_d] = raw[:self.out_dim]
                if self.debug_mode:
                    print(f"Freq #{i_f}, direction: {i_d}:")
                    print("***Sampling:")
                    print(samp_time, samp_sig)
                    print("***Matrix:")
                    print(matrix)
                    print("***Solution:")
                    print(raw)
        return out

    def get_Unit_of_Measure(self,i_order): # not sure if this is a general or specific method - let it here for the moment
        ratio = SVCMm12VMm1 / AU2VMm1 # From AU to statVolt/cm ...is there a better way than this?
        if i_order == 0:
            return np.power(ratio, 1, dtype=np.float64)
        return np.power(ratio, i_order - 1, dtype=np.float64)

    @abstractmethod
    def output_analysis(self,out,to_file):
        pass

    @abstractmethod
    def reconstruct_signal(self,out,to_file):
        pass

    def append_runinfo(self,T):
        s = "# Field details:"
        s+= "# Type field            "+str(self.efield["name"])+"\n"
        s+= "# Field intensity       "+str(self.efield["intensity"])+"\n"
        s+= "# Field versor          "+str(self.efield["versor"])+"\n"
        s+= "# Analysis details:"
        s+="# Max harmonic order   : "+str(self.X_order)+"\n"
        s+="# Solver               : "+str(self.solver)+"\n"
        s+="# Sampling points      : "+str(self.nsamp)+"\n"
        s+="# Start sampling time  : "+str(T/C.FsToAu)+" [fs] \n"
        return s
        

### some maths auxiliary methods
    def IsSquare(self,m):
        return m.shape[0] == m.shape[1]

    def IsWellConditioned(self,m): # with this I am trying to avoid inverting a matrix
        return np.linalg.cond(m) < 1/sys.float_info.epsilon

    def residuals_func(x,M,S_i):
        x_cmplx=x[0:int(x.size/2)] + 1j * x[int(x.size/2):x.size]
        return np.linalg.norm(np.dot(M, x_cmplx) - S_i)
    
#
#
# Derived class for monochromatic signal
#    
class Xn_from_sine(Xn_from_signal):
        
        def set_defaults(self):
            EFIELDS = ["SIN","SOFTSIN"]
            if self.efield["name"] not in EFIELDS:
                raise ValueError(f"Invalid electric field for frequency mixing analysis. Expected one of: {EFIELDS}")
            for i_n in range(len(self.pumps)):
                if self.pumps[i_n]["name"] != 'none':
                    raise ValueError("This analysis is for one monochromatic field only.")
            self.out_dim = self.X_order + 1  # Why not 2*X_order + 1? Only positive frequencies? And in case which is meaning of the zero order term? 
            if self.nsamp <= 0: 
                self.nsamp = 2*self.X_order + 1
            return

        def set_sampling(self,ifrq):
            T_period = 2.0 * np.pi / self.freqs[ifrq] # From E= h/T, E=freqs and h=2*pi in au
            T_range = self.update_time_range(T_period) # set the time of a period of the signal
            return T_range
        
        def update_time_range(self,T_period): # not sure if this is a general or specific method - let it here for the moment
            T_range = self.T_urange
            if T_range[0] <= 0.0:
                T_range[0] = self.time[-1] - T_period
            
            if T_range[1] > 0.0:
                T_range[1] = T_range[0] + T_period
            else:
                T_range[1] = self.time[-1]                
            
            if T_range[1] > self.time[-1]:
                T_range[1] = self.time[-1]
                T_range[0] = T_range[1] - T_period
                print(f'User range redifined for frequency {self.freqs[ifrq]*C.HaToeV:.3e} [eV]')
                print(f"Time range: {T_range[0] / C.FsToAu:.3f} - {T_range[1] / C.FsToAu:.3f} [fs]")
            return T_range
        
        def define_matrix(self,T_i,ifrq):
            M_size = len(T_i)
            M = np.zeros((M_size, M_size), dtype=np.cdouble)
            M[:, 0] = 1.0
            W = self.freqs[ifrq]
            for i_n in range(1, self.X_order+1):
                exp_neg = np.exp(-1j * i_n* W * T_i, dtype=np.cdouble)
                exp_pos = np.exp(1j * i_n * W * T_i, dtype=np.cdouble)
                M[:, i_n] = exp_neg
                M[:, i_n +self.X_order] = exp_pos
            return M

        def output_analysis(self,out):
            for i_order in range(self.X_order + 1):
                T = 10000.0
                for i_f in range(self.n_runs):
                    out[i_order, i_f, :] *= Divide_by_the_Field(self.efields[i_f], i_order)
                    T_period = 2.0 * np.pi / self.freqs[i_f]
                    Trange = self.update_time_range(T_period)
                    T = min(T,Trange[0])
                out[i_order,:,:]*=self.get_Unit_of_Measure(i_order)
                run_info = self.append_runinfo(T)
                print(run_info)
            return (self.freqs, out)

        def reconstruct_signal(self,out,to_file=False):
            Seff = np.zeros((self.n_runs, 3, len(self.time)), dtype=np.cdouble)
            for i_f in tqdm(range(self.n_runs)):
                for i_d in range(3):
                    for i_order in range(self.X_order + 1):
                        freq_term = np.exp(-1j * i_order * self.freqs[i_f] * self.time)
                        Seff[i_f, i_d, :] += out[i_order, i_f, i_d] * freq_term
                        Seff[i_f, i_d, :] += np.conj(out[i_order, i_f, i_d]) * np.conj(freq_term)    
            for i_f in tqdm(range(self.n_runs)):
                values = np.column_stack((self.time/C.FsToAu, Seff[i_f, 0, :].real, Seff[i_f, 1, :].real, Seff[i_f, 2, :].real))
            return values


#
# Derived class for frequency-mixed signal
#    
class Xn_from_freqmix(Xn_from_signal):
    #
    def set_defaults(self):
        EFIELDS = ["SIN","SOFTSIN"]
        if self.efield["name"] not in EFIELDS:
            raise ValueError(f"Invalid electric field for frequency mixing analysis. Expected one of: {EFIELDS}")
        EFIELDS = ["","SIN","SOFTSIN"]
        if self.pumps[0]["name"] not in EFIELDS:
            raise ValueError(f"Invalid pump field for frequency mixing analysis. Expected one of: {EFIELDS}")
        self.pump_freq = 0.0
        if self.pumps[0]["name"] != "none":
            self.pump_freq = self.pumps[0]["freq_range"][0]
        if len(self.pumps)>1 and self.pumps[1]["name"] != "none":
            raise ValueError("Only mixing of two fields implemented.")
        if isinstance (self.X_order,int): # for frequency mixing it expects 2 orders, if one is given the same is given
            self.X_order = [self.X_order,self.X_order]
            if self.pump_freq == 0.0:
                self.X_order[1] = 0
        self.out_dim = (2*self.X_order[0] + 1)*(2*self.X_order[1] + 1)
    #
    def update_time_range(self):
        T_range = self.T_urange
        if T_range[0] <= 0.0:
            T_range[0] = self.T_deph
        if T_range[1] <= 0.0:
            T_range[1]=self.time[-1]
        return T_range
    #    
    def set_sampling(self,ifrq):
        samp_order = self.out_dim
        if self.nsamp == -1 or self.nsamp < samp_order:
            self.nsamp = int(samp_order**2)
        #
        T_range = self.update_time_range()
        return T_range

    def define_matrix(self,T_i,ifrq):
        NX,MX = self.X_order[:]
        W1 = self.freqs[ifrq]
        W2 = self.pump_freq
        M = np.zeros((self.nsamp, self.out_dim), dtype=np.cdouble)
        for i_t in range(self.nsamp):
            for i_c,(i_n,i_m) in enumerate(itertools.product(range(-NX, NX+1),range(-MX, MX+1))):
                M[i_t, i_c] = np.exp(-1j * (i_n*W1+i_m*W2) * T_i[i_t],dtype=np.cdouble)
        return M

    def output_analysis(self,out,to_file=True):
        #
        NX,MX = self.X_order[:]
        T, _ = self.update_time_range()
        run_info = self.append_runinfo(T)
        #
        for i_f in range(self.n_runs):
            for i_c,(i_n,i_m) in enumerate(itertools.product(range(-NX, NX+1),range(-MX, MX+1))):
                field_fac = 1.0
                if i_n !=0: field_fac *= Divide_by_the_Field(self.efields[i_f],abs(i_n)) #check this
                if self.pump_freq !=0:
                    if i_m !=0: field_fac *= Divide_by_the_Field(self.pumps[0],abs(i_m)) #check this! what is nldb.Efield2?
                out[i_c,i_f,:] *= field_fac
                out[i_c,i_f,:] *= self.get_Unit_of_Measure(abs(i_n)+abs(i_m))
        # if (to_file):
        #     for i_c,(i_n,i_m) in enumerate(itertools.product(range(-NX, NX+1),range(-MX, MX+1))):
        #         output_file = f'o{self.prefix}.YamboPy-X_probe_order_{i_n}_{i_m}'
        #         iorder = (i_n+i_m,i_n,i_m)
        #         header = "E[eV] " + " ".join([f"X{iorder}/Im_{d} X{iorder}/Re_{d}" for d in ('x','y','z')])
        #         values = np.column_stack((self.freqs * ha2ev, out[i_c, :, 0].imag, out[i_c, :, 0].real,
        #                 out[i_c, :, 1].imag, out[i_c, :, 1].real,
        #                 out[i_c, :, 2].imag, out[i_c, :, 2].real))
        #         np.savetxt(output_file, values, header=header, delimiter=' ', footer=f"Frequency mixing analysis with pump frequency {self.pump_freq*ha2ev} eV")
        #         outf = open(output_file, "a")
        #         outf.writelines(run_info)
        #         outf.close()
        else:
            print(run_info)
            return (self.freqs, out)

    def reconstruct_signal(self,out,to_file=True):
        NX,MX = self.X_order[:]
        Seff = np.zeros((self.n_runs, 3, len(self.time)), dtype=np.cdouble)
        for i_f in tqdm(range(self.n_runs)):
            for i_d in range(3):
                for i_c,(i_n,i_m) in enumerate(itertools.product(range(-NX, NX+1),range(-MX, MX+1))):
                    freq_term = np.exp(-1j * (i_n * self.freqs[i_f] + i_m*self.pump_freq) * self.time)
                    Seff[i_f, i_d, :] += out[i_c, i_f, i_d] * freq_term
        # if (to_file):
        #     for i_f in tqdm(range(self.n_runs)):
        #         values = np.column_stack((self.time / fs2aut, Seff[i_f, 0, :].real, Seff[i_f, 1, :].real, Seff[i_f, 2, :].real))
        #         output_file = f'o{self.prefix}.YamboPy-pol_reconstructed_F{i_f + 1}'
        #         header="[fs] Px Py Pz"
        #         np.savetxt(output_file, values, header=header, delimiter=' ', footer="Reconstructed signal")
        else:
            return values
    