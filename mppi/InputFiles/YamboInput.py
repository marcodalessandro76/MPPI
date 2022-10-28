
"""
Class to create and manipulate the yambo input files.
The class is partially inspired from the YamboIn class of YamboPy. In this implementation
the input object inherit from dict, so all the standard methods for python dictionaries can
be used to modify the attribute of the input.
"""

from subprocess import Popen, PIPE
import os, re
from sys import exit

class YamboInput(dict):

    #Regular expressions
    _variaexp   = '([A-Za-z\_0-9]+(?:\_[A-Za-z]+)?)' #variables names
    _numexp     = '([+-]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)' #number
    _spacexp    = '(?:[ \t]+)?' #space
    _stringexp  = '["\']([a-zA-Z0-9_ ]+?)["\']' #string
    _arrayexp   = '%'+_spacexp+_variaexp+'\s+(?:\#.+)?((?:(?:\s|\.|[+-]?\d)+?\|)+)\s+([a-zA-Z]+)?' #arrays
    _complexexp = '\('+_spacexp+_numexp+_spacexp+','+_spacexp+_numexp+_spacexp+'\)' #complex numbers
    _runexp     = '([a-zA-Z0-9_]+)'
    # list of available runlevels to be stored in the arguments array.
    # Also the 'options' like RmTimeRev or DephCVonly were included in the _runlevels list but they have
    # been removed since otherwise it seems that these values are always included in the arguments list
    # after the parsing of the input file.
    _runlevels  = ['rim_cut','RIM_W','chi','em1s','bse','optics','bsk','bss','em1d','gw0','HF_and_locXC','setup',
                   'ppa','cohsex','life','collisions','negf','el_ph_scatt','el_el_scatt','excitons',
                   'wavefunction','fixsyms','QPDBs', 'QPDB_merge','RealTime','RT_X','RToccDos',
                   'RToccBnd','RToccEner','RToccTime','RTlifeBnd','amplitude','bzgrids','Random_Grid',
                   'gkkp','el_ph_corr','WRbsWF','Select_energy', 'RTDBs','photolum','kpts_map',
                   'RTtime','RToccupations','RTfitbands']

    def __init__(self,args='',folder='.',filename='yambo.in'):
        """
        Initalize the class
        """
        dict.__init__(self,args=args,folder=folder,filename=filename)

        if args != '': # if args is not empty call yambo to generate the filename input file
            workdir = os.getcwd()
            os.chdir(folder)
            os.system('rm -f %s'%filename)
            args+= ' -F %s'%filename # add -F filename so yambo generates filename with the chosen args
            yambo = Popen(args, stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True)
            yambo.wait()
            os.chdir(workdir)
            self.read_file(os.path.join(folder,filename))
        else: # otherwise directly read the filename input file
            self.read_file(os.path.join(folder,filename))

    def read_file(self,file):
        """
        Open filename and run parseInputFile to parse the input into the dictionary
        """
        try:
            yambofile = open(file,'r')
        except IOError:
            print('Could not read the file %s'%filename)
            print('ERROR: yambo did not create the input file or the file you are trying to read does not exist')
            print('command: %s'%self['args'])
            print('folder:  %s/'%self['folder'])
            exit()
        self.parseInputFile(yambofile.read())
        yambofile.close()

    def write(self,folder,filename,reformat=True):
        """
        Write the yambo input on file. If the args of the object is not empty and
        the reformat variable is True run yambo to recover the original format of
        the yambo input.
        """
        f = open(os.path.join(folder,filename),'w')
        f.write(self.convert_string())
        f.close()
        if self['args'] != '' and reformat:
            action = 'cd %s; %s -F %s'%(folder,self['args'],filename)
            os.system(action)

    def parseInputFile(self,file):
        """
        Read the arguments and variables from the input file
        """
        arguments = []
        variables = {}

        var_real     = re.findall(self._variaexp + self._spacexp + '='+ self._spacexp +
                                  self._numexp + self._spacexp + '([A-Za-z]+)?',file)
        var_string   = re.findall(self._variaexp + self._spacexp + '='+ self._spacexp + self._stringexp, file)
        var_array    = re.findall(self._arrayexp,file)
        var_complex  = re.findall(self._variaexp + self._spacexp + '='+ self._spacexp +
                                  self._complexexp + self._spacexp + '([A-Za-z]+)?', file)
        var_runlevel = re.findall(self._runexp + self._spacexp, file)

        def clean(a):
            """
            clean the variables according to the type of data
            """
            a = a.strip()
            if a.replace('.','',1).isdigit():
                if "." in a: return float(a)
                else:        return int(a)
            return a

        # Determination of the arguments
        for key in self._runlevels:
            if key in var_runlevel:
                arguments.append(key)

        #float variables
        for var in var_real:
            name, value, unit = var
            variables[name] = [float(value),unit]

        #string variables
        for var in var_string:
            name, string = var
            variables[name] = string

        #complex variables
        for var in var_complex:
            name, real, imag, unit = var
            variables[name] = [complex(float(real),float(imag)),unit]

        #array variables
        for var in var_array:
            name, array, unit = var
            array = [clean(val) for val in array.split('|')[:-1]]
            variables[name] = [array,unit]

        self['arguments'] = arguments
        self['variables'] = variables

    def convert_string(self):
        """
        Convert the input object into a string
        """
        s  = ""
        s += "\n".join(self['arguments'])+'\n'

        for key,value in self['variables'].items():
            if type(value)==bytes or type(value)==str:
                s+= "%s = %10s\n"%(key,"'%s'"%value)
                continue
            if type(value[0])==float:
                val, unit = value
                if val > 1e-6:
                    s+="%s = %lf %s\n"%(key,val,unit)
                else:
                    s+="%s = %e %s\n"%(key,val,unit)
                continue
            if type(value[0])==int:
                val, unit = value
                s+="%s = %d %s\n"%(key,val,unit)
                continue
            if type(value[0])==list:
                array, unit = value
                if type(array[0])==list:
                    s+='%% %s\n'%key
                    for l in array:
                        s+="%s \n"%(" | ".join(map(str,l))+' | ')
                    s+='%s'%unit
                    s+='%\n'
                else:
                    s+="%% %s\n %s %s \n%%\n"%(key," | ".join(map(str,array))+' | ',unit)
                continue
            if type(value[0])==str:
                array = value
                s+="%% %s\n %s \n%%\n"%(key," | ".join(map(lambda x: "'%s'"%x.replace("'","").replace("\"",""),array))+' | ')
                continue
            if type(value[0])==complex:
                value, unit = value
                s+="%s = (%lf,%lf) %s\n"%(key,value.real,value.imag,unit)
                continue
            raise ValueError( "Unknown type %s for variable: %s" %( type(value), key) )
        return s

    # Set methods useful for Yambo inputs

    def set_extendOut(self):
        """
        Activate the ExtendOut option to print all the variable in the output file.
        """
        self['arguments'].append('ExtendOut')

    def set_kRange(self,first_k,last_k):
        """
        Set the the kpoint interval in the variable QPkpoint.
        """
        bands = self['variables']['QPkrange'][0][2:4]
        kpoint_bands = [first_k,last_k] + bands
        self['variables']['QPkrange'] = [kpoint_bands,'']

    def set_bandRange(self,first_band,last_band):
        """
        Set the the band interval in the variable QPkpoint.
        """
        kpoint = self['variables']['QPkrange'][0][0:2]
        kpoint_bands = kpoint + [first_band,last_band]
        self['variables']['QPkrange'] = [kpoint_bands,'']

    def set_GbndRange(self,first_band,last_band):
        """
        Set the the band interval in the variable GbndRnge.
        This variable specifies the range of bands entering in the sum over states
        in the correlation part of the self energy.
        """
        bands = [first_band,last_band]
        self['variables']['GbndRnge'] = [bands,'']

    def set_BndsRnXp(self,first_band,last_band):
        """
        Set the the band interval in the variable BndsRnXp.
        This variable specifies the number of bands entering in the sum over states
        in the RPA response function.
        """
        bands = [first_band,last_band]
        self['variables']['BndsRnXp'] = [bands,'']

    def activate_RIM_W(self):
        """
        Activate the RIM_W option to perform to the random integration method
        on the effective potential.
        """
        if 'RIM_W' not in self['arguments'] :
            self['arguments'].append('RIM_W')

    def deactivate_RIM_W(self):
        """
        Remove the RIM_W option from the runlevel list.
        """
        if 'RIM_W' in self['arguments'] :
            self['arguments'].remove('RIM_W')

    # Set methods useful for Yambo_rt inputs

    def set_rt_field(self,index=1,int=1e3,int_units='kWLm2',fwhm=100.,fwhm_units='fs',
                    freq=1.5,freq_units='eV',kind='QSSIN',polarization='linear',
                    direction=[1.,0.,0.],direction_circ=[0.,1.,0.],tstart=0.,tstart_units='fs'):
        """
        Set the parameters of the field.
        The index parameter is an integer that defines the name of the Field$index.
        Useful to set more than one field

        """
        field_name = 'Field'+str(index)
        self['variables'][field_name+'_Int'] = [int,int_units]
        self['variables'][field_name+'_FWHM'] = [fwhm,fwhm_units]
        self['variables'][field_name+'_Freq'] = [[freq,freq],freq_units]
        self['variables'][field_name+'_kind'] = kind
        self['variables'][field_name+'_pol'] = polarization
        self['variables'][field_name+'_Dir'] = [direction,'']
        self['variables'][field_name+'_Dir_circ'] = [direction_circ,'']
        self['variables'][field_name+'_Tstart'] = [tstart,tstart_units]

    def set_rt_bands(self,bands=None,scissor=0.,damping_valence=0.05,damping_conduction=0.05):
        """
        Set the bands, the scissor and the damping parameters for the RT analysis
        """
        if bands is not None:
            self['variables']['RTBands'] = [bands,'']
        self['variables']['GfnQP_E'] = [[scissor, 1.0, 1.0], '']
        self['variables']['GfnQP_Wv'] = [[damping_valence, 0.0, 0.0], '']
        self['variables']['GfnQP_Wc'] = [[damping_conduction, 0.0, 0.0], '']

    def set_rt_simulationTimes(self,time_step=10,step_units='as',sim_time=1000.,time_units='fs',
                                io_time = [1.,5.,1.],io_units='fs'):
        """
        Set the time parameters of the simulation
        """
        self['variables']['RTstep'] = [time_step,step_units]
        self['variables']['NETime'] = [sim_time,time_units]
        self['variables']['IOtime'] = [io_time,io_units]

    def set_rt_cpu(self,k=1,b=1,q=1,qp=1):
        """
        Set the parallelization roles of the run
        """
        self['variables']['RT_CPU'] = '%s.%s.%s.%s'%(k,b,qp,q)
        self['variables']['RT_ROLEs'] = 'k.b.qp.q'

    # Set methods useful for Ypp inputs

    def removeTimeReversal(self):
        """
        Remove the time reversal symmetry
        """
        self['arguments'].append('RmTimeRev')

    def set_ypp_extFields(self, Efield1 = [1.,0.,0.], Efield2 = [0.,1.,0.]):
        """
        Set the direction of the external electric field(s)
        """
        self['variables']['Efield1'] = [Efield1,'']
        self['variables']['Efield2'] = [Efield2,'']
