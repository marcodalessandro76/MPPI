
# YamboIn class of YamboPy, Modified for compliance with python 3

from subprocess import Popen, PIPE
import os, json, re
from sys import exit

class YamboIn():
    """
    Class to read, write, create and manipulate yambo input files with python.

    Examples of use:

    Initialize an empty input file:

        .. code-block:: python

            y = YamboIn(filename='somefile.in')
            print y

    Call yambo to initialize the input file with variables according to the runlevel,
    parse the input file and store the variables:

        .. code-block:: python

            y = YamboIn('yambo -o c',folder='ip')
            print y

    If the argument ``args`` was used then the filename should be left as ``yambo.in`` because that's the default input filename that yambo will create.

    Call ypp to initialize the input file:

        .. code-block:: python

            y = YamboIn('yyp -e w'args=,filename='ypp.in')
            print y

    **Arguments:**

        ``args``:     if specified yambopy will run yambo, read the generated input file and initialize the class with those variables.

        ``folder``:   the folder where the SAVE directory is located

        ``vim``:      if yambo is compiled using vim as an editor this variable should be set as True because then `yambopy` will close vim.
        In newer versions an argument for yambo '-Q' tells it to not call vim

        ``filename``: the name of the input file to be read

    """
    #Regular expressions
    _variaexp   = '([A-Za-z\_0-9]+(?:\_[A-Za-z]+)?)' #variables names
    _numexp     = '([+-]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)' #number
    _spacexp    = '(?:[ \t]+)?' #space
    _stringexp  = '["\']([a-zA-Z0-9_ ]+?)["\']' #string
    _arrayexp   = '%'+_spacexp+_variaexp+'\s+(?:\#.+)?((?:(?:\s|\.|[+-]?\d)+?\|)+)\s+([a-zA-Z]+)?' #arrays
    _complexexp = '\('+_spacexp+_numexp+_spacexp+','+_spacexp+_numexp+_spacexp+'\)' #complex numbers
    _runexp     = '([a-zA-Z0-9_]+)' #runlevels
    # list of available runlevels to be stored in the arguments array
    _runlevels  = ['rim_cut','chi','em1s','bse','optics','bsk','bss',
                   'em1d','gw0','HF_and_locXC','setup','ppa','cohsex','life',
                   'collisions','negf','el_ph_scatt','el_el_scatt','excitons','wavefunction','fixsyms',
                   'QPDBs', 'QPDB_merge','RealTime','RT_X','RToccDos','RToccBnd','RToccEner',
                   'RToccTime','RTlifeBnd','amplitude','bzgrids','Random_Grid','gkkp','el_ph_corr','WRbsWF','Select_energy', 'RTDBs','photolum','kpts_map',
                   'RTtime','RToccupations','RTfitbands']

    def __init__(self,args='',folder='.',vim=False,filename='yambo.in'):
        """
        Initalize the class
        """
        self.folder = folder
        self.yamboargs = args

        #the type of the variables is determined from the type of variable in this dictionary
        self.variables = {} #here we will store the values of the variables
        self.arguments = [] #here we will store the arguments

        # if we initalize the class with arguments we call yambo to generate the input file
        if args != '':
            workdir = os.getcwd()
            os.chdir(folder)
            os.system('rm -f %s'%filename)
            yambo = Popen(args, stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True)
            # if yambo calls vim we have to close it. We just want the generic input file
            # that yambo generates.
            #if vim: yambo.stdin.write(":q!\n")
            yambo.wait()
            os.chdir(workdir)
            self.read_file(filename="%s/%s"%(folder,filename))
        else:
            if filename:
                self.read_file(filename="%s/%s"%(folder,filename))

    def __getitem__(self,key):
        """ Get the value of a variable in the input file
        """
        return self.variables[key]

    def __setitem__(self,key,value):
        """ Set the value of a variable in the input file
        """
        #if the units are not specified, add them
        if type(value) == list and str not in map(type,value):
            value = [value,'']
        if type(value) in [int,float,complex]:
            value = [value,'']
        self.variables[key] = value

    def __delitem__(self,key):
        """ Remove a keyword from the dicitonary
        """
        del self.variables[key]

    def read_file(self,filename='yambo.in'):
        """ Read the variables from a file
        """
        try:
            yambofile = open(filename,"r")
        except IOError:
            print('Could not read the file %s'%filename)
            print('Something is wrong, yambo did not create the input file. Or the file you are trying to read does not exist')
            print('command: %s'%self.yamboargs)
            print('folder:  %s/'%self.folder)
            exit()
        inputfile = self.read_string(yambofile.read())
        yambofile.close()

    def add_dict(self,variables):
        """
        Add a dictionary containing variables to the current inputfile
        """
        self.variables.update(variables)

    def read_string(self,inputfile):
        """
        Read the input variables from a string
        """
        var_real     = re.findall(self._variaexp + self._spacexp + '='+ self._spacexp +
                                  self._numexp + self._spacexp + '([A-Za-z]+)?',inputfile)
        var_string   = re.findall(self._variaexp + self._spacexp + '='+ self._spacexp + self._stringexp, inputfile)
        var_array    = re.findall(self._arrayexp,inputfile)
        var_complex  = re.findall(self._variaexp + self._spacexp + '='+ self._spacexp + self._complexexp + self._spacexp + '([A-Za-z]+)?', inputfile)
        var_runlevel = re.findall(self._runexp + self._spacexp, inputfile)

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
                self.arguments.append(key)

        #float variables
        for var in var_real:
            name, value, unit = var
            self[name] = [float(value),unit]

        #string variables
        for var in var_string:
            name, string = var
            self[name] = string

        #complex variables
        for var in var_complex:
            name, real, imag, unit = var
            self[name] = [complex(float(real),float(imag)),unit]

        #array variables
        for var in var_array:
            name, array, unit = var
            array = [clean(val) for val in array.split('|')[:-1]]
            self[name] = [array,unit]

        return {"arguments": self.arguments, "variables": self.variables}

    def write(self,filename='yambo.in'):
        """
        Write a yambo input file
        """
        f = open(filename,"w")
        f.write(str(self))
        f.close()

    def pack(self,filename):
        """
        Pack all the data of this structure in a `.json` file
        """
        f = open(filename,'w')
        json.dump(f,[self.arguments,self.real,self.string,self.complex,self.array],indent=5)
        f.close()

    def __str__(self):
        """
        Returns the input file as a string
        """
        s  = ""

        #arguments
        s += "\n".join(self.arguments)+'\n'

        for key,value in self.variables.items():
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
