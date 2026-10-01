
"""
Class to create and manipulate the yambo input files.
The class is partially inspired from the YamboIn class of YamboPy. In this implementation
the input object inherit from dict, so all the standard methods for python dictionaries can
be used to modify the attribute of the input.
"""

from subprocess import run
import os, re

# regular expressions used to parse the lines of a yambo input file
_comment_exp = re.compile(r'''("[^"]*"|'[^']*')|#.*''')
_runlevel_exp = re.compile(r'[A-Za-z_][A-Za-z0-9_]*')
_number_exp = re.compile(r'[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?')
_complex_exp = re.compile(r'\(\s*(\S+)\s*,\s*(\S+)\s*\)\s*(.*)')

def strip_comment(line):
    """
    Remove the comment (introduced by #) from a line of the input file, preserving
    the # characters inside quoted strings.
    """
    return _comment_exp.sub(lambda m: m.group(1) or '', line).strip()

def convert_value(token):
    """
    Convert a token of the input file into an int, a float or a string (without quotes).
    """
    token = token.strip()
    if token[:1] in ('"',"'"):
        return token.strip('"\'')
    if _number_exp.fullmatch(token):
        if re.fullmatch(r'[+-]?\d+',token):
            return int(token)
        return float(token.replace('d','e').replace('D','e'))
    return token

def parse_scalar(rhs):
    """
    Parse the right hand side of a `name = value [units]` line.

    Returns:
        the string value (for quoted strings) or the list [value,units] (for numbers and complex numbers)

    """
    rhs = rhs.strip()
    if rhs[:1] in ('"',"'"):
        return rhs[1:rhs.find(rhs[0],1)]
    complex_match = _complex_exp.fullmatch(rhs)
    if complex_match:
        real, imag, units = complex_match.groups()
        return [complex(float(real),float(imag)),units.strip()]
    value, _, units = rhs.partition(' ')
    return [convert_value(value),units.strip()]

def parse_array_row(row):
    """
    Parse a row of an array block, with the structure `v1 | v2 | ... | [units]`.

    Returns:
        :py:class:`tuple` : the list of the values of the row and the units (empty string if not given)

    """
    fields = row.split('|')
    values = [convert_value(v) for v in fields[:-1]]
    return values, fields[-1].strip()

def format_number(value):
    """
    Format a number of the input file. Floats are written with the shortest representation
    that preserves their value.
    """
    if isinstance(value,complex):
        return '( %s , %s )'%(format_number(value.real),format_number(value.imag))
    if isinstance(value,str):
        return '"%s"'%value
    return str(value)

def format_scalar(name, value, units):
    """
    Format a scalar variable as `name= value units`.
    """
    return ('%s= %s %s'%(name,format_number(value),units)).rstrip()

def format_array(name, values, units):
    """
    Format an array variable as a yambo array block. values can be a list (a single row) or a
    list of lists (a matrix). The units are written at the end of the last row.
    """
    rows = values if len(values) > 0 and isinstance(values[0],list) else [values]
    lines = ['%% %s'%name]
    for row in rows:
        if len(row) > 0:
            lines.append(' | '.join(format_number(v) for v in row) + ' |')
    if units != '' and len(lines) > 1:
        lines[-1] += ' ' + units
    lines.append('%')
    return '\n'.join(lines)

def format_variable(name, value):
    """
    Convert a variable of the input object into the yambo syntax. The variable can be
    a string, a list [value,units] (where value is a number, a complex number or a list) or a
    list of strings (an array of strings without units).
    """
    if isinstance(value,str):
        return '%s= "%s"'%(name,value)
    if isinstance(value,(list,tuple)) and len(value) == 2 and isinstance(value[1],str) \
        and not isinstance(value[0],str):
        val, units = value
        if isinstance(val,(list,tuple)):
            return format_array(name,list(val),units)
        return format_scalar(name,val,units)
    if isinstance(value,(list,tuple)) and all(isinstance(v,str) for v in value):
        return format_array(name,list(value),'')
    raise ValueError('Unknown type %s for variable: %s'%(type(value),name))

class YamboInput(dict):
    """
    Class to create and manipulate the input files of yambo (and of the other executables of the
    package, like ypp, yambo_rt, yambo_nl, ...).

    The object is a dictionary with the keys:

    * `args`, `folder`, `filename` : the command used to generate the input, the folder and the name of
      the input file
    * `arguments` : list with the active runlevels and flags (e.g. ['HF_and_locXC'])
    * `variables` : dictionary with the variables. A string variable is stored as a string, while
      numbers, complex numbers and arrays are stored as a list [value,units], for instance
      `'EXXRLvcs' : [5985,'RL']` and `'QPkrange' : [[1,32,1,8],'']`. An array with more rows is stored
      as a list of lists

    Args:
        args (:py:class:`string`) : command line used to generate the input (e.g. 'yambo -x -V rl').
            If not empty yambo is executed in the folder (that must contain the SAVE folder) to write the
            input file. If empty an existing input file is read
        folder (:py:class:`string`) : folder of the input file
        filename (:py:class:`string`) : name of the input file

    """

    def __init__(self,args='',folder='.',filename='yambo.in'):
        dict.__init__(self,args=args,folder=folder,filename=filename)
        if args != '': # call yambo to generate the input file with the chosen args
            file = os.path.join(folder,filename)
            if os.path.isfile(file): os.remove(file)
            run('%s -F %s'%(args,filename),shell=True,cwd=folder,capture_output=True)
        self.read_file(os.path.join(folder,filename))

    def read_file(self,file):
        """
        Open filename and run parseInputFile to parse the input into the dictionary
        """
        try:
            with open(file,'r') as yambofile:
                self.parseInputFile(yambofile.read())
        except IOError:
            raise IOError('Could not read the file %s. Yambo did not create the input file or the file '
                          'you are trying to read does not exist (command: %s, folder: %s)'
                          %(file,self['args'],self['folder']))

    def write(self,folder,filename,reformat=True):
        """
        Write the yambo input on file. If the args of the object is not empty and
        the reformat variable is True run yambo to recover the original format of
        the yambo input.
        """
        with open(os.path.join(folder,filename),'w') as f:
            f.write(self.convert_string())
        if self['args'] != '' and reformat:
            run('%s -F %s'%(self['args'],filename),shell=True,cwd=folder,capture_output=True)

    def parseInputFile(self,file):
        """
        Read the arguments and variables from the content of an input file. The lines are
        classified as:

        * array blocks, that start with `% name` and end with a line with `%`
        * variables, with the structure `name = value [units]`
        * runlevels and flags, given as a single word

        Comments (introduced by #) are ignored.

        Args:
            file (:py:class:`string`) : the content of the input file

        """
        arguments = []
        variables = {}
        lines = iter(file.splitlines())
        for line in lines:
            line = strip_comment(line)
            if line == '':
                continue
            if line.startswith('%'):
                name = line[1:].strip()
                rows, units = [], ''
                for row in lines:
                    row = strip_comment(row)
                    if row == '%': break
                    if row == '': continue
                    values, row_units = parse_array_row(row)
                    rows.append(values)
                    units = row_units or units
                values = rows[0] if len(rows) == 1 else rows
                variables[name] = [values,units]
            elif '=' in line:
                name, rhs = line.split('=',1)
                variables[name.strip()] = parse_scalar(rhs)
            elif _runlevel_exp.fullmatch(line):
                arguments.append(line)
        self['arguments'] = arguments
        self['variables'] = variables

    def convert_string(self):
        """
        Convert the input object into a string with the syntax of the yambo input files
        """
        lines = list(self['arguments'])
        lines += [format_variable(name,value) for name,value in self['variables'].items()]
        return '\n'.join(lines) + '\n'

    # Set methods useful for Yambo inputs

    def set_array_variables(self,units='',**kwargs):
        """
        Add to the `variables` key of the input dictionary
        the elements kwargs[key] = [kwargs[value],units] for all the
        elements of the kwargs provided as input.

        Args:
            units (:py:class:`string`) : the units associated to (all the) variables
            kwargs : variable(s) added in the form name = value. All the
                variables must have the same units

        """
        for name,value in kwargs.items():
            self['variables'][name] = [value,units]

    def set_scalar_variables(self,**kwargs):
        """
        Add to the `variables` key of the input dictionary
        the elements kwargs[key] = kwargs[value] for all the
        elements of the kwargs provided as input.

        Args:
            kwargs : variable(s) added in the form name = value

        """
        for name,value in kwargs.items():
            self['variables'][name] = value

    def set_extendOut(self):
        """
        Activate the ExtendOut option to print all the variable in the output file.
        """
        if 'ExtendOut' not in self['arguments']:
            self['arguments'].append('ExtendOut')

    def set_kRange(self,first_k,last_k):
        """
        Set the the kpoint interval in the variable QPkrange.
        """
        bands = self['variables']['QPkrange'][0][2:4]
        kpoint_bands = [first_k,last_k] + bands
        self['variables']['QPkrange'] = [kpoint_bands,'']

    def set_bandRange(self,first_band,last_band):
        """
        Set the the band interval in the variable QPkrange.
        """
        kpoint = self['variables']['QPkrange'][0][0:2]
        kpoint_bands = kpoint + [first_band,last_band]
        self['variables']['QPkrange'] = [kpoint_bands,'']

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

    # Set methods useful for yambo_rt inputs

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

    def set_rt_bands(self,bands=None,scissor=0.,stretch_con=1.,stretch_val=1.,damping_valence=0.05,damping_conduction=0.05):
        """
        Set the bands, the scissor and the damping parameters for the RT analysis
        """
        if bands is not None:
            self['variables']['RTBands'] = [bands,'']
        self['variables']['GfnQP_E'] = [[scissor, stretch_con, stretch_val], '']
        self['variables']['GfnQP_Wv'] = [[damping_valence, 0.0, 0.0], '']
        self['variables']['GfnQP_Wc'] = [[damping_conduction, 0.0, 0.0], '']

    def set_rt_simulationTimes(self,time_step=10,step_units='as',sim_time=1000.,time_units='fs',
                                io_time = [1.,5.,1.],io_cache_time = [1.,1.],io_units='fs'):
        """
        Set the time parameters of the simulation
        """
        self['variables']['RTstep'] = [time_step,step_units]
        self['variables']['NETime'] = [sim_time,time_units]
        self['variables']['IOtime'] = [io_time,io_units]
        self['variables']['IOCachetime'] = [io_cache_time,io_units]

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
        if 'RmTimeRev' not in self['arguments']:
            self['arguments'].append('RmTimeRev')

    def set_ypp_extFields(self, Efield1 = [1.,0.,0.], Efield2 = None):
        """
        Set the direction of the external electric field(s). The second field (useful for the
        circular polarization) is added only is the value is not `None`.
        """
        self['variables']['Efield1'] = [Efield1,'']
        if Efield2 is not None :
            self['variables']['Efield2'] = [Efield2,'']
