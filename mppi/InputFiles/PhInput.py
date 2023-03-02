"""
This module manages the input file for ph.x computations of QuantumESPRESSO.
The input can be created from scratch or can be initialized from an existing input file.
"""
import os
from copy import deepcopy

def fortran_bool(boolean):
    return {True:'.true.',False:'.false.'}[boolean]

class PhInput(dict):
    """
    Class to generate an manipulate the QuantumESPRESSO ph.x input files.
    Can be initialized either reading from a file or starting from scratch.


    Note that the default parameters for the init of the class set the ``qplot=True``
    and ``ldisp=False`` so ph.x performs the calculation on the list of q-points
    in the input in the form

        nqs \n
        xq(1,i)    xq(2,i)    xq(3,1)    nq(1) \n
        ...                                    \n
        xq(1,nqs)  xq(2,nqs)  xq(3,nqs)  nq(nqs) \n

    where nqs are the number of points, xq(j,i) is the j-th coordinate of the i-th point,
    in units of 2pi/a0 (a0 = lattice parameter), and nq(i) is the weigth of the i-th point

    """

    namelist = ['inputph']
    cards =['kpoints']

    def __init__(self,file=None,tr2_ph=1e-12,trans=True,qplot=True,ldisp=False,prefix='ph',outdir='./'):
        """
        Initialize the keys of the object and update the dictionaries with the kwargs passed as
        input parameters.. If an input file is provided it is parsed and the 'file' key is added
        to the object dictionary.

        Args:
            file (:py:class:`string`) : name of an exsistent input file, used to
                initialize the dictionaries of the object
            **kwargs : keyword arguments used to initialize the dictionaries of the
                object

        """
        dict.__init__(self)

        for key in self.namelist:
            self[key] = dict()
        for key in self.cards:
            self[key] = dict()
        default = {'inputph':dict(tr2_ph=tr2_ph,trans=trans,qplot=qplot,ldisp=ldisp,prefix=prefix,outdir=outdir)}
        self.update(default)

        if file is not None:
            self['file'] = file
            self.parseInputFile(file)

    def parseInputFile(self,file):
        """
        Read the arguments and variables from the input file.

        Args:
            file (:py:class:`string`) : name of an exsistent input file, used to
                initialize the dictionary of the object

        """
        f = open(file,"r")

        self.file_lines = f.readlines()
        for group in self.namelist:
            self.store(group)
        self.read_kpoints()

    def write(self,file):
        """
        Write the QE input on file.

        Args:
            file (:py:class:`string`) : name of the file
        """
        f = open(file,'w')
        f.write(self.convert_string())
        f.close()

    def convert_string(self):
        """
        Convert the input object into a string
        """
        lines = []

        for group in self.namelist:
            if self[group] != {}:
                lines.append(self.stringify_group(group))

        self.stringify_kpoints(lines)

        return '\n'.join(lines)

    def slicefile(self, group):
        """
        Return a list that contains the variables associated to the group
        key of the input file
        """
        import re
        lines = re.findall('&%s(?:.?)+\n((?:.+\n)+?)(?:\s+)?\/'%group,"".join(self.file_lines),re.MULTILINE)
        return lines

    def store(self,group):
        """
        Look for the namelist (control, system, electrons,...) in the file and
        attribute the associated variables in the dictionary
        """
        import re
        from mppi.Utilities import Utils
        for file_slice in self._slicefile(group):
            for key, value in re.findall('([a-zA-Z_0-9_\(\)]+)(?:\s+)?=(?:\s+)?([a-zA-Z/\'"0-9_.-]+)',file_slice):
                self[group][key.strip()]=Utils.convertTonumber(value.strip())

    def read_kpoints(self):
        """
        Read the kpoints from theinput file and attribute the associated variables
        in the dictionary....to be implemented
        """
        print('The parser of the phonon kpoints has not been implemented!')

    def stringify_group(self,group):
        variables = self[group]
        string='&%s\n' %group
        for var in sorted(variables):
            string += "%20s = %s\n" % (var, variables[var])
        string += "/&end"
        return string

    def stringify_kpoints(self,line):
        if 'values' in self['kpoints'] :
            line.append('%s'%len(self['kpoints']['values']))
            for kpt in self['kpoints']['values']:
                line.append( ('%12.8lf %12.8lf %12.8lf %.0f ')%tuple(kpt))

    def set_kpoints(self,klist):
        """
        Set the kpoints on which phonons are computed. Actually only the method associated
        to  ``qplot=True`` and ``ldisp=False`` is implemented.

        Args:
            klist (:py:class:`list`) : list with the coordinates and the weights of the kpoints
            kweigth (:py:class:`list`) : array with the weigth of each kpoint. If is None
                a uniform weight equal to 1 is attributed to each kpoint

        """
        self['kpoints'] = {'values' : klist}

    def set_prefix(self,prefix):
        """
        Set the value of prefix

        Args:
            prefix (:py:class:`string`) : value of the prefix

        """
        self['inputph']['prefix'] = "'%s'"%prefix

    def set_outdir(self,outdir):
        """
        Set the value of outdir

        Args:
            outdir (:py:class:`string`) : value of the outdir

        """
        self['inputph']['outdir'] = "'%s'"%outdir

    def set_inputph_variables(inp,**kwargs):
        """
        Add to the `inputph` key of the input dictionary
        the elements kwargs[key] = kwargs[value] for all the
        elements of the kwargs provided as input.

        Args:
            kwargs : variable(s) added in the form name = value

        """
        for name,value in kwargs.items():
            inp['inputph'][name] = [value,units]

    # Get methods

    def get_prefix(self):
        """
        Get the value of prefix.

        Returns:
            :py:class:`string` : The value of the prefix key of the input dictionary.
            If the key is not present return the default value 'pwscf'

        """
        pref = self['inputph'].get('prefix','pwscf')
        return pref.strip("'")


    def get_outdir(self):
        """
        Get the value of outdir.

        Returns:
            :py:class:`string` : The value of the outdir key of the input dictionary.
                If the key is not present return the default value '.'

        """
        outdir = self['inputph'].get('outdir','.')
        return outdir.strip("'")
