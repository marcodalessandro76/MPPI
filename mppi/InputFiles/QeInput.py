
"""
Module to create and manipulate an input file for QuantumESPRESSO computations.
The input can be created from scratch or can be initialized from an existing input file.
"""

#import os
import re

class QeInput(dict):
    """
    Class to generate an manipulate Quantum Espresso input files.
    Can be initialized either reading from a file or starting from scratch.
    """

    _basedict = {'control' : {'verbosity' : "'high'"},
                 'system' : {'force_symmorphic' : '.true', 'occupations' : "'fixed'"}} #,
                 #'atomic_pos_type' : 'crystal',
                 #'cell_units' : 'bohr'}

    # Define the group variables that are possible keys of self
    _groups = ['control','system','electrons','ions','cell']

    def __init__(self,filename=None):
        """
        Initalize the class
        """
        dict.__init__(self,filename=filename)

        # init the basic keys of self. The keys associated to the 'ions' and
        # 'cell' groups are initialized only if needed
        self['control'] = dict()
        self['system'] = dict()
        self['electrons'] = dict()
        # Variables associated to (some) QuantumESPRESSO input CARDS
        self['atomic_species'] = dict()
        self['atomic_positions'] = []

        # These keys has to be defined before the object is written to file
        #self['atomic_pos_type'] = ''
        #self['ktype'] = ''
        #self['cell_units'] = ''

        #self['cell_parameters'] = [] #init da read_cell_parameters

        #self['klist'] = []

        # set some keys to the values of _basedict
        #self.update(self._basedict)

        #in case we start from a reference file
        if filename:
            self.parseInputFile(filename)

    def parseInputFile(self,filename):
        """
        Read the arguments and variables from the input file
        """
        f = open(filename,"r")

        self.file_lines = f.readlines()
        self.store('control')
        self.store('system')
        self.store('electrons')
        self.store('ions')
        self.store('cell')

        self.read_atomicspecies()
        self.read_atoms()
        self.read_cell_parameters()
        self.read_kpoints()

    def write(self,filename):
        """
        Write the QE input on file.
        """
        f = open(filename,'w')
        f.write(self.convert_string())
        f.close()

    def convert_string(self):
        """
        Convert the input object into a string
        """
        lines = []

        for group in self._groups:
            if group in self.keys():
                lines.append(self.stringify_group(group))

        self.stringify_atomic_species(lines)
        self.stringify_atomic_positions(lines)
        self.stringify_kpoints(lines)
        self.stringify_cell_parameters(lines)

        return '\n'.join(lines)

    def slicefile(self, keyword):
        lines = re.findall('&%s(?:.?)+\n((?:.+\n)+?)(?:\s+)?\/'%keyword,"".join(self.file_lines),re.MULTILINE)
        return lines

    def store(self,group):
        """
        Look for the group (control, system, electrons,...) in the file and
        attribute the associated variables in the dictionary
        """
        for file_slice in self.slicefile(group):
            for key, value in re.findall('([a-zA-Z_0-9_\(\)]+)(?:\s+)?=(?:\s+)?([a-zA-Z/\'"0-9_.-]+)',file_slice):
                self[group][key.strip()]=value.strip()

    def read_atomicspecies(self):
        """
        Look for the ATOMIC_SPECIES field in the input file and attribute the
        associated variables in the dictionary
        """
        lines = iter(self.file_lines)
        for line in lines:
            if "ATOMIC_SPECIES" in line:
                for i in range(int(self['system']['ntyp'])):
                    atype, mass, psp = next(lines).split()
                    self['atomic_species'][atype] = [mass,psp]

    def read_atoms(self):
        """
        Look for the ATOMIC_POSITIONS field in the input file and attribute the
        associated variables in the dictionary
        """
        lines = iter(self.file_lines)
        for line in lines:
            if "ATOMIC_POSITIONS" in line:
                atomic_pos_type = line
                self['atomic_pos_type'] = re.findall('([A-Za-z]+)',line)[-1]
                for i in range(int(self['system']['nat'])):
                    atype, x,y,z = next(lines).split()
                    self['atomic_positions'].append([atype,[float(i) for i in (x,y,z)]])
        self['atomic_pos_type'] = atomic_pos_type.replace('{','').replace('}','').strip().split()[1]

    def read_cell_parameters(self):
        """
        If ibrav == 0 look for the CELL_PARAMETERS field in the input file and
        attribute the associated variables in the dictionary
        """
        import numpy as np
        def rmchar(string,symbols): return ''.join([c for c in string if c not in symbols])
        ibrav = int(self['system']['ibrav'])
        if ibrav == 0:
            if 'celldm(1)' in list(self['system'].keys()):
                a = float(self['system']['celldm(1)'])
            else:
                a = 1
            lines = iter(self.file_lines)
            for line in lines:
                if "CELL_PARAMETERS" in line:
                    units = rmchar(line.strip(),'{}()').split()
                    cell_parameters = [[],[],[]]
                    if len(units) > 1:
                        self['cell_units'] = units[1]
                    else:
                        self['cell_units'] = 'bohr'
                    for i in range(3):
                        cell_parameters[i] = [ float(x)*a for x in next(lines).split() ]
            if self['cell_units'] == 'angstrom' or self['cell_units'] == 'bohr':
                if 'celldm(1)' in self['system']: del self['system']['celldm(1)']
            if 'celldm(1)' not in list(self['system'].keys()):
                a = np.linalg.norm(cell_parameters[0])
            self['cell_parameters'] = cell_parameters

    def read_kpoints(self):
        """
        Look for the K_POINTS field in the input file and attribute the
        associated variables in the dictionary
        """
        lines = iter(self.file_lines)
        for line in lines:
            if "K_POINTS" in line:
                if "automatic" in line:
                    self['ktype'] = 'automatic'
                    vals = list(map(float, next(lines).split()))
                    self['kpoints'], self['shiftk'] = vals[0:3], vals[3:6]
                else:
                    nkpoints = int(lines.__next__().split()[0])
                    self['ktype'] = line.split()[2]
                    self['klist'] = []
                    try:
                        lines_list = list(lines)
                        for n in range(nkpoints):
                            vals = lines_list[n].split()[:4]
                            self['klist'].append( list(map(float,vals)) )
                    except IndexError:
                        print('wrong k-points list format')
                        exit()

    # Define the methods that convert the object attributes into a string with the
    # correct QuantumESPRESSO format.
    def stringify_group(self,group):
        """
        Convert the elements of a group (control, system, electrons,...) in a string
        with the proper QE input file format
        """
        variables = self[group]
        if variables != {}:
            string='&%s\n' %group
            for var in sorted(variables):
                string += "%20s = %s\n" % (var, variables[var])
            string += "/&end"
            return string
        else:
            return ''

    def stringify_atomic_species(self,line):
        line.append( 'ATOMIC_SPECIES' )
        for atype in self['atomic_species']:
            line.append(" %3s %8s %20s" % (atype, self['atomic_species'][atype][0],
            self['atomic_species'][atype][1]))

    def stringify_atomic_positions(self,line):
        line.append('ATOMIC_POSITIONS { %s }'%self['atomic_pos_type'])
        for atom in self['atomic_positions']:
            line.append('%3s %14.10lf %14.10lf %14.10lf'%
            (atom[0], atom[1][0], atom[1][1], atom[1][2]))

    def stringify_kpoints(self,line):
        line.append("K_POINTS { %s }" % self['ktype'])
        if self['ktype'] == 'automatic':
            line.append(("%3d"*6)%(tuple(self['kpoints']) + tuple(self['shiftk'])))
        else:
            line.append( "%d" % len(self['klist']))
            for kpt in self['klist']:
                line.append( ("%12.8lf "*4)%tuple(kpt))

    def stringify_cell_parameters(self,line):
        if int(self['system']['ibrav']) == 0:
            line.append("CELL_PARAMETERS %s"%self['cell_units'])
            line.append(("%14.10lf "*3)%tuple(self['cell_parameters'][0]))
            line.append(("%14.10lf "*3)%tuple(self['cell_parameters'][1]))
            line.append(("%14.10lf "*3)%tuple(self['cell_parameters'][2]))
