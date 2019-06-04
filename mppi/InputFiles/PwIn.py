
"""
PwIn class of YamboPy, Modified for compliance with python 3.
This module allows us to create and manipulate an input file for pw.x.
The input object can be created from scratch or can be init from an existing
file
"""

import os
import re
from math import sqrt

class PwIn():
    """
    Class to generate an manipulate Quantum Espresso input files
    Can be initialized either reading from a file or starting from a new file.
    """

    def __init__(self, filename=None):

        # init dictionaries
        self.control = dict()
        self.system = dict()
        self.electrons = dict()
        self.ions = dict()
        self.cell = dict()
        self.atypes = dict()
        self.atoms = []
        self.cell_parameters = []
        self.klist = []

        # set some variables to 'standard' values
        self.cell_units = 'angstrom'
        self.control['verbosity'] = "'high'"
        self.ktype = 'automatic'
        self.kpoints = [1,1,1]
        self.shiftk = [0,0,0]

        #in case we start from a reference file
        if filename:
            f = open(filename,"r")
            self.file_name = filename #set filename
            self.file_lines = f.readlines() #set file lines
            self.store(self.control,"control")     #read &control
            self.store(self.system,"system")      #read &system
            self.store(self.electrons,"electrons")   #read &electrons
            self.store(self.ions,"ions")        #read &ions
            self.store(self.cell,"cell")        #read &ions
            self.read_atomicspecies()  #read ATOMIC_SPECIES
            self.read_atoms()  #read ATOMIC_POSITIONS
            self.read_kpoints()   #read K_POINTS
            self.read_cell_parameters()   #read CELL_PARAMETERS

    # get methods
    def get_symmetry_spglib(self):
        """
        get the symmetry group of this system using spglib
        """
        import spglib

        lat, positions, atypes = self.get_atoms()
        lat = np.array(lat)

        at = np.unique(atypes)
        an = dict(zip(at,xrange(len(at))))
        atypes = [an[a] for a in atypes]

        cell = (lat,positions,atypes)

        spacegroup = spglib.get_spacegroup(cell,symprec=1e-5)
        return spacegroup

    def get_masses(self):
        """ Get an array with the masses of all the atoms
        """
        masses = []
        for atom in self.atoms:
            atype = self.atypes[atom[0]]
            mass = float(atype[0])
            masses.append(mass)
        return masses

    def get_atoms(self):
        """ Get the lattice parameters, postions of the atoms and chemical symbols
        """
        self.read_cell_parameters()
        cell = self.cell_parameters
        sym = [atom[0] for atom in self.atoms]
        pos = [atom[1] for atom in self.atoms]
        if self.atomic_pos_type == 'bohr':
            pos = car_red(pos,cell)
        return cell, pos, sym

    # set methods
    def set_kpoints(self,**kwargs):
        """
        Set the numbers of K_POINTS
        If ktpye = automatic needs kpoints = [kx,ky,kz]
        if ktype = tpiba needs klist = [P1+[weight1],P2+[weight2],...]
        where P-i is the i-th point of the grid and weight-i is the associated weight.
        In a non-scf calculation, weights do not affect the results.
        If ktype = tpiba_b needs klist = [P1+[numStep1],P2+[numStep2],...]
        where P-i is the i-th point of the path and numStep-i is the number of step
        from P-i to P-i+1
        """
        ktype = kwargs['ktype']
        self.ktype = ktype
        if ktype == 'automatic':
            self.kpoints = kwargs['kpoints']
            if 'shiftk' in kwargs :
                self.shiftk = kwargs['shiftk']
        if ktype in ['tpiba','tpiba_b']:
            self.klist = kwargs['klist']

    def set_atoms_type(self,n):
        """
        Set the number of atomic species
        """
        self.system['ntyp'] = n

    def set_atoms_number(self,n):
        """
        Set the number of atoms in the cell
        """
        self.system['nat'] = n

    def set_atoms_position(self,pos_type='angstrom',pos_list=[]):
       """
       Set the position of atoms in the cell
       Args:
            pos_type (str) : unit in which the position are given (angstrom, crystal,...)
            pos_list (list) : given in the form [['atoms1',[x1,y1,z1]],['atoms2'],[x2,y2,z2],...]
       """
       self.atomic_pos_type = pos_type
       if 'nat' in self.system:
           if len(pos_list) == self.system['nat']:
               self.atoms=pos_list
           else : print('position list does not match with the number of atoms')
       else : print('provide the number of atoms in the cell')

    def set_energy_cutoff(self,e):
        """
        Set the energy cutoff value (in Ry)
        """
        self.system['ecutwfc'] = e

    def set_num_bands(self,nb):
        """
        Set the number of bands
        """
        self.system['nbnd'] = nb

    def set_convergence_thr(self,thr):
        """
        Set the threshold level to seek convergence
        """
        self.electrons['conv_thr'] = thr

    def set_prefix(self,prefix):
        """
        Set the prefix pf the .save folder.
        The path starts from the position of the input file.
        """
        self.control['prefix'] = "'"+prefix+"'"

    def set_calculation(self,calc):
        """
        Set the calculation type (scf,nscf,...).
        """
        self.control['calculation'] = "'"+calc+"'"

    def set_occupations(self,occ):
        """
        Set the electron occupatopns (fixed,smearing,...).
        """
        self.system['occupations'] = "'"+occ+"'"

    def set_pseudo_dir(self,pseudo_dir):
        """
        Set the location of the folder with the pseudopotential.
        The path starts from the position of the input file.
        """
        self.control['pseudo_dir'] = "'"+pseudo_dir+"'"

    def set_path(self,path):
        self.klist = path.get_klist()

    def read_atomicspecies(self):
        lines = iter(self.file_lines)
        #find ATOMIC_SPECIES keyword in file and read next line
        for line in lines:
            if "ATOMIC_SPECIES" in line:
                for i in range(int(self.system["ntyp"])):
                    atype, mass, psp = next(lines).split()
                    self.atypes[atype] = [mass,psp]

    def read_atoms(self):
        lines = iter(self.file_lines)
        #find READ_ATOMS keyword in file and read next lines
        for line in lines:
            if "ATOMIC_POSITIONS" in line:
                atomic_pos_type = line
                self.atomic_pos_type = re.findall('([A-Za-z]+)',line)[-1]
                for i in range(int(self.system["nat"])):
                    atype, x,y,z = next(lines).split()
                    self.atoms.append([atype,[float(i) for i in (x,y,z)]])
        self.atomic_pos_type = atomic_pos_type.replace('{','').replace('}','').strip().split()[1]

    def read_cell_parameters(self):
        ibrav = int(self.system['ibrav'])
        if ibrav == 0:
            if 'celldm(1)' in self.system.keys():
                a = float(self.system['celldm(1)'])
            else:
                a = 1
            lines = iter(self.file_lines)
            for line in lines:
                if "CELL_PARAMETERS" in line:
                    #self.cell_units = line.translate(None, '{}()').split()[1]
                    self.cell_units = line.translate(str.maketrans('','','{}()')).split()[1]
                    self.cell_parameters = [[1,0,0],[0,1,0],[0,0,1]]
                    for i in range(3):
                        self.cell_parameters[i] = [ float(x)*a for x in next(lines).split() ]
            if self.cell_units == 'angstrom' or self.cell_units == 'bohr':
                if 'celldm(1)' in self.system: del self.system['celldm(1)']
        elif ibrav == 4:
            a = float(self.system['celldm(1)'])
            c = float(self.system['celldm(3)'])
            self.cell_parameters = [[   a,          0,  0],
                                    [-a/2,sqrt(3)/2*a,  0],
                                    [   0,          0,c*a]]
        elif ibrav == 2:
            a = float(self.system['celldm(1)'])
            self.cell_parameters = [[ -a/2,   0, a/2],
                                    [    0, a/2, a/2],
                                    [ -a/2, a/2,   0]]
        elif ibrav == 6:
            a = float(self.system['celldm(1)'])
            c = float(self.system['celldm(3)'])
            self.cell_parameters = [[  a,   0,   0],
                                    [  0,   a,   0],
                                    [  0,   0, c*a]]
        else:
            print('ibrav = %d not implemented'%ibrav)
            exit(1)

    def read_kpoints(self):
        lines = iter(self.file_lines)
        #find K_POINTS keyword in file and read next line
        for line in lines:
            if "K_POINTS" in line:
                #check if the type is automatic
                if "automatic" in line:
                    self.ktype = "automatic"
                    vals = list(map(float, next(lines).split()))
                    self.kpoints, self.shiftk = vals[0:3], vals[3:6]
                #otherwise read a list
                else:
                    #read number of kpoints
                    nkpoints = int(lines.next().split()[0])
                    self.klist = []
                    self.ktype = ""
                    try:
                        lines_list = list(lines)
                        for n in range(nkpoints):
                            vals = lines_list[n].split()[:4]
                            self.klist.append( map(float,vals) )
                    except IndexError:
                        print('wrong k-points list format')
                        exit()

    def slicefile(self, keyword):
        lines = re.findall('&%s(?:.?)+\n((?:.+\n)+?)(?:\s+)?\/'%keyword,"".join(self.file_lines),re.MULTILINE)
        return lines

    def store(self,group,name):
        """
        Save the variables specified in each of the groups on the structure
        """
        for file_slice in self.slicefile(name):
            for keyword, value in re.findall('([a-zA-Z_0-9_\(\)]+)(?:\s+)?=(?:\s+)?([a-zA-Z/\'"0-9_.-]+)',file_slice):
                group[keyword.strip()]=value.strip()

    def stringify_group(self, keyword, group):
        if group != {}:
            string='&%s\n' % keyword
            for keyword in group:
                string += "%20s = %s\n" % (keyword, group[keyword])
            string += "/&end\n"
            return string
        else:
            return ''

    def remove_key(self,group,key):
        """ if a certain key exists in the group, remove it
        """
        if key in group.items():
            del group[key]

    def write(self,filename):
        f = open(filename,'w')
        f.write(str(self))
        f.close()

    def __str__(self):
        """
        Output the file in the form of a string
        """
        string = ''
        string += self.stringify_group("control",self.control) #print control
        string += self.stringify_group("system",self.system) #print system
        string += self.stringify_group("electrons",self.electrons) #print electrons
        string += self.stringify_group("ions",self.ions) #print ions
        string += self.stringify_group("cell",self.cell) #print ions

        #print atomic species
        string += "ATOMIC_SPECIES\n"
        for atype in self.atypes:
            string += " %3s %8s %20s\n" % (atype, self.atypes[atype][0], self.atypes[atype][1])
        #print atomic positions
        string += "ATOMIC_POSITIONS { %s }\n"%self.atomic_pos_type
        for atom in self.atoms:
            string += "%3s %14.10lf %14.10lf %14.10lf\n" % (atom[0], atom[1][0], atom[1][1], atom[1][2])
        #print kpoints
        if self.ktype == "automatic":
            string += "K_POINTS { %s }\n" % self.ktype
            string += ("%3d"*6+"\n")%tuple(self.kpoints + self.shiftk)
        elif self.ktype == "crystal":
            string += "K_POINTS { %s }\n" % self.ktype
            string += "%d\n" % len(self.klist)
            for i in self.klist:
              string += ('%12.8lf '*4+'\n') % tuple(i)
        else:
            string += "K_POINTS { %s }\n" % self.ktype
            string += "%d\n" % len(self.klist)
            for i in self.klist:
                string += (("%12.8lf "*4)+"\n")%tuple(i)
        if self.system['ibrav'] == 0 or self.system['ibrav'] == '0':
            string += "CELL_PARAMETERS %s\n"%self.cell_units
            for i in range(3):
                string += ("%14.10lf "*3+"\n")%tuple(self.cell_parameters[i])
        return string
