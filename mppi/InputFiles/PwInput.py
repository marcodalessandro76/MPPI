
"""
Class to create and manipulate an input file for pw.x computations of QuantumESPRESSO.
The input can be created from scratch or can be initialized from an existing input file.
"""

def fortran_bool(boolean):
    return {True:'.true.',False:'.false.'}[boolean]

class PwInput(dict):
    """
    Class to generate an manipulate the pw.x input files.
    Can be initialized either reading from a file or starting from scratch.
    All the elements that belong to the possible namelist of pw are parsed.
    Instead, among the possible pw cards the actual implementation considers
    only : ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS.
    """

    _namelist = ['control','system','electrons','ions','cell']
    _cards =['atomic_species','atomic_positions','kpoints','cell_parameters']

    def __init__(self,file=None,**kwargs):
        """
        Initialized the keys of self with the namelist and cards and update the
        dictionaries with the kwargs passed as input parameters. If the file is
        provided parse it and add the 'file' key.
        """
        dict.__init__(self)

        for key in self._namelist:
            self[key] = dict()
        for key in self._cards:
            self[key] = dict()
        self.update(kwargs)

        if file is not None:
            self['file'] = file
            self.parseInputFile(file)

    def parseInputFile(self,file):
        """
        Read the arguments and variables from the input file.
        All the elements that belong to the possible namelist of pw are parsed.
        Instead, among the possible pw cards the actual implementation considers
        only : ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS
        """
        f = open(file,"r")

        self.file_lines = f.readlines()
        for group in self._namelist:
            self._store(group)

        self._read_atomic_species()
        self._read_atomic_positions()
        self._read_cell_parameters()
        self._read_kpoints()

    def write(self,file):
        """
        Write the QE input on file.

        Args:
            file (string) : the file name
        """
        f = open(file,'w')
        f.write(self.convert_string())
        f.close()

    def convert_string(self):
        """
        Convert the input object into a string
        """
        lines = []

        for group in self._namelist:
            if self[group] != {}:
                lines.append(self._stringify_group(group))

        self._stringify_atomic_species(lines)
        self._stringify_atomic_positions(lines)
        self._stringify_kpoints(lines)
        self._stringify_cell_parameters(lines)

        return '\n'.join(lines)

    def _slicefile(self, group):
        """
        Return a list that contains the variables associated to the group
        key of the input file
        """
        import re
        lines = re.findall('&%s(?:.?)+\n((?:.+\n)+?)(?:\s+)?\/'%group,"".join(self.file_lines),re.MULTILINE)
        return lines

    def _store(self,group):
        """
        Look for the namelist (control, system, electrons,...) in the file and
        attribute the associated variables in the dictionary
        """
        import re
        for file_slice in self._slicefile(group):
            for key, value in re.findall('([a-zA-Z_0-9_\(\)]+)(?:\s+)?=(?:\s+)?([a-zA-Z/\'"0-9_.-]+)',file_slice):
                self[group][key.strip()]=value.strip()

    def _read_atomic_species(self):
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

    def _read_atomic_positions(self):
        """
        Look for the ATOMIC_POSITIONS field in the input file and attribute the
        associated variables in the dictionary
        """
        lines = iter(self.file_lines)
        for line in lines:
            if "ATOMIC_POSITIONS" in line:
                type = line.replace('{','').replace('}','').strip().split()[1]
                self['atomic_positions']['type'] = type
                self['atomic_positions']['values'] = []
                for i in range(int(self['system']['nat'])):
                    atype, x,y,z = next(lines).split()
                    self['atomic_positions']['values'].append([atype,[float(i) for i in (x,y,z)]])

    def _read_cell_parameters(self):
        """
        If ibrav == 0 look for the CELL_PARAMETERS field in the input file and
        attribute the associated variables in the dictionary.
        If the units are angstrom or bohr alat is the norm of the first lattice
        vector, so celldm(1) is removed. If the units are alat then the lattice
        vectors are in units of the lattice parameter celldm(1).

        Check the behavior of the code when 'celldm(1)' is provided.
        It seems that the most safe approach is to do not set 'celldm(1)' and
        express the units in angstrom or bohr.
        """
        def rmchar(string,symbols): return ''.join([c for c in string if c not in symbols])
        ibrav = int(self['system']['ibrav'])
        if ibrav == 0:
            if 'celldm(1)' in list(self['system'].keys()):
                a = float(self['system']['celldm(1)'])
            else:
                a = 1
            cell = self['cell_parameters']
            lines = iter(self.file_lines)
            for line in lines:
                if "CELL_PARAMETERS" in line:
                    units = rmchar(line.strip(),'{}()').split()
                    cell_parameters = [[],[],[]]
                    if len(units) > 1:
                        cell['type'] = units[1]
                    else:
                        cell['type'] = 'bohr'
                    for i in range(3):
                        cell_parameters[i] = [ float(x)*a for x in next(lines).split() ]
            if cell['type'] == 'angstrom' or cell['type'] == 'bohr':
                if 'celldm(1)' in self['system']: del self['system']['celldm(1)']
            cell['values'] = cell_parameters

    def _read_kpoints(self):
        """
        Look for the K_POINTS field in the input file and attribute the
        associated variables in the dictionary
        """
        lines = iter(self.file_lines)
        for line in lines:
            if "K_POINTS" in line:
                kp = self['kpoints']
                if "automatic" in line:
                    kp['type'] = 'automatic'
                    vals = list(map(float, next(lines).split()))
                    kp['values'] = (vals[0:3],vals[3:6])
                else:
                    nkpoints = int(lines.__next__().split()[0])
                    kp['type'] = line.split()[2]
                    kp['values'] = []
                    try:
                        lines_list = list(lines)
                        for n in range(nkpoints):
                            vals = lines_list[n].split()[:4]
                            kp['values'].append( list(map(float,vals)) )
                    except IndexError:
                        print('wrong k-points list format')
                        exit()

    # Methods that convert the object attributes into a string with the
    # correct QuantumESPRESSO format.

    def _stringify_group(self,group):
        variables = self[group]
        string='&%s\n' %group
        for var in sorted(variables):
            string += "%20s = %s\n" % (var, variables[var])
        string += "/&end"
        return string

    def _stringify_atomic_species(self,line):
        if self['atomic_species'] != {}:
            line.append( 'ATOMIC_SPECIES' )
            for atom in self['atomic_species']:
                line.append(" %3s %8s %20s" % (atom, self['atomic_species'][atom][0],
                self['atomic_species'][atom][1]))

    def _stringify_atomic_positions(self,line):
        if self['atomic_positions'] != {}:
            line.append('ATOMIC_POSITIONS { %s }'%self['atomic_positions']['type'])
            for atom in self['atomic_positions']['values']:
                line.append('%3s %14.10lf %14.10lf %14.10lf'%
                (atom[0], atom[1][0], atom[1][1], atom[1][2]))

    def _stringify_kpoints(self,line):
        if self['kpoints'] != {}:
            line.append("K_POINTS { %s }" % self['kpoints']['type'])
            if self['kpoints']['type'] == 'automatic':
                line.append(("%3d"*6)%(tuple(self['kpoints']['values'][0]) +
                    tuple(self['kpoints']['values'][1])))
            else:
                line.append( "%d" % len(self['kpoints']['values']))
                for kpt in self['kpoints']['values']:
                    line.append( ("%12.8lf "*4)%tuple(kpt))

    def _stringify_cell_parameters(self,line):
        if self['cell_parameters'] != {}:
            if int(self['system']['ibrav']) == 0:
                line.append("CELL_PARAMETERS %s"%self['cell_parameters']['type'])
                line.append(("%14.10lf "*3)%tuple(self['cell_parameters']['values'][0]))
                line.append(("%14.10lf "*3)%tuple(self['cell_parameters']['values'][1]))
                line.append(("%14.10lf "*3)%tuple(self['cell_parameters']['values'][2]))

    # Set methods

    def set_energy_cutoff(self,e):
        """
        Set the value of the kinetic energy cutoff (Ry) for wavefunctions
        """
        self['system']['ecutwfc'] = e

    def set_prefix(self,prefix):
        """
        Set the value of prefix
        """
        self['control']['prefix'] = "'%s'"%prefix

    def set_pseudo_dir(self,pseudo_dir='../pseudos'):
        """
        Set the value of the pseudo_folder. The default is '../pseudos'
        because usually the pw input is written in the run_dir folder by
        the QeCalculator or Dataset.
        """
        self['control']['pseudo_dir'] = "'%s'"%pseudo_dir

    def set_occupations(self,occupations='fixed',smearing='fermi-dirac',degauss=50.):
        """
        Set the type of orbital of occupations of the ks states. If a smearing is
        present set also the type of smearing and the value of the degauss parameter.

        Args:
            occupations (string) : type of occupation of the ks states (fixed, smearing,...)
            smearing (string) : type of smearing (gaussian, fermi-dirac,...)
            degauss (float) : value of the gaussian spreading (meV) for brillouin-zone
                integration in metals
        """
        from mppi.Utilities import HaToeV
        self['system']['occupations'] = "'"+occupations+"'"
        if occupations == 'smearing':
            self['system']['smearing'] = "'"+smearing+"'"
            self['system']['degauss'] = degauss/(0.5*1e3*HaToeV)

    def set_scf(self,conv_thr=1e-8,diago_full_acc=False,
                force_symmorphic=False,verbosity='high'):
        """
        Set the variables for a scf calculation

        Args:
            conv_thr(float): the convergence threshold value
            diago_full_acc(boolean)
            force_symmorphic(boolean)
            verbosity(string)
        """
        self['control']['calculation'] = "'scf'"
        self['control']['verbosity'] = "'"+verbosity+"'"
        self['electrons']['conv_thr'] = conv_thr
        self['electrons']['diago_full_acc'] = fortran_bool(diago_full_acc)
        self['system']['force_symmorphic'] = fortran_bool(force_symmorphic)

    def set_nscf(self,nbnd,conv_thr=1e-8,diago_full_acc=False,
                 force_symmorphic=False,verbosity='high'):
        """
        Set the variables for a nscf calculation

        Args:
            nbnd(int): number of bands
            conv_thr(float): the convergence threshold value
            diago_full_acc(boolean)
            force_symmorphic(boolean)
            verbosity(string)
        """
        self['control']['calculation'] = "'nscf'"
        self['control']['verbosity'] = "'"+verbosity+"'"
        self['system']['nbnd'] = nbnd
        self['electrons']['conv_thr'] = conv_thr
        self['electrons']['diago_full_acc'] = fortran_bool(diago_full_acc)
        self['system']['force_symmorphic'] = fortran_bool(force_symmorphic)

    def set_bands(self,nbnd,conv_thr=1e-8,diago_full_acc=False,
                 force_symmorphic=False,verbosity='high'):
        """
        Set the variables for a bands calculation

        Args:
            nbnd(int): number of bands
            conv_thr(float): the convergence threshold value
            diago_full_acc(boolean)
            force_symmorphic(boolean)
            verbosity(string)
        """
        self['control']['calculation'] = "'bands'"
        self['control']['verbosity'] = "'"+verbosity+"'"
        self['system']['nbnd'] = nbnd
        self['electrons']['conv_thr'] = conv_thr
        self['electrons']['diago_full_acc'] = fortran_bool(diago_full_acc)
        self['system']['force_symmorphic'] = fortran_bool(force_symmorphic)

    def add_atom(self,atom,pseudo_name,mass = '1.0'):
        """
        Update the self['atomic_species'] dictionary

        Args:
            atom(str)
            mass(str): is used only for molecular dynamics run
            pseudo_name(str)
        """
        at = {atom : [mass,pseudo_name]}
        self['atomic_species'].update(at)

    def set_atoms_number(self,nat):
        """
        Set the self['control']['ntyp'] using the self['atomic_species']
        dictionary and set the number of atoms in the cell.
        If the number of atoms in the cell is lower than the number of atomic
        species write an alert and set nat = ntyp.

        Args:
            nat(int): number of atoms in the cell
        """
        ntyp = len(set(self['atomic_species'].keys()))
        if nat < ntyp:
            print('Number of atoms in the cell cannot be lower than number of atomic species')
            nat = ntyp
        self['system']['ntyp'] = str(ntyp)
        self['system']['nat'] = str(nat)

    def set_atomic_positions(self,positions,type='alat'):
        """
        Define the the positions of the atoms in the lattice.

        Args:
            type(str) : units for the positions (alat,angstrom,crystal,...)
            positions(list) : a list with the structure
                    [['atom1',[x,y,z]],['atom2',[x,y,z]],...]}
        """
        pos = {'type' : type, 'values' : positions}
        self['atomic_positions'] = pos

    def set_lattice(self,ibrav,cell_vectors=None,cell_units='angstrom',
                    celldm1=None,celldm2=None,celldm3=None,
                    celldm4=None,celldm5=None,celldm6=None):

        """
        Set the lattice structure using the typical QuantumESPRESSO input variables.

        Args:
            ibrav (int) : Bravais-lattice index. 0 : Free, 1 : Cubic, 2 : Fcc, 3 : bcc,
                4 : Hexagonal
            cell_vectors (:py:class:`list`) : the list with the crystal lattice vectors in the
                form [[v1_x,v1_y,v1_z],[v2_x,v2_y,v2_z],[v3_x,v3_y,v3_z]]
            cell_units (string) : units used for the cell vectors (angstrom, BOHR or alat).
                If alat is used the lattice vectors are in units of the lattice parameter celldm(1)
            celldms (**kwargs) : Crystallographic constants. Only needed values
                (depending on "ibrav") must be specified. alat = celldm(1) is the lattice parameter (in BOHR).
                If ibrav==0, only celldm(1) is used if present and cell vectors are
                read from card CELL_PARAMETERS

        """
        if ibrav == 0 and cell_vectors is None:
            raise ValueError('ibrav = 0 implies that the cell_parameters variable is set')
        if cell_vectors:
            self['cell_parameters'] = {'type' : cell_units, 'values' : cell_vectors}
        self['system']['ibrav'] = ibrav
        if celldm1 is not None: self['system']['celldm(1)'] = celldm1
        if celldm2 is not None: self['system']['celldm(2)'] = celldm2
        if celldm3 is not None: self['system']['celldm(3)'] = celldm3
        if celldm4 is not None: self['system']['celldm(4)'] = celldm4
        if celldm5 is not None: self['system']['celldm(5)'] = celldm5
        if celldm6 is not None: self['system']['celldm(6)'] = celldm6

    def set_kpoints(self,type='automatic',points=[1.,1.,1.],shift=[0.,0.,0.],klist=[]):
        """
        Define the sampling of the Brillouin zone.

        Args:
            type(str) : type of sampling (automatic, tpiba, tpiba_b,...)
            points(list) : number of kpoints in the x,y,z directions. Used only is
                       type is automatic
            shift(list) : shifts in the x,y,z directions. Used only is
                       type is automatic
            klist(list) : list with the structure:
                       [[k1x,k1y,k1z,w1],[k2x,k2y,k2z,w2],....]
                       Used only if type is not automatic
        """
        if type == 'automatic':
            k = {'type' : type, 'values' : (points,shift)}
        else:
            k = {'type' : type, 'values' : klist}
        self['kpoints'] = k

    def set_spinorbit(self):
        """
        Set the lspinorb and noncolin to True.
        """
        self['system']['lspinorb'] = '.true.'
        self['system']['noncolin'] = '.true.'

    def set_spinpolarized(self):
        """
        Set the value of nspin to 2. Useful to treat collinear spin polarized
        systems.
        """
        self['system']['lspinorb'] = '.false.'
        self['system']['noncolin'] = '.false.'
        self['system']['nspin']    = 2
