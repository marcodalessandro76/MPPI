
"""
This module manages the input file for pw.x computations of QuantumESPRESSO.
The input can be created from scratch or can be initialized from an existing input file.
"""
import os
from copy import deepcopy

def fortran_bool(boolean):
    return {True:'.true.',False:'.false.'}[boolean]

class PwInput(dict):
    """
    Class to generate an manipulate the QuantumESPRESSO pw.x input files.
    Can be initialized either reading from a file or starting from scratch.
    All the elements that belong to the possible namelist of pw are parsed.
    Instead, among the possible pw cards the actual implementation considers
    only : ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS.

    """

    namelist = ['control','system','electrons','ions','cell']
    cards =['atomic_species','atomic_positions','kpoints','cell_parameters']

    default = {'control':{
                   'calculation':"'scf'",
                   'verbosity':"'high'",
                   'prefix':"'pwscf'",
                   'outdir':"'./'"},
               'system':{
                    'force_symmorphic':fortran_bool(False)},
                'electrons':{
                    'diago_full_acc':fortran_bool(False)}
               }

    def __init__(self,file=None,**kwargs):
        """
        Initialize the keys of the object with the namelist and cards and update the
        dictionaries with the kwargs passed as input parameters. Some keys have a
        default value specified bt the class member dictionary` default`. If an input file
        is provided it is parsed and the 'file' key is added to the object dictionary.

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
        self.update(deepcopy(self.default))
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

        Args:
            file (:py:class:`string`) : name of an exsistent input file, used
                initialize the dictionaries of the object

        """
        f = open(file,"r")

        self.file_lines = f.readlines()
        for group in self.namelist:
            self._store(group)

        self._read_atomic_species()
        self._read_atomic_positions()
        self._read_cell_parameters()
        self._read_kpoints()

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
        from mppi.Utilities import Utils
        for file_slice in self._slicefile(group):
            for key, value in re.findall('([a-zA-Z_0-9_\(\)]+)(?:\s+)?=(?:\s+)?([a-zA-Z/\'"0-9_.-]+)',file_slice):
                self[group][key.strip()]=Utils.convertTonumber(value.strip())

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
                    self['atomic_species'][atype] = [float(mass),psp]

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
                        cell_parameters[i] = [ float(x) for x in next(lines).split() ]
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

    def set_energy_cutoff(self,energy):
        """
        Set the value of the kinetic energy cutoff (Ry) for wavefunctions

        Args:
            energy (:py:class:`float`) : value of the energy cutoff in Ry

        """
        self['system']['ecutwfc'] = energy

    def set_prefix(self,prefix):
        """
        Set the value of prefix

        Args:
            prefix (:py:class:`string`) : value of the prefix

        """
        self['control']['prefix'] = "'%s'"%prefix

    def set_outdir(self,outdir):
        """
        Set the value of outdir

        Args:
            outdir (:py:class:`string`) : value of the outdir

        """
        self['control']['outdir'] = "'%s'"%outdir

    def set_pseudo_dir(self, pseudo_dir='pseudos', abs_path=False):
        """
        Set the position of the folder with the pseudo-potentials.
        If `abs_path` is True the path is converted in a absolute path. In this way it is possible to provide
        a relative path (expressed from the root of the folder where the notebook is located) and the pseudo
        location can be found from an arbitrary folder.

        Args:
            pseudo_dir (:py:class:'string') : (relative) path of the folder with the pseduopotentials

        Note:
            If the folder tree contains blank spaces, QuantumESPRESSO cannot be able to find the pseudo, in this
            cas it is safer to provide a relative path (expressed from the folder where the input file is
            written)

        """
        if abs_path: pseudo_dir = os.path.abspath(pseudo_dir)
        self['control']['pseudo_dir'] = "'%s'"%pseudo_dir

    def set_occupations(self,occupations='fixed',smearing='fermi-dirac',degauss=50.):
        """
        Set the type of orbital of occupations of the ks states. If a smearing is
        present set also the type of smearing and the value of the degauss parameter.

        Args:
            occupations (:py:class:`string`) : type of occupation of the ks states (fixed, smearing,...)
            smearing (:py:class:`string`) : type of smearing (gaussian, fermi-dirac,...)
            degauss (:py:class:`float`) : value of the gaussian spreading (meV) for brillouin-zone
                integration in metals
        """
        from mppi.Utilities.Constants import HaToeV
        self['system']['occupations'] = "'"+occupations+"'"
        if occupations == 'smearing':
            self['system']['smearing'] = "'"+smearing+"'"
            self['system']['degauss'] = degauss/(0.5*1e3*HaToeV)

    def set_scf(self,conv_thr=1e-8,diago_full_acc=False,
                force_symmorphic=False,verbosity='high'):
        """
        Set the variables for a scf calculation.

        Args:
            conv_thr (:py:class:`string`) : the convergence threshold value
            diago_full_acc (:py:class:`bool`)
            force_symmorphic (:py:class:`bool`)
            verbosity (:py:class:`string`)

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
            nbnd (:py:class:`int`) : number of bands
            conv_thr (:py:class:`float`) : the convergence threshold value
            diago_full_acc (:py:class:`bool`)
            force_symmorphic (:py:class:`bool`)
            verbosity (:py:class:`string`)

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
            nbnd (:py:class:`int`) : number of bands
            conv_thr (:py:class:`float`) : the convergence threshold value
            diago_full_acc (:py:class:`bool`)
            force_symmorphic (:py:class:`bool`)
            verbosity (:py:class:`bool`)

        """
        self['control']['calculation'] = "'bands'"
        self['control']['verbosity'] = "'"+verbosity+"'"
        self['system']['nbnd'] = nbnd
        self['electrons']['conv_thr'] = conv_thr
        self['electrons']['diago_full_acc'] = fortran_bool(diago_full_acc)
        self['system']['force_symmorphic'] = fortran_bool(force_symmorphic)

    def set_num_bnds(self,nbnd):
        """
        Set the value of the variable nbnd. This method is useful if we need to
        set the number of bands in a scf calculation to include some empty bands.

        Args:
            nbnd (:py:class:`int`) : number of bands

        """
        self['system']['nbnd'] = nbnd

    def add_atom(self,atom,pseudo_name,mass = '1.0'):
        """
        Update the self['atomic_species'] dictionary

        Args:
            atom (:py:class:`string`)
            mass (:py:class:`string`) : is used only for molecular dynamics run
            pseudo_name (:py:class:`string`)
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
            nat (:py:class:`int`) : number of atoms in the cell
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
            type (:py:class:`string`) : units for the positions (alat,angstrom,crystal,...)
            positions (:py:class:`list`) : a list with the structure
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
            ibrav (:py:class:`int`) : Bravais-lattice index. 0 : Free, 1 : Cubic, 2 : Fcc, 3 : bcc,
                4 : Hexagonal
            cell_vectors (:py:class:`list`) : the list with the crystal lattice vectors in the
                form [[v1_x,v1_y,v1_z],[v2_x,v2_y,v2_z],[v3_x,v3_y,v3_z]]
            cell_units (:py:class:`string`) : units used for the cell vectors (angstrom, BOHR or alat).
                If alat is used the lattice vectors are in units of the lattice parameter celldm(1)
            celldms (:py:class:`float`) : Crystallographic constants. Only some values depending on
                `ibrav` have to be specified. For instance, alat = celldm(1) is the lattice parameter
                (in BOHR). If ibrav==0, only celldm(1) is used if present and cell vectors are
                read from card CELL_PARAMETERS

        """
        if ibrav == 0 and cell_vectors is None:
            raise ValueError('ibrav = 0 implies that the cell_parameters variable is set')
        if cell_vectors is not None:
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
            type (:py:class:`string`) : type of sampling (automatic, tpiba, tpiba_b,...)
            points (:py:class:`list`) : number of kpoints in the x,y,z directions. Used only if
                       the type variable is set to `automatic`
            shift (:py:class:`list`) : shifts in the x,y,z directions. Used only if the
                       type varible is set to `automatic`
            klist (:py:class:`list`) : list with the structure:
                       [[k1x,k1y,k1z,w1],[k2x,k2y,k2z,w2],....]
                       Used if type variable is not se to `automatic`

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

    # Get methods

    def get_prefix(self):
        """
        Get the value of prefix.

        Returns:
            :py:class:`string` : The value of the prefix key of the input dictionary.
            If the key is not present return the default value 'pwscf'

        """
        pref = self['control'].get('prefix','pwscf')
        return pref.strip("'")

    def get_outdir(self):
        """
        Get the value of outdir.

        Returns:
            :py:class:`string` : The value of the outdir key of the input dictionary.
            If the key is not present return the default value '.'

        """
        outdir = self['control'].get('outdir','.')
        return outdir.strip("'")
