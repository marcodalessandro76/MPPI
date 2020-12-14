"""
Class to perform the parsing of a QuantumESPRESSO XML file. Makes usage of the
data-file-schema.xml file that is found in the run_dir/prefix.save folder.
"""

from mppi.Utilities import HaToeV
from mppi.Parsers import ParsersUtils as U

class PwParser():
    """
    Initialize the data member of the class.

    Args:
        file (str): The name, including the path of the data-file-schema.xml
        verbose (bool) : set the amount of information written on terminal

    Attributes:
        natoms : number of atoms in the cell
        natypes : number of atomic species
        atomic_positions : list with the position of each atom
        atomic_species : dictionary with mass and pseudo for each species
        alat : lattice parameter (in a.u.)
        lattice : array with the lattice vectors. The i-th row represents the
            i-th lattice vectors in cartesian units
        syms : list with the symmetries of the lattice
        num_electrons : number of electrons
        nkpoints : numer of kpoints
        nbands : number of bands
        nbands_full : number of occupied bands (for systems with a gap)
        nbands_empty : number of empty bands (for systems with a gap)
        occupations_kind : type of occupation (fixed or smearing)
        kpoints : array with the kpoints
        occupations : array with the bands occupation for each kpoint
        weights : array with the weight of each kpoint, each element is a
            one-dimensional array.
        energy : total energy of the system (in Hartree)
        evals : array of the ks energies for each kpoint (in Hartree)
        lsda : True if collinear spin is activated
        noncolin : True if noncollinear spin calculation is activated
        spinorbit : True if spin-orbit couping is present
        spin_degen : 1 if lsda or non collinear spin is activated, 2 otherwise

    """

    def __init__(self,file,verbose=True):
        self.file = file
        if verbose: print('Parse file : %s'%self.file)
        try:
            self.parseXML(self.file)
        except TypeError: #FileNotFoundError or TypeError:
            if verbose: print('Failed to read %s'%self.file)
            self.data = None

    def parseXML(self,file):
        """
        Read the data from the xml file in the new format of QuantumESPRESSO.
        Some variable are extracted from the XML file and stored in the attributes
        of the object.
        """
        import xml.etree.ElementTree as ET
        import numpy as np

        self.data = ET.parse(file).getroot()

        # Add some relevant attributes to the object

        #units used in the XML data file
        self.units = self.data.attrib['Units']

        #atomic number and positions
        self.natoms = int(self.data.find("output/atomic_structure").get('nat'))
        self.atomic_positions = []
        atoms = self.data.findall("output/atomic_structure/atomic_positions/atom")
        for i in range(self.natoms):
            atype = atoms[i].get('name')
            pos = [float(x) for x in atoms[i].text.strip().split()]
            self.atomic_positions.append([atype,pos])

        #atomic species
        self.natypes = int(self.data.find("output/atomic_species").get('ntyp'))
        atypes = self.data.findall("output/atomic_species/species")
        self.atomic_species = {}
        for i in range(self.natypes):
            atype_string = atypes[i].get('name')
            atype_mass = atypes[i].findall('mass')[0].text
            atype_pseudo = atypes[i].findall('pseudo_file')[0].text
            self.atomic_species[atype_string]=[atype_mass,atype_pseudo]

        #lattice properties
        self.alat = float(self.data.find("output/atomic_structure").get('alat'))
        lattice = []
        for i in [1,2,3]:
            lat_vect = self.data.findall("output/atomic_structure/cell/a%d"%i)[0]
            lat_vect = [float(x) for x in lat_vect.text.strip().split()]
            lattice.append(lat_vect)
        self.lattice = np.array(lattice)
        symmetries = self.data.findall('output/symmetries/symmetry')
        self.syms = []
        for sym in symmetries:
            rotation = [float(x) for x in sym.find('rotation').text.split()]
            self.syms.append(np.array(rotation).reshape(3,3))

        #number of electrons
        self.num_electrons = float(self.data.find('output/band_structure/nelec').text)

        #number of kpoints and bands
        self.nkpoints = int(self.data.find('output/band_structure/nks').text)
        self.nbands = int(self.data.find('output/band_structure/nbnd').text)

        #spin related properties and spin-orbit coupling
        lsda = self.data.find('output/band_structure/lsda').text
        self.lsda = True if lsda == 'true' else False
        noncolin = self.data.find('output/band_structure/noncolin').text
        self.noncolin = True if noncolin == 'true' else False
        spinorbit = self.data.find('output/band_structure/spinorbit').text
        self.spinorbit = True if spinorbit == 'true' else False
        self.spin_degen = 1 if self.lsda or self.noncolin else 2

        #number of occupied and empty bands (for systems with a gap)
        self.nbands_full = int(self.num_electrons/self.spin_degen)
        self.nbands_empty = self.nbands-self.nbands_full

        #total energy
        self.energy = float(self.data.find('output/total_energy/etot').text)

        #occupations kind (fixed or smearing)
        self.occupations_kind = self.data.find('output/band_structure/occupations_kind').text

        kstates = self.data.findall('output/band_structure/ks_energies')

        #arrays with the kpoints, the associated weights, the ks energies and the occupations
        self.kpoints = []
        self.weights = []
        self.evals = []
        self.occupations = []
        for k in kstates:
            kpoint = [float(x) for x in k.find('k_point').text.split()]
            weight = [float(x) for x in k.find('k_point').get('weight').split()]
            eval = [float(x) for x in k.find('eigenvalues').text.split()]
            occ = [float(x) for x in k.find('occupations').text.split()]
            self.kpoints.append(kpoint)
            self.weights.append(weight)
            self.evals.append(eval)
            self.occupations.append(occ)
        self.kpoints = np.array(self.kpoints)
        self.weights = np.array(self.weights)
        self.evals = np.array(self.evals)
        self.occupations = np.array(self.occupations)

    def get_energy(self,convert_eV = True):
        """
        Return the total energy the system. If convert_eV is True the energy
        is provided in eV other the Hartree units are used.

        """
        if convert_eV:
            return HaToeV*self.energy
        else:
            return self.energy

    def get_fermi(self,convert_eV = True):
        """
        Return the fermi energy of the system (if present in the xml file).
        If convert_eV is True the fermi energy is provided in eV, otherwise the Hartree units are used.

        """
        fermi = self.data.find('output/band_structure/fermi_energy')
        if fermi is not None:
            fermi = float(fermi.text)
            if convert_eV: fermi *= HaToeV
            return fermi
        else:
            print('Fermi energy attribute not found in the ouput file. Maybe `fixed` occupation type is used?')
            return None

    def eval_lattice_volume(self):
        """
        Compute the volume of a lattice (in a.u.)

        Returns:
            :py:class:`float` : lattice volume in a.u.
        """
        return U.eval_lattice_volume(self.lattice)

    def get_lattice(self, rescale = False):
        """
        Compute the lattice vectors. If rescale = True the vectors are expressed in units
        of the lattice constant.

        Args:
            rescale (:py:class:`bool`)  : if True express the lattice vectors in units alat

        Returns:
            :py:class:`array` : array with the lattice vectors a_i as rows

        """
        return U.get_lattice(self.lattice,self.alat,rescale=rescale)

    def get_reciprocal_lattice(self, rescale = False):
        """
        Compute the reciprocal lattice vectors. If rescale = False the vectors are normalized
        so that np.dot(a_i,b_j) = 2*np.pi*delta_ij, where a_i is a basis vector of the direct
        lattice. If rescale = True the reciprocal lattice vectors are expressed in units of
        2*np.pi/alat.

        Args:
            rescale (:py:class:`bool`)  : if True express the reciprocal vectors in units of 2*np.pi/alat

        Returns:
            :py:class:`array` : array with the reciprocal lattice vectors b_i as rows

        """
        return U.get_reciprocal_lattice(self.lattice,self.alat,rescale=rescale)

    def get_evals(self, set_scissor = None, set_gap = None, set_direct_gap = None, verbose = True):
        """
        Return the ks energies for each kpoint (in eV). The top of the valence band is used as the
        reference energy value. It is possible to shift the energies of the empty bands by setting an arbitrary
        value for the gap (direct or indirect) or by adding an explicit scissor.
        Implemented only for semiconductors, the energy shift of empty bands does not update their occupation levels.

        Args:
            set_scissor (:py:class:`float`) : set the value of the scissor (in eV) that is added to the empty bands.
                If a scissor is provided the set_gap and set_direct_gap parameters are ignored
            set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system. If set_gap is provided
                the set_direct_gap parameter is ignored
            set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.

        Return:
            :py:class:`numpy.array`  : an array with the ks energies for each kpoint

        """
        evals = U.get_evals(self.evals,self.nbands,self.nbands_full,
                set_scissor=set_scissor,set_gap=set_gap,set_direct_gap=set_direct_gap,verbose=verbose)
        return evals

    def get_transitions(self, initial = 'full', final = 'empty',set_scissor = None, set_gap = None, set_direct_gap = None):
        """
        Compute the (vertical) transitions energies. For each kpoint compute the transition energies, i.e.
        the (positive) energy difference (in eV) between the final and the initial states.

        Args:
            initial (string or list) : specifies the bands from which electrons can be extracted. It can be set to `full` or
                `empty` to select the occupied or empty bands, respectively. Otherwise a list of bands can be
                provided
            final  (string or list) : specifies the final bands of the excited electrons. It can be set to `full` or
                `empty` to select the occupied or empty bands, respectively. Otherwise a list of bands can be
                provided
            set_scissor (:py:class:`float`) : set the value of the scissor (in eV) that is added to the empty bands.
                If a scissor is provided the set_gap and set_direct_gap parameters are ignored
            set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system. If set_gap is provided
                the set_direct_gap parameter is ignored
            set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.

        Return:
            :py:class:`numpy.array`  : an array with the transition energies for each kpoint

        """
        transitions = U.get_transitions(self.evals,self.nbands,self.nbands_full,initial=initial,final=final,
                      set_scissor=set_scissor,set_gap=set_gap,set_direct_gap=set_direct_gap)
        return transitions

    def get_gap(self, verbose = True):
        """
        Compute the energy gap of the system (in eV). The method check if the gap is direct or
        indirect. Implemented and tested only for semiconductors.

        Return:
            :py:class:`dict` : a dictionary with the values of direct and indirect gaps and the positions
            of the VMB and CBM

        """
        gap = U.get_gap(self.evals,self.nbands_full,verbose=verbose)
        return gap
