"""
Class to perform the parsing of a QuantumESPRESSO XML file. Makes usage of the
data-file-schema.xml file that is found in the run_dir/prefix.save folder.
"""

from mppi.Utilities import HaToeV

class PwParser():
    """
    Initialize the data member of the class. The init method is able to manage a
    TypeError since, if the QeCalculator simulation fails the name of the output
    file is set to None.

    Args:
        file (str): The name, including the path of the data-file-schema.xml
        verbose (bool) : set the amount of information written on terminal

    Attributes:
        natoms : number of atoms in the cell
        natypes : number of atomic species
        atomic_positions : list with the position of each atom
        atomic_species : dictionary with mass and pseudo for each species
        nkpoints : numer of kpoints
        nbands : number of bands
        occupations : list with the bands occupation for each kpoint
        weights : list with thw weight of each kpoint
        energy : total energy of the system (in Hartree)
        fermi : fermi energy (in Hartree)
        evals : list of the ks energies for each kpoint (in Hartree)
        lsda : True if collinear spin is activated
        noncolin : True if noncollinear spin calculation is activated
        spinorbit : True if spin-orbit couping is present
        spin_degen : 1 if lsda or non collinear spin is activated, 0 otherwise

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
        self.natoms = int(self.data.findall("output/atomic_structure")[0].get('nat'))
        self.atomic_positions = []
        atoms = self.data.findall("output/atomic_structure/atomic_positions/atom")
        for i in range(self.natoms):
            atype = atoms[i].get('name')
            pos = [float(x) for x in atoms[i].text.strip().split()]
            self.atomic_positions.append([atype,pos])

        #atomic species
        self.natypes = int(self.data.findall("output/atomic_species")[0].get('ntyp'))
        atypes = self.data.findall("output/atomic_species/species")
        self.atomic_species = {}
        for i in range(self.natypes):
            atype_string = atypes[i].get('name').strip()
            atype_mass = atypes[i].findall('mass')[0].text.strip()
            atype_pseudo = atypes[i].findall('pseudo_file')[0].text.strip()
            self.atomic_species[atype_string]=[atype_mass,atype_pseudo]

        #number of kpoints and bands
        self.nkpoints = int(self.data.findall('output/band_structure/nks')[0].text.strip())
        self.nbands = int(self.data.findall('output/band_structure/nbnd')[0].text.strip())

        #total energy
        self.energy = float(self.data.findall('output/total_energy/etot')[0].text)

        #fermi level
        self.fermi = float(self.data.find('output/band_structure/highestOccupiedLevel').text)

        kstates = self.data.findall('output/band_structure/ks_energies')

        #lists with the kpoints and the associated weights,
        #arrays with the ks energies and with the occupations
        self.kpoints = []
        self.weights = []
        self.evals = []
        self.occupations = []
        for k in kstates:
            kpoint = [float(x) for x in k.findall('k_point')[0].text.strip().split()]
            weight = [float(x) for x in k.findall('k_point')[0].get('weight').strip().split()]
            eval = [float(x) for x in k.findall('eigenvalues')[0].text.strip().split()]
            occ = [float(x) for x in k.findall('occupations')[0].text.strip().split()]
            self.kpoints.append(kpoint)
            self.weights.append(weight)
            self.evals.append(eval)
            self.occupations.append(occ)
        self.evals = np.array(self.evals)
        self.occupations = np.array(self.occupations)

        #spin related properties and spin-orbit coupling
        lsda = self.data.findall('output/band_structure/lsda')[0].text.strip()
        self.lsda = True if lsda == 'true' else False
        noncolin = self.data.findall('output/band_structure/noncolin')[0].text.strip()
        self.noncolin = True if noncolin == 'true' else False
        spinorbit = self.data.findall('output/band_structure/spinorbit')[0].text.strip()
        self.spinorbit = True if spinorbit == 'true' else False
        self.spin_degen = 1 if self.lsda or self.noncolin else 2

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
        Return the fermi energy the system. If convert_eV is True the energy
        is provided in eV other the Hartree units are used.
        """
        if convert_eV:
            return HaToeV*self.fermi
        else:
            return self.fermi

    def get_num_occupied_bands(self):
        """
        Compute the number of occupied bands of the system. The method check if the number is
        equal for all the kpoints.

        """
        num_occupied = []
        for occupation in self.occupations:
            num_occupied.append((occupation == 1.).sum())
        if all(occ == num_occupied[0] for occ in num_occupied):
            return num_occupied[0]
        else:
            print('number of occupied bands is k-dependent')
            return None

    def get_evals(self, set_gap=None, set_direct_gap=None):
        """
        Return the ks energies for each kpoint (in eV). The top of the valence band is used as the
        reference energy value. It is possible to shift the energies of the empty bands by setting an arbitrary
        value for the gap (direct or indirect). Implemented and tested only for semiconductors.

        Args:
            set_gap (float) : target gap in eV
            set_direct_gap (float) : target direct_gap in eV. If set_gap is provided this
                parameter is ignored

        """
        evals = self.evals*HaToeV
        homo_band = evals[:,self.get_num_occupied_bands()-1]
        VBM = homo_band.max()
        evals -= VBM

        if self.get_num_occupied_bands() == self.nbands: # only occupied bands are present
            return evals
        else: # shift the energy level of the empty bands if needed
            lumo_index = self.get_num_occupied_bands()
            gap = self.get_gap(verbose=False)
            scissor = 0.
            if set_gap is not None: scissor = set_gap - gap['gap']
            elif set_direct_gap is not None: scissor = set_direct_gap - gap['direct_gap']
            print('Apply a scissor of',scissor,'eV')
            evals[:,lumo_index:] += scissor
            return evals

    def get_gap(self,verbose=True):
        """
        Compute the energy gap of the system (in eV). The method check if the gap is direct or
        indirect. Implemented and tested only for semiconductors.

        Return:
            :py:class:`dict` : a dictionary with the values of direct and indirect gaps and the positions
            of the VMB and CBM

        """
        import numpy as np
        homo_band = self.evals[:,self.get_num_occupied_bands()-1]
        lumo_band = self.evals[:,self.get_num_occupied_bands()]

        VBM = homo_band.max()
        position_vbm = homo_band.argmax()
        CBM = lumo_band.min()
        position_cbm = lumo_band.argmin()
        gap = (CBM-VBM)*HaToeV
        direct_gap = (lumo_band[position_vbm]-homo_band[position_vbm])*HaToeV

        # If there are several copies of the same point on the path it can happen that
        # the system is recognized as an indirect gap for numerical noise, so we check
        if np.allclose(gap,direct_gap,atol=1e-5,rtol=1e-5):
            gap = direct_gap
            position_cbm = position_vbm
        if verbose:
            if position_cbm == position_vbm:
                print('Direct gap system')
                print('=================')
                print('Gap :',gap,'eV')
            else :
                print('Indirect gap system')
                print('===================')
                print('Gap :',gap,'eV')
                print('Direct gap :',direct_gap,'eV')
        return {'gap':gap,'direct_gap':direct_gap,'position_cbm':position_cbm,'positon_vbm':position_vbm}

    def Dos(self,Emin=-20, Emax=20, deltaE=0.001, deg=0.00):
        """
        Compute the DOS.
        DOS(E) = sum_{n,K} [delta(E - E_{n}(K)) * weight(K)]

        Todo:
            Transform this function into a class. In this way one could compute
            a dos also starting from something different from a PwParser object.
            Moreover, it could be possible to combine the data of several different
            object into a single dos.

        Note:
            What happens if we do not specify the weight of the kpoints????

        Args:
            Emin:   Starting energy for the DOS (in eV)
            Emax:   Final energy for the DOS (in eV)
            deltaE: Tick separation on the X axis
            deg:    Sigma to be used for a gaussian broadening.
		            Default = 0.0: Does not apply any broadening.

        Return:
            :py:class:`array`: The first column runs over the energies, the second
            one contain the corresponding value of the DOS
        """
        import numpy as np
        #from ..tools.broad import broad
        res = np.linspace(Emin, Emax, (Emax-Emin)/deltaE+1).reshape(1,-1)
        res = np.pad(res, ((0,1),(0,0)), 'constant')
        energies = self.evals*HaToeV - self.get_fermi()

        for n,egv in enumerate(energies):
            i = np.floor((egv - Emin) / deltaE +0.5).astype(dtype='int')
            i = i[np.where( (0 <= i) & (i < res[0].size))]
            res[1,i] += self.weights[n]
        res[1:] /= deltaE

        if deg > 0:
            res = broad(res, t='gauss', deg=deg, axis=1)

        return res.T
