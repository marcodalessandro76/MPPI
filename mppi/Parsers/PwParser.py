"""
Class to perform the parsing of a QuantumESPRESSO XML file. Makes usage of the
data-file-schema.xml file that is found in the run_dir/prefix.save folder.
"""

HaToeV = 27.211386

class PwParser():

    def __init__(self,file,verbose=True):
        """
        Initialize the data member of the class. If the parsing of the XML file
        does not succeeds the data attribute is set to None.

        Args:
            file(str): The name, including the pathm of the data-file-schema.xml
        """
        self.file = file

        if verbose: print("Parse file : %s"%self.file)
        try:
            self.parseXML(self.file)
        except FileNotFoundError:
            if verbose: print('Failed to read %s'%self.file)
            self.data = None

    def parseXML(self,file):
        """
        Read the data from the xml file in the new format of quantum espresso.
        Some variable are extracted from the XML file and stored in the attribute
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

    def get_kpath(self):
        """
        Assumes that self.kpoints contains a sampling of kpoints along a path.
        Compute the curvilinear ascissa along the path.

        Returns:
        kpath(array) : array the value of the curvilinear ascissa along the path
        """
        import numpy as np
        kpoints = np.array(self.kpoints)
        kpath = [0]
        distance = 0
        for nk in range(1,len(kpoints)):
            distance += np.linalg.norm(kpoints[nk-1]-kpoints[nk])
            kpath.append(distance)
        return np.array(kpath)

    def get_bands(self,convert_eV=True):
        """
        Convert the array self.evals into the the array bands, where bands[i] gives
        the energies of the i-th band along the path.
        The fermi level is used as the reference energy and the results are
        expressed in eV if convert_eV = True
        """
        import numpy as np
        bands = []
        for b in range(len(self.evals[0])): #number of bands
            bands.append(self.evals[:,b])
        bands = np.array(bands) - self.fermi
        if convert_eV: bands *= HaToeV
        return bands

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
            What happens if I do not specify the weight of the kpoints????

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
