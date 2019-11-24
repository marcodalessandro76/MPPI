"""
Class to perform the parsing of a QuantumESPRESSO XML file. Makes usage of the
data-file-schema.xml file that is found in the run_dir/prefix.save folder.
"""

HaToEv = 27.211386

class PwParser():

    def __init__(self,prefix,path='.',verbose=True):
        """
        Initialize the data member of the class. If the parsing of the XML file
        does not succeeds data contains the name of the XML file.
        """
        self.file = "%s/%s.save/data-file-schema.xml"%(path, prefix)

        if verbose: print("Parse file : %s"%self.file)
        try:
            self.parseXML(self.file)
        except FileNotFoundError:
            self.data = self.file
            print('Failed to read data-file-schema.xml in %s/%s.save'%(path,prefix))

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
        self.nkpoints = int(self.data.findall("output/band_structure/nks")[0].text.strip())
        self.nbands = int(self.data.findall("output/band_structure/nbnd")[0].text.strip())

        #total energy
        self.energy = float(self.data.findall('output/total_energy/etot')[0].text)

        #fermi level
        self.fermi = float(self.data.find('output/band_structure/highestOccupiedLevel').text)

        kstates = self.data.findall('output/band_structure/ks_energies')

        #list with the kpoints
        self.kpoints = []
        for i in range(self.nkpoints):
            kpoint = [float(x) for x in kstates[i].findall('k_point')[0].text.strip().split()]
            self.kpoints.append(kpoint)

        #list with the ks energies for each kpoint
        self.evals = []
        for k in kstates:
            eval = [float(x) for x in k.findall('eigenvalues')[0].text.strip().split()]
            self.evals.append(eval)
        self.evals = np.array(self.evals)

        #occupations of the ks states for each kpoint
        self.occupations = []
        for k in kstates:
            occ = [float(x) for x in k.findall('occupations')[0].text.strip().split()]
            self.occupations.append(occ)
        self.occupations = np.array(self.occupations)

    def get_energy(self,convert_eV = True):
        """
        Return the total energy the system. If convert_eV is True the energy
        is provided in eV other the Hartree units are used.
        """
        if convert_eV:
            return HaToEv*self.energy
        else:
            return self.energy

    def get_fermi(self,convert_eV = True):
        """
        Return the fermi energy the system. If convert_eV is True the energy
        is provided in eV other the Hartree units are used.
        """
        if convert_eV:
            return HaToEv*self.fermi
        else:
            return self.fermi 
