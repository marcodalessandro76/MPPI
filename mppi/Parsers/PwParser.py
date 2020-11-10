"""
Class to perform the parsing of a QuantumESPRESSO XML file. Makes usage of the
data-file-schema.xml file that is found in the run_dir/prefix.save folder.
"""

# def _compute_transitions(bands,in_list,fin_list):
#     """
#     Compute the (positive) transition energies for the bands (on a single kpoint)
#
#     Args:
#         bands (list) : list with the energies
#         in_list (list) : indexes of the bands used as starting points of the transitions
#         fin_list (list) : indexes of the bands used as final points of the transitions
#
#     Returns:
#         transitions (list) : list with the transition energies for each possible couple
#         of (distinct) in and out bands
#     """
#     transitions = []
#     for v in in_list:
#         for c in fin_list:
#             if c > v:
#                 transitions.append(bands[c]-bands[v])
#     return transitions

from mppi.Utilities import HaToeV
from mppi.Parsers import Functions as F

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
        num_electrons : number of electrons
        nkpoints : numer of kpoints
        nbands : number of bands
        nbands_valence : number of occupied bands (for systems with a gap)
        nbands_conduction : number of empty bands (for systems with a gap)
        occupations_kind : type of occupation (fixed or smearing)
        kpoints : list of the kpoints
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
        self.nbands_valence = int(self.num_electrons/self.spin_degen)
        self.nbands_conduction = self.nbands-self.nbands_valence

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

    def get_evals(self, set_gap = None, set_direct_gap = None):
        """
        Return the ks energies for each kpoint (in eV). The top of the valence band is used as the
        reference energy value. It is possible to shift the energies of the empty bands by setting an arbitrary
        value for the gap (direct or indirect).
        Implemented only for semiconductors, the energy shift of empty bands does not update their occupation levels.

        Args:
            set_gap (float) : set the value of the gap (in eV) of the system
            set_direct_gap (float) : set the value of the direct gap (in eV) of the system. If set_gap
                            is provided this parameter is ignored

        Return:
            :py:class:`numpy.array`  : an array with the ks energies for each kpoint

        """
        evals = F.get_evals(self.evals,self.nbands,self.nbands_valence,set_gap=set_gap,set_direct_gap=set_direct_gap)
        return evals

    # def get_evals(self, set_gap = None, set_direct_gap = None):
    #     """
    #     Return the ks energies for each kpoint (in eV). The top of the valence band is used as the
    #     reference energy value. It is possible to shift the energies of the empty bands by setting an arbitrary
    #     value for the gap (direct or indirect).
    #     Implemented only for semiconductors, the energy shift of empty bands does not update their occupation levels.
    #
    #     Args:
    #         set_gap (float) : set the value of the gap (in eV) of the system
    #         set_direct_gap (float) : set the value of the direct gap (in eV) of the system. If set_gap
    #                         is provided this parameter is ignored
    #
    #     Return:
    #         :py:class:`numpy.array`  : an array with the ks energies for each kpoint
    #
    #     """
    #     evals = self.evals*HaToeV
    #     homo_band = evals[:,self.nbands_valence-1]
    #     VBM = homo_band.max()
    #     evals -= VBM
    #
    #     if set_gap == None and set_direct_gap == None:
    #         return evals
    #     elif self.nbands_valence == self.nbands: # only occupied bands are present
    #         print('There are no empty bands. `Set gap` has not been applied')
    #         return evals
    #     else: # shift the energy level of the empty bands if needed
    #         gap = self.get_gap(verbose=False)
    #         scissor = 0.
    #         if set_gap is not None: scissor = set_gap - gap['gap']
    #         elif set_direct_gap is not None: scissor = set_direct_gap - gap['direct_gap']
    #         print('Apply a scissor of',scissor,'eV')
    #         evals[:,self.nbands_valence:] += scissor
    #         return evals

    def get_transitions(self, initial = 'full', final = 'empty',set_gap = None, set_direct_gap = None):
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
            set_gap (float) : set the value of the gap (in eV) of the system
            set_direct_gap (float) : set the value of the direct gap (in eV) of the system. If set_gap
                            is provided this parameter is ignored

        Return:
            :py:class:`numpy.array`  : an array with the transition energies for each kpoint

        """
        transitions = F.get_transitions(self.evals,self.nbands,self.nbands_valence,initial=initial,final=final,set_gap=set_gap,set_direct_gap=set_direct_gap)
        return transitions

    # def get_transitions(self, initial = 'full', final = 'empty',set_gap = None, set_direct_gap = None):
    #     """
    #     Compute the (vertical) transitions energies. For each kpoint compute the transition energies, i.e.
    #     the (positive) energy difference (in eV) between the final and the initial states.
    #
    #     Args:
    #         initial (string or list) : specifies the bands from which electrons can be extracted. It can be set to `full` or
    #             `empty` to select the occupied or empty bands, respectively. Otherwise a list of bands can be
    #             provided
    #         final  (string or list) : specifies the final bands of the excited electrons. It can be set to `full` or
    #             `empty` to select the occupied or empty bands, respectively. Otherwise a list of bands can be
    #             provided
    #         set_gap (float) : set the value of the gap (in eV) of the system
    #         set_direct_gap (float) : set the value of the direct gap (in eV) of the system. If set_gap
    #                         is provided this parameter is ignored
    #
    #     Return:
    #         :py:class:`numpy.array`  : an array with the transition energies for each kpoint
    #
    #     """
    #     import numpy as np
    #     if initial == 'full':
    #         in_list = [ind for ind in range(self.nbands_valence)]
    #     elif initial == 'empty':
    #         in_list = [ind for ind in range(self.nbands_valence,self.nbands)]
    #     else:
    #         in_list = initial
    #     if final == 'full':
    #         fin_list = [ind for ind in range(self.nbands_valence)]
    #     elif final == 'empty':
    #         fin_list = [ind for ind in range(self.nbands_valence,self.nbands)]
    #     else:
    #         fin_list = final
    #
    #     evals = self.get_evals(set_gap=set_gap,set_direct_gap=set_direct_gap)
    #     transitions = []
    #     for bands in evals:
    #         transitions.append(_compute_transitions(bands,in_list,fin_list))
    #     transitions = np.array(transitions)
    #
    #     return transitions

    def get_gap(self, verbose = True):
        """
        Compute the energy gap of the system (in eV). The method check if the gap is direct or
        indirect. Implemented and tested only for semiconductors.

        Return:
            :py:class:`dict` : a dictionary with the values of direct and indirect gaps and the positions
            of the VMB and CBM

        """
        gap = F.get_gap(self.evals,self.nbands_valence,verbose=verbose)
        return gap

    # def get_gap(self, verbose = True):
    #     """
    #     Compute the energy gap of the system (in eV). The method check if the gap is direct or
    #     indirect. Implemented and tested only for semiconductors.
    #
    #     Return:
    #         :py:class:`dict` : a dictionary with the values of direct and indirect gaps and the positions
    #         of the VMB and CBM
    #
    #     """
    #     import numpy as np
    #     homo_band = self.evals[:,self.nbands_valence-1]
    #     try:
    #         lumo_band = self.evals[:,self.nbands_valence]
    #     except IndexError:
    #         print('There are no empty states. Gap cannot be computed.')
    #         return None
    #
    #     VBM = homo_band.max()
    #     position_vbm = homo_band.argmax()
    #     CBM = lumo_band.min()
    #     position_cbm = lumo_band.argmin()
    #     gap = (CBM-VBM)*HaToeV
    #     direct_gap = (lumo_band[position_vbm]-homo_band[position_vbm])*HaToeV
    #
    #     # If there are several copies of the same point on the path it can happen that
    #     # the system is recognized as an indirect gap for numerical noise, so we check
    #     if np.allclose(gap,direct_gap,atol=1e-5,rtol=1e-5):
    #         gap = direct_gap
    #         position_cbm = position_vbm
    #     if verbose:
    #         if position_cbm == position_vbm:
    #             print('Direct gap system')
    #             print('=================')
    #             print('Gap :',gap,'eV')
    #         else :
    #             print('Indirect gap system')
    #             print('===================')
    #             print('Gap :',gap,'eV')
    #             print('Direct gap :',direct_gap,'eV')
    #     return {'gap':gap,'direct_gap':direct_gap,'position_cbm':position_cbm,'positon_vbm':position_vbm}
