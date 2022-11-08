"""
Module that manages the parsing of a Yambo o- file(s).
"""

import numpy as np

# Specifies the name of the columns of the o- files for various type of runs. There are
# two distint dictionaries because the qp labels depend on the extendOut option.

rt_column_names = {
    'carriers' : ['time','dnhmne','dnh','dne'],
    'currents' : ['time','j_x','j_y','j_z'],
    'polarization' : ['time','Pol_x','Pol_y','Pol_z'],
    'spin_magnetization' :
        ['time','Ms_x','Ms_y','Ms_z','Mv_x','Mv_y','Mv_z','Mc_x','Mc_y','Mc_z'],
    'orb_magnetization' :
        ['time','Ml_x','Ml_y','Ml_z','Mi_x','Mi_y','Mi_z'],
    'external_field' :
        ['time','Ex_Re','Ey_Re','Ez_Re','Ex_Im','Ey_Im','Ez_Im','Profile','Intensity','Fluence']
}
hf_column_names = {
    'hf' : ['kpoint','band','E0','Ehf','Vxc','Vnlxc']
}
qp_column_names = {
    'qp' : ['kpoint','band','E0','EmE0','sce0']
}
qp_column_names_extendOut = {
    'qp' : ['kpoint','band','E0','E','EmE0','Dft','hf','sce0','sce','dsc_dwe0','z_Re','z_Im','width_mev','width_fs']
}
eps_column_names = {
    'eps_q1_ip' : ['energy','eps_imag','eps_real'],
    'eps_q1_diago_bse' : ['energy','eps_imag','eps_real','eps_o_imag','eps_o_real'],
    'eps_q1_haydoc_bse' : ['energy','eps_imag','eps_real','eps_o_imag','eps_o_real','eps_p_imag','eps_p_real']

}
reference_column_names = {**rt_column_names,**hf_column_names,**qp_column_names,**eps_column_names}
reference_column_names_extendOut = {**rt_column_names,**hf_column_names,**qp_column_names_extendOut,**eps_column_names}

def file_to_list(filename,skip='#'):
    """
    Read the filename and append all the lines that do not start
    with the skip string, to a list.

    Args:
        filename (str): name of the file
        skip (str): first elements of the skipped lines
    """
    lines = []
    with open(filename) as f:
        for l in f:
            if not l.startswith(skip): lines.append(l)
    return lines

def _floats_from_string(line):
  """
  Split a string using blank spaces and convert the elements to float. If an element
  cannot be converted it is skipped.
  """
  line_float = []
  for value in line.split():
      try: line_float.append(float(value))
      except ValueError: pass
  return line_float

def build_columns(lines):
    """
    Split each line of the output of file_to_list into a list and convert
    its elements to float. The procedure deletes the values that cannot be converted
    to float, for istance the string that specifies the high-symmetry points in the
    ypp bands_interpolated post-processing.
    Then transpose the array so that each element is a column of the data of the file.
    """
    splitted = []
    for line in lines:
         splitted.append(_floats_from_string(line))

    columns = np.array(splitted).transpose()
    return columns

def make_dict(columns,suffix,extendOut):
    """
    Create a dictionary from the columns array. If the suffix is found in the
    ref dictionary attribute to the keys the associated names, otherwise
    associate string value 'col'+str(ind), where ind is the column index starting
    from zero. The choice of the ref dictionary depends on the value of extendOut.

    Args:
        columns (:py:class:`array`) : array with the data sorted in columns
        suffix (string) : specifies the run level
        extendOut (bool) : specifies which dictionary has to be used as reference
            values of the columns names
    """
    if extendOut:
        ref = reference_column_names_extendOut
    else:
        ref = reference_column_names
    data = {}
    for ind,col in enumerate(columns):
        if suffix in ref:
            key = ref[suffix][ind]
        else:
            key = 'col'+str(ind)
        data[key] = col
    return data

def parseYamboOutput(file,suffix,extendOut):
    """
    Read the data from the o- file. Data of the file are stored as key : values
    in the self[suffix] dictionary. The names of the keys are taken from the
    reference_column_names or from the reference_column_names_extendOut (depending on
    the value of the boolean extendOut), if the suffix is recognized.
    """
    lines = file_to_list(file)
    columns = build_columns(lines)
    return make_dict(columns,suffix,extendOut)

class YamboOutputParser(dict):
    """
    Class that performs the parsing of a Yambo o- file(s). The class ineriths from :py:class:`dict`
    and the instance of the class is a dictionary with the data. The keys correspond to the extension
    of the parsed files.

    Args:
        output (:py:class:`list`): Dictionary with the structure of the output of :py:meth:`get_output_files`
            of the `YamboCalculator` module.
        verbose (:py:class:`boolean`) : Determine the amount of information provided on terminal
        extendOut (:py:class:`boolean`) : Determine which dictionary is used as reference for the
                        names of the variables
    """

    def __init__(self, output, verbose=True,extendOut=True):
        """
        Initialize the data member of the class.
        """
        dict.__init__(self)
        for suffix,file in output.items():
            if verbose: print('Parse file',file)
            self[suffix] = parseYamboOutput(file,suffix,extendOut)

    @classmethod
    def from_path(cls, path, verbose = True, extendOut = True):
        """
        Init the a :class:`YamboOutputParser` instance using the 'o-' files found inside the path. If several replica
        of the output files are found the method select the one associated to the last Yambo computation.

        Args:
            path (:py:class:`string`): name of the folder that contains the 'o-' files
            verbose (:py:class:`boolean`) : Determine the amount of information provided on terminal
            extendOut (:py:class:`boolean`) : Determine which dictionary is used as reference for the
                            names of the variables
        """
        from mppi.Calculators.YamboCalculator import get_output_files

        files = get_output_files(path)
        return cls(files,verbose=verbose,extendOut=extendOut)

    @classmethod
    def from_file(cls, file, verbose = True, extendOut = True):
        """
        Init the a :class:`YamboOutputParser` instance using a single 'o-' file. The key of the dictionary
        built by the parser is computed using the :py:meth:`type_identifier` method of the `YamboCalculator` module.

        Args:
            file (:py:class:`string`): name of the 'o-' file
            verbose (:py:class:`boolean`) : Determine the amount of information provided on terminal
            extendOut (:py:class:`boolean`) : Determine which dictionary is used as reference for the
                            names of the variables
        """
        from mppi.Calculators.YamboCalculator import type_identifier

        key = type_identifier(file)
        output =  {key : file}
        return cls(output,verbose=verbose,extendOut=extendOut)

    def get_info(self):
        """
        Provide information on the keys structure of the instance of the class
        """
        print('YamboOutputParser variables structure')
        for key,value in self.items():
            print('suffix',key,'with',value.keys())

    def get_energy(self,k,bnd,verbose=False):
        """
        Compute the energy (in eV) of the selected state with k and bnd indexes.
        The method is implemented for the 'hf' and 'qp' runleves. In the first case
        it seeks for the 'Ehf' variable while in the former it looks for the 'E' variable
        (written only with the extendOut option enabled!).
        In the other cases a warning is created.

        Args:
            k (:py:class:`int`): k-point index
            bnd (:py:class:`int`) : band index

        Return:
            (:py:class:`float`) : the value of the energy in eV

        """
        energy_value = {'hf':'Ehf','qp':'E'}
        runlevel = list(self.keys())[0]
        if runlevel not in energy_value:
            print('Actual runlevel %s not implemented for energy computation!'%runlevel)
            return None
        if verbose : print('Compute energy of the state k=%s,b=%s for the %s runlevel.'%(k,bnd,runlevel))
        kpoint = self[runlevel]['kpoint']
        band = self[runlevel]['band']
        index = -1
        for ind, t  in enumerate(zip(kpoint,band)):
            if t == (k,bnd) : index = ind
        if index == -1:
            print('k-points and/or band indexes not found!')
            return None
        energy = self[runlevel][energy_value[runlevel]][index]
        return energy

    def get_gap(self,k_full,band_full,verbose=False,**kwargs):
        """
        Compute the energy gap of the selected (k_full,band_full) and
        (k_empty,band_empty) couples (in eV).

        Args:
            k_full (:py:class:`int`): k-point index of the full state
            band_full (:py:class:`int`) : band index of the full state
            **kwargs : these parameters allow the user to set the k_empty and
                band_empty parameters. If not provided the values k_empty=k_full
                and band_empty=band_full+1 are used.

        Return:
            (:py:class:`float`) : the value of the band gap in eV

        """
        if 'k_empty' in kwargs : k_empty = kwargs['k_empty']
        else: k_empty = k_full
        if 'band_empty' in kwargs : band_empty = kwargs['band_empty']
        else: band_empty = band_full+1
        E_full = self.get_energy(k_full,band_full,verbose=verbose)
        E_empty = self.get_energy(k_empty,band_empty,verbose=verbose)
        energy_gap = E_empty-E_full
        if verbose :
            print('(kpoint,band) indexes of the full state:',(k_full,band_full))
            print('(kpoint,band) indexes of the empty state:',(k_empty,band_empty))
            print('Energy gap = %s eV'%energy_gap)
        return energy_gap
