"""
Class to perform the parsing of a Yambo o- file(s).
"""

import numpy as np

# Specifies the name of the columns of the o- files for various type of runs. There are
# two distint dictionaries depending if the ExtendOut option has been activated or not.

reference_column_names_extendOut = {
    'hf' : ['kpoint','band','e0','ehf','dft','hf'],
    'qp' : ['kpoint','band','e0','e','eme0','dft','hf','sce0','sce','dsc_dwe0','z_Re','z_Im','width_mev','width_fs'],
    # complete the list for this run levels.....
    'carriers' : ['time','dnhmne','dnh','dne'],
    'currents' : ['time','j_x','j_y','j_z'],
    'polarization' : ['time','Pol_x','Pol_y','Pol_z'],
    'spin_magnetization' :
        ['time','Ms_x','Ms_y','Ms_z','Mv_x','Mv_y','Mv_z','Mc_x','Mc_y','Mc_z'],
    'orb_magnetization' :
        ['time','Ml_x','Ml_y','Ml_z','Mi_x','Mi_y','Mi_z'],
    'external_field' :
        ['time','Ex_Re','Ey_Re','Ez_Re','Ex_Im','Ey_Im','Ez_Im','Intensity','Fluence']
}

reference_column_names = {
    'hf' : ['kpoint','band','e0','ehf','dft','hf'],
    'qp' : ['kpoint','band','e0','eme0','sce0'],
    'carriers' : ['time','dnhmne','dnh','dne'],
    'currents' : ['time','j_x','j_y','j_z'],
    'polarization' : ['time','Pol_x','Pol_y','Pol_z'],
    'spin_magnetization' :
        ['time','Ms_x','Ms_y','Ms_z','Mv_x','Mv_y','Mv_z','Mc_x','Mc_y','Mc_z'],
    'orb_magnetization' :
        ['time','Ml_x','Ml_y','Ml_z','Mi_x','Mi_y','Mi_z'],
    'external_field' :
        ['time','Ex_Re','Ey_Re','Ez_Re','Ex_Im','Ey_Im','Ez_Im','Intensity','Fluence']
}

def file_to_list(filename,skip='#'):
    """
    Read the filename and append all the lines, that do not start
    with the skip string, to a list.

    Args:
        filename(str) : name of the file
        skip(str) : first elements of the skipped lines
    """
    lines = []
    with open(filename) as f:
        for l in f:
            if not l.startswith(skip): lines.append(l)
    return lines

def floats_from_string(line):
  """
  Split a string using blank spaces and convert the elements to float.
  """
  line_float = []
  for value in line.split():
      try: line_float.append(float(value))
      except ValueError: pass
  return line_float

def eval_columns(lines):
    """
    Split each line of the output of file_to_list into a list and convert
    its elements to float. The procedure deletes the values that cannot be converted
    to float, for istance the string that specifies the high-symmetry points in the
    ypp bands_interpolated post-processing.
    Then transpose the array so that each element is a column of the data of the file.
    """
    splitted = []
    # simple implementation that fails if some elements cannot be converted to float
    # for l in lines: # run over lines
    #     splitted.append(list(map(float,l.split())))

    for l in lines: #run over lines
        l_float = [] #each line is converted in a string of float
        for value in l.split():
            # if some element cannot converted to float are ignored
            try: l_float.append(float(value))
            except ValueError: pass
        splitted.append(l_float)

    # for line in lines:
    #     splitted.append(floats_from_string(line))

    columns = np.array(splitted).transpose()
    return columns

def make_dict(columns,suffix,extendOut):
    """
    Create a dictionary from the columns array. If the suffix is found in the
    ref dictionary attribute to the keys the associated names, otherwise
    associate string value 'col'+str(ind), where ind is the column index starting
    from zero. The choice of the ref dictionary depends on the value of extendOut.

    Args:
        columns (:numpy:class:`array`) : array with the data sorted in columns
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

class YamboParser(dict):

    def __init__(self,files,verbose=False,extendOut=True):
        """
        Initialize the data member of the class.

        Args:
            files(list): The list of strings with the names of the file to be parsed.
            verbose (bool) : Determine the amount of information provided on terminal.
            extendOut (bool) : Determine which dictionary is used as reference for the
                            name of the variables
        """
        dict.__init__(self)
        for file in files:
            suffix = file.rsplit('.')[-1]
            self[suffix] = {}
            if verbose: print('Parse file',file)
            self.parseYamboOutput(file,suffix,extendOut)

    def parseYamboOutput(self,file,suffix,extendOut):
        """
        Read the data from the o- file. Data of the file are stored as key : values
        in the self[suffix] dictionary. The names of the keys are taken from the
        reference_column_names or from the reference_column_names_extendOut (depending on
        the value of the boolean extendOut), if the suffix is recognized.
        """
        lines = file_to_list(file)
        columns = eval_columns(lines)
        self[suffix] = make_dict(columns,suffix,extendOut)
