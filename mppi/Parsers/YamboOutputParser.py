"""
Module that manages the parsing of a Yambo o- file(s).
"""

import numpy as np

# Specifies the name of the columns of the o- files for various type of runs. There are
# two distint dictionaries depending if the ExtendOut option has been activated or not.

# The rt outputs are not modified by the extendOut option
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

reference_column_names_extendOut = {
    'hf' : ['kpoint','band','e0','ehf','dft','hf'],
    'qp' : ['kpoint','band','e0','e','eme0','dft','hf','sce0','sce','dsc_dwe0','z_Re','z_Im','width_mev','width_fs'],
}
reference_column_names_extendOut.update(rt_column_names)

reference_column_names = {
    'hf' : ['kpoint','band','e0','ehf','dft','hf'],
    'qp' : ['kpoint','band','e0','eme0','sce0'],
}
reference_column_names.update(rt_column_names)

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

def files_from_folder(path):
    """
    Scan the files in the folder and build a list with the names of all the files
    that contain the 'o-' term in their name.

    Args:
        path (string) : name of the folder
    """
    import os
    listdir= os.listdir(path)
    ofiles = []
    for file in listdir:
        if 'o-' in file:
            ofiles.append(os.path.join(path,file))
    return ofiles

class YamboOutputParser(dict):
    """
    Class that performs the parsing of a Yambo o- file(s). The class ineriths from :py:class:`dict`
    and the instance of the class is a dictionary with the data. The keys correspond to the extension
    of the parsed files

    Args:
        files (:py:class:`list`): The list of strings with the names of the file to be parsed
        verbose (:py:class:`boolean`) : Determine the amount of information provided on terminal
        extendOut (:py:class:`boolean`) : Determine which dictionary is used as reference for the
                        names of the variables

    """

    def __init__(self,files,verbose=True,extendOut=True):
        """
        Initialize the data member of the class.
        """
        dict.__init__(self)
        for file in files:
            suffix = file.rsplit('.')[-1]
            if verbose: print('Parse file',file)
            self.parseYamboOutput(file,suffix,extendOut)
            self[suffix] = self.parseYamboOutput(file,suffix,extendOut)

    @classmethod
    def from_path(cls,path,verbose = False, extendOut = True):
        """
        Init the a :class:`YamboOutputParser` instance using all the 'o-' files found inside the path.

        Args:
            path (:py:class:`string`): name of the folder that contains the 'o-' files
            verbose (:py:class:`boolean`) : Determine the amount of information provided on terminal
            extendOut (:py:class:`boolean`) : Determine which dictionary is used as reference for the
                            names of the variables
        """
        files = files_from_folder(path)
        return cls(files,verbose=verbose,extendOut=extendOut)

    def parseYamboOutput(self,file,suffix,extendOut):
        """
        Read the data from the o- file. Data of the file are stored as key : values
        in the self[suffix] dictionary. The names of the keys are taken from the
        reference_column_names or from the reference_column_names_extendOut (depending on
        the value of the boolean extendOut), if the suffix is recognized.
        """
        lines = file_to_list(file)
        columns = build_columns(lines)
        return make_dict(columns,suffix,extendOut)

    def get_info(self):
        """
        Provide information on the keys structure of the instance of the class
        """
        print('YamboOutputParser variables structure')
        for key,value in self.items():
            print('suffix',key,'with',value.keys())
