"""
Class to perform the parsing of a Yambo o- file(s).
"""

import numpy as np

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
    Read the filename and append all the lines that do not start
    with skip to a list

    Args:
        filename(str) : name of the file
        skip(str) : first elements of the skipped lines
    """
    lines = []
    with open(filename) as f:
        for l in f:
            if not l.startswith(skip): lines.append(l)
    return lines

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

    columns = np.array(splitted).transpose()
    return columns

def make_dict(columns,suffix):
    """
    Create a dictionary from the columns array. If the suffix is found in the
    reference_column_names attribute to the keys the associated names, otherwise
    associate string value 'col'+str(ind), where .
    """
    data = {}
    for ind,col in enumerate(columns):
        if suffix in reference_column_names:
            key = reference_column_names[suffix][ind]
        else:
            key = 'col'+str(ind)
            #key = str(ind)
        data[key] = col
    return data

class YamboParser(dict):

    def __init__(self,files,verbose=False):
        """
        Initialize the data member of the class.

        Args:
            files(list): The list of strings with the names of the file to be parsed.
        """
        dict.__init__(self)
        for file in files:
            suffix = file.rsplit('.')[-1]
            self[suffix] = {}
            if verbose: print('Parse file',file)
            self.parseYamboOutput(file,suffix)

    def parseYamboOutput(self,file,suffix):
        """
        Read the data from the o- file. Data of the file are stored as key : values
        in the self[suffix] dictionary. The names of the keys are taken from the
        reference_column_names, if the suffix is recognized.
        """
        lines = file_to_list(file)
        columns = eval_columns(lines)
        self[suffix] = make_dict(columns,suffix)
