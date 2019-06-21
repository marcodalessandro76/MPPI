"""
This module performs a parsing pf the Yambo output file. The method dict_parser
builds a dictionary with the results. The keys are read from the line that contains
'K-point' (for hf or qp output files) or '|k|' for the ypp outputs.
The class AttributeDict convert the dictionary in a object and allows us to access
to its attribute in the form AttributeDict.attr.
"""

def _parserArrayFromFile(fname):
    """"
    Build a list that contains the lines of fname avoiding the ones that start
    with #
    """
    lines = []
    with open(fname) as f:
        for l in f:
            if not l.startswith('#'):
                lines.append(l)
    #split each line in a list (of strings)
    larray = [[] for i in range(len(lines))]
    for ind,l in enumerate(lines):
        larray[ind] = l.split()
    #convert the string to double. If some elements is a string (it can happen in
    #the 4.4 Yambo version in the output file of a bands calculation) remove it
    for row in range(len(larray)):
        for col in range(len(larray[row])):
            try:
                larray[row][col] = float(larray[row][col])
            except ValueError:
                del larray[row][col]
    return larray

def _build_keys(fname):
    """
    Seek for the line that contains the proper string and build a list with the
    names of the columns of data. The value of the string to look for depend on
    the suffix of fname
    """
    suffix = fname.split('.')[-1]
    hook_string = None
    # hf or qp yambo computation
    if suffix in ['hf','qp']: hook_string = 'K-point'
    # yambo optics computation
    if suffix in ['eps_q1_ip','eel_q1_ip']: hook_string = 'E/ev'
    # ypp computation
    if suffix in  ['bands_interpolated','magnetization_x',\
                   'magnetization_y','magnetization_z']: hook_string = '|k|'
    # real time yambo_rt computation
    if suffix in ['carriers','external_field','current',\
                  'polarization'] : hook_string = 'Time[fs]'
    if suffix in ['magnetization'] : hook_string = 'Ms_x'
    # ypp_rt computations
    if suffix in ['YPP-RT_occupations_DATA','YPP-RT_occupations_dn_DATA',\
                  'YPP-RT_occupations_up_DATA'] : hook_string = 'E [eV]'
    line_keys = ''
    with open(fname) as f:
        for l in f:
            if hook_string in l :
                line_keys = l
                break
    keys = line_keys.split()
    return keys

def _clean_keys(keys):
    """
    Remove some spurious elements (if present) from the list that contains the
    keys. Some names are changed for better usability.
    """
    purge_list = ['#','(a.u.)','(rlu)','(alat)','(cc)','[eV]','@','[fs]']
    for k in reversed(keys):
        if k in purge_list:
            keys.remove(k)
    old_names = ['|k|','Time[fs]','E/ev[1]','EPS-Im[2]','EPS-Re[3]','EEL-Im[2]','EEL-Re[3]']
    new_names = ['k','time','E','Im','Re','Im','Re']
    for ind,k in enumerate(keys):
        if k in old_names:
            keys[ind] = new_names[old_names.index(k)]
    return keys

def dict_parser(fname):
    """
    Build a dictionary from the yambo output file in the form key:value.
    """
    keys = _build_keys(fname)
    keys = _clean_keys(keys)
    larray = _parserArrayFromFile(fname)
    results = {}
    for ind,key in enumerate(keys):
        results[key] = []
        for line in larray:
            results[key].append(line[ind])
    return(results)

class AttributeDict(object):
    """
    A class to convert a nested Dictionary into an object with key-values
    accessibly using attribute notation (AttributeDict.attribute) instead of
    key notation (Dict["key"]). This class recursively sets Dicts to objects,
    allowing you to recurse down nested dicts (like: AttributeDict.attr.attr)
    """
    def __init__(self, **entries):
        self.add_entries(**entries)

    def add_entries(self, **entries):
        for key, value in entries.items():
            if type(value) is dict:
                self.__dict__[key] = AttributeDict(**value)
            else:
                self.__dict__[key] = value

    def getAttributes(self):
        """
        Return all the attributes of the object
        """
        return self.__dict__.keys()
