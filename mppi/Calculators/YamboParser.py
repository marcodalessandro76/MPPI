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
    Seek for the line that contains the string 'K-points' or '|k|' and build a
    list with the names of the columns of data.
    """
    line_keys = ''
    with open(fname) as f:
        for l in f:
            if 'K-point' in l or '|k|' in l:
                line_keys = l
                break
    keys = line_keys.split()
    return keys

def _purge_values(keys):
    """
    Remove the strings #,(a.u.),(rlu),(alat),(cc) (if present) from the
    list tha contains the keys.
    """
    purge_list = ['#','(a.u.)','(rlu)','(alat)','(cc)']
    for k in reversed(keys):
        if k in purge_list:
            keys.remove(k)
    return keys

def dict_parser(fname):
    """
    Build the dictionary from the output file in the form key:value.
    """
    keys = _build_keys(fname)
    keys = _purge_values(keys)
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
