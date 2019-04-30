"""
This module performs a parsing pf the Yambo output file. It builds a dictionary with the results.
The keys are read from the line that contains "K-point"
"""

def _parserArrayFromFile(fname):
    """"
    Build a list that contains the lines of fname avoiding the ones that start with #
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
    #convert the string to double. If some elements is a string (it can happen in the 4.4 Yambo
    #version in the output file of a bands calculation) remove it
    for row in range(len(larray)):
        for col in range(len(larray[row])):
            try:
                larray[row][col] = float(larray[row][col])
            except ValueError:
                del larray[row][col]
    return larray

def _build_keys(fname):
    """
    Seek for the line that contains the string 'K-points' and build a list with the names of the
    columns of data.
    """
    line_keys = []
    with open(fname) as f:
        for l in f:
            if 'K-point' in l:
                line_keys.append(l)
                break
    keys = line_keys[0].split()[1:]
    return keys

def YamboParser(fname):
    """
    Build the dictionary from the output file in the form key:value.
    """
    keys = _build_keys(fname)
    larray = _parserArrayFromFile(fname)
    results = {}
    for ind,key in enumerate(keys):
        results[key] = []
        for line in larray:
            results[key].append(line[ind])
    return(results)
