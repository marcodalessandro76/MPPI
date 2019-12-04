"""
This file contains some low-level useful functions.
The module can be loaded in the notebook in one of the following way

import mppi.Utilities.Utils as U
from mppi.Utilities import Utils as U
from mppi.Utilities.Utils import build_kpath,...
"""

def fortran_bool(boolean):
    return {True:'.true.',False:'.false.'}[boolean]

def build_kpath(*kpoints,numstep=40):
    """
    Build a list of kpoints to be passed to the set_kpoints methods of the
    :class:`PwInput` for computing the band structure along a path.

    Args:
        klist(*list): identifies the high symmetry points of the k-path
        numstep(int): specifies the number of intermediate points used to build
                      the path
    """
    klist = []
    for k in kpoints[:-1]:
        klist.append(k+[numstep])
    klist.append(kpoints[-1]+[0])
    return klist

def name_from_id(id):
    """
    Convert the id into a run name. If id is a string, set name = id, if it is a
    dictionary build the name string of the run from the id dictionary.

    Args:
        id : id associated to the run
    Returns:
       name (str): name of the run associated to the dictionary ``id``
    """
    if type(id) is str :
        name = id
    elif type(id) is dict :
        keys=sorted(id.keys())
        name=''
        for k in keys:
            name += k+'_'+str(id[k])+'-'
        name = name.rstrip('-')
    else :
        print('id type not recognized')
        name = None
    return name

def names_from_id(id):
    """
    Convert the id into a list of run names to search with the function id_in_names
    and add the separator ',' to have the proper value of a key.
    """
    if id is None:
        return ['']
    else:
        return [name_from_id({k: v})+',' for k, v in id.items()]

def build_SAVE(source_dir,run_dir,command = 'p2y -a 2'):
    """
    Build the SAVE folder for a yambo computation.

    The method build the run_dir folder if it does not exists. Then check if the
    SAVE folder exists, if not runs the command in the source_dir and copy the
    SAVE folder in the run_dir. The option -a 2 ensures that labelling of the
    high-symmetry kpoints is consistent in both QE and Yambo.
    Lastly, executes yambo (without arguments) in the run_dir to build the r_setup.

    Note:
        Maybe is better to link the SAVE instead f copy it, to save space on disk.
    """
    import os
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)
        print('Create folder %s'%run_dir)
    if not os.path.isdir(source_dir):
        print('source folder',source_dir,'not found')
    elif os.path.isdir(os.path.join(run_dir,'SAVE')):
        print('SAVE folder already present in %s'%run_dir)
    else:
        # create the SAVE folder
        comm_str = 'cd %s; %s'%(source_dir,command)
        print('Executing command:', comm_str)
        os.system(comm_str)
        # copy the SAVE folder
        src = os.path.join(source_dir,'SAVE')
        comm_str = 'cp -r %s %s'%(src,run_dir)
        print('Executing command:', comm_str)
        os.system(comm_str)
        # build the r_setup
        comm_str = 'cd %s;OMP_NUM_THREADS=1 yambo'%run_dir
        print('Executing command:', comm_str)
        os.system(comm_str)

class Dict_to_attribute(object):
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
