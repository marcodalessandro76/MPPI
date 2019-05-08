"""
This module contains all the allowed pre_processing functions that can be called
before running a dataset.
"""

import os
#from .Datasets import *

# This function is already defined in Datasets. Its definition is repetead here
# with the underscore symbol because there is a problem of visibility from the
# Datasets module
def _name_from_id(id):
    """
    Hash the id into a run name
    Construct the name of the run from the id dictionary
    Args:
        id (dict): id associated to the run
    Returns:
       str: name of the run associated to the dictionary ``id``
    """
    keys=sorted(id.keys())
    name=''
    for k in keys:
        name += k+':'+str(id[k])+','
    return name.rstrip(',')

pre_processing_list = ['scf','nscf','yambo','break_sym']

def scf_pre_processing(run_dir):
    """
    Define the pre_processing function for a scf dataset.
    The method build the run_dir folder if it does not exists.
    """
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)
        print('Create folder %s'%run_dir)
    else :
        print('Folder %s already exists'%run_dir)

def nscf_pre_processing(run_dir,**kwargs):
    """
    Define the pre_processing function for a nscf dataset.
    The method build the run_dir folder if it does not exists. Then, for each id of
    the dataset copy the source_dir with name name_from_id(id).save.
    """
    source = kwargs['source_dir']
    ids = kwargs['ids']
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)
        print('Create folder %s'%run_dir)
    else :
        print('Folder %s already exists'%run_dir)

    for id in ids:
        dest = run_dir + '/' + _name_from_id(id) + '.save'
        # check if the source folder exists
        if not os.path.isdir(source): print('source folder : ',source,'not found')
        else :
            # copy the source folder only if the dest is not present
            if not os.path.isdir(dest):
                string = 'cp -r %s %s'%(source,dest)
                print('execute : ',string)
                os.system(string)
            else : print('SAVE folder already exsists. Source %s NOT COPIED'%source)

def yambo_pre_processing(run_dir,**kwargs):
    """
    Define the pre_processing function for a Yambo dataset.
    The method build the run_dir folder if it does not exists.Then check if the SAVE
    folder exists, if not runs p2y in the source_dir and copy the SAVE folder in the run_dir.
    Lastly, executes yambo (without arguments) in the run_dir to build the r_setup
    """
    source = kwargs['source_dir']
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)
        print('Create folder %s'%run_dir)
    if not os.path.isdir(source):
        print('source folder : ',source,' not found')
    elif os.path.isdir(run_dir +'/SAVE'):
        print('SAVE folder already present in %s'%run_dir)
    else:
        # run p2 y
        string = 'cd %s;p2y'%source
        print('execute : ',string)
        os.system(string)
        # copy the SAVE folder
        savef = source + '/SAVE'
        string = 'cp -r %s %s'%(savef,run_dir)
        print('execute : ',string)
        os.system(string)
        # build the r_setup
        string = 'cd %s;OMP_NUM_THREADS=1 yambo'%run_dir
        print('execute : ',string)
        os.system(string)
