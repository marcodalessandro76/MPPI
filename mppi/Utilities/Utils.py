"""
This file contains some low-level useful functions.
The module can be loaded in the notebook in one of the following way

from mppi import Utilities as U

from mppi.Utilities import build_kpath,...
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
