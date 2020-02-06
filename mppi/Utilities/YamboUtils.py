"""
This file contains some useful functions to perform computation with Yambo.
The module can be loaded in the notebook in one of the following way

>>> from mppi import Utilities as U

>>> U.build_SAVE

or to load directly some elements

>>> from mppi.Utilities import build_SAVE

>>> build_SAVE

"""

def build_SAVE(source_dir,run_dir,command = 'p2y -a 2',make_link = True):
    """
    Build the SAVE folder for a yambo computation.

    The function creates the SAVE folder in the source_dir using the command provided
    as the command parameter (the option -a 2 ensures that labelling of the
    high-symmetry kpoints is consistent in both QE and Yambo) and create a symbolic
    link (or a copy) of the SAVE folder in the run_dir. This procedure is performed only if the SAVE
    folder is not already found in the run_dir.
    If the source_dir is not found an exception is raised.

    Args:
        source_dir (str) : name of the folder with the source nscf QuantumESPRESSO computation
        run_dir (st) : folder where the SAVE folder is linked or copied
        command (str) : command for generation of the SAVE Folder. Default is 'p2y -a 2'
        make_link (bool) : if True create a symbolic link

    """
    import os
    # check if the source_dir exists
    if not os.path.isdir(source_dir):
        raise ValueError('The source directory', source_dir,
                         ' does not exists.')
    # create run_dir
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)
        print('Create folder %s'%run_dir)
    # check if the SAVE folder is already present in the run_dir
    if os.path.isdir(os.path.join(run_dir,'SAVE')):
        print('SAVE folder already present in %s'%run_dir)
    else: # actions performed if the SAVE folder does not exists
        comm_str = 'cd %s; %s'%(source_dir,command)
        print('Executing command:', comm_str)
        os.system(comm_str)
        # copy (or create a symbolik link) of the SAVE folder in the run_dir
        src = os.path.abspath(os.path.join(source_dir,'SAVE'))
        dest = os.path.abspath(os.path.join(run_dir,'SAVE'))
        if make_link:
            os.symlink(src,dest,target_is_directory=True)
            print('Create a symlink of %s in %s'%(src,run_dir))
        else:
            from shutil import copytree
            copytree(src,dest)
            print('Create a copy of %s in %s'%(src,run_dir))
        # build the r_setup
        comm_str = 'cd %s;OMP_NUM_THREADS=1 yambo'%run_dir
        print('Executing command:', comm_str)
        os.system(comm_str)
        
