"""
This module contains some functions that implement some operations useful to deal
with Yambo and QuantumESPRESSO calculators.
The module can be loaded in the notebook in one of the following way

>>> from mppi.Utilities import Tools

>>> Tools.init_yambo_run_dir

or as

>>> import mppi.Utilities.Tools

>>> Tools.init_yambo_run_dir

or to load directly some elements

>>> from mppi.Utilities.Tools import init_yambo_run_dir

>>> init_yambo_run_dir

"""
import numpy as np
import os

def find_string_file(file,string):
    """
    Look for a string in the lines of a file.

    Args:
        file (:py:class:`string`) : name of the file, including the path
        string (:py:class:`string`) : name of the string

    Return:
        :py:class:`string` : return the first occurence of the line that match
        the search. If no line is found or the file does not exsists return None

    """
    line = None
    if os.path.isfile(file):
        for l in open(file,'r'):
            if string in l:
                line = l
                break
    return line

def build_r_setup(run_dir, overwrite_if_found = False, yambo_command = 'yambo'):
    """
    Create the `r_setup` file by executing yambo with no arguments in the ``run_dir``.

    Args:
        run_dir (:py:class:`string`) : folder with the SAVE directory
        overwrite_if_found (:py:class:`bool`) : if True delete previous istance of the r_setup file (if found)
        yambo_command (:py:class:`string`) : command for generation the r_setup file. Default is 'yambo'

    """
    rsetup_file = os.path.join(run_dir,'r_setup')
    if os.path.isfile(rsetup_file) and overwrite_if_found:
        print('Remove previous instance of the r_setup file')
        comm_str = 'rm %s'%rsetup_file
        os.system(comm_str)
    if not os.path.isfile(rsetup_file):
        print('Build the r_setup in the run_dir path %s'%run_dir)
        comm_str = 'cd %s; %s'%(run_dir,yambo_command)
        os.system(comm_str)

def init_yambo_run_dir(source_dir, run_dir ='.', make_link = True, overwrite_if_found = False, yambo_command = 'yambo') :
    """
    Create and initialize the run directory where Yambo computations can be performed. The function creates the run_dir (if it
    does not exists), then perform a copy (or a link) of the source_dir into the run_dir and run the ``build_r_setup`` function.

    Args:
        source_dir (:py:class:`string`) : location of the SAVE folder with the Yambo core databases
        run_dir (:py:class:`string`) : location of the run_dir. The run_dir can be
            a nested directory path and if the path is not found is created by the function using the
            :py:meth:`os.makedirs`
        make_link (:py:class:`bool`) : if True create a symbolic link of the SAVE folder, otherwise the SAVE folder is copied in
            the run_dir
        overwrite_if_found (:py:class:`bool`) : if True delete the SAVE folder in the run_dir and the r_setup and l_setup (if found)
            and build them again (also the 'yambo.in' input file is deleted to build the `r_setup` file).
            Note that the SAVE folder in the run_dir is erased using the command ``rm -r run_dir/SAVE``, without the final `/`.
            This ensures that, if the SAVE in the run_dir is a symbolic link, the source folder is not erased
        yambo_command (:py:class:`string`) : command for generation the r_setup file. Default is 'yambo'

    """
    if not os.path.isdir(source_dir):
        raise ValueError('The SAVE directory', source_dir,
                         ' does not exists.')
    if not os.path.isdir(run_dir):
        os.makedirs(run_dir)
        print('Create folder path %s'%run_dir)
    save_dir = os.path.join(run_dir,'SAVE')
    # Evaluate if the save_dir folder has to be removed if found
    if os.path.isdir(save_dir):
        if overwrite_if_found:
            print('clean the run_dir %s to build a new SAVE folder'%run_dir)
            rl_setup_file = os.path.join(run_dir,'*_setup')
            yambo_in_file =  os.path.join(run_dir,'yambo.in')
            comm_str = 'rm -r %s %s %s'%(save_dir,rl_setup_file,yambo_in_file)
            print('Executing command:', comm_str)
            os.system(comm_str)
        else:
            print('SAVE folder already present in %s. No operations performed.'%run_dir)
    # Actions performed if the save_dir is not present (or if it has been removed)
    if not os.path.isdir(save_dir):
        src = os.path.abspath(source_dir)
        dest = os.path.abspath(save_dir)
        if make_link:
            os.symlink(src,dest,target_is_directory=True)
            print('Create a symlink of %s in %s'%(source_dir,run_dir))
        else:
            from shutil import copytree
            copytree(src,dest)
            print('Create a copy of %s in %s'%(src,run_dir))
    # Build the r_setup file
    build_r_setup(run_dir,overwrite_if_found=overwrite_if_found,yambo_command=yambo_command)

def make_p2y(source_dir, p2y_command = 'p2y', overwrite_if_found = False):
    """
    Run p2y the build the SAVE folder with the yambo core databases. The function creates the
    SAVE Folder in the source_dir using the command provided as the ``p2y_command`` parameter.
    If a SAVE folder is already found in the source_dir no operations are performed, unless the
    overwrite_if_found option is `True`.

    Args:
        source_dir (:py:class:`string`) : name of the folder with the source nscf QuantumESPRESSO computation
        p2y_command (:py:class:`string`) : command for generation of the SAVE Folder. Default is 'p2y'
        overwrite_if_found (:py:class:`bool`) : if True delete the SAVE folder

    Returns:
        :py:class:`list` : path of the SAVE folder

    """
    save_dir = os.path.join(source_dir,'SAVE')
    if not os.path.isdir(source_dir): # check if the source_dir exists
        raise ValueError('The source directory', source_dir,
                         ' does not exists.')
    # Evaluate if the save_dir folder has to be removed if found
    if os.path.isdir(save_dir):
        if overwrite_if_found:
            print('Delete the SAVE folder %s'%save_dir)
            comm_str = 'rm -r %s'%save_dir
            print('Executing command:', comm_str)
            os.system(comm_str)
        else:
            print('SAVE folder %s already present. No operations performed.'%save_dir)
    # Actions performed if the save_dir is not present (or if it has been removed)
    if not os.path.isdir(save_dir):
        comm_str = 'cd %s; %s'%(source_dir,p2y_command)
        print('Executing command:', comm_str)
        os.system(comm_str)
    return save_dir

def build_FixSymm_input(run_dir, polarization= 'linear', Efield1 = [1.,0.,0.], Efield2 = [0.,1.,0.],
                removeTimeReversal = True):
    """
    Create the input file for the fixSymm procedure to remove the symmetries broken by the electric field.

    Args:
        run_dir (:py:class:`string`) : folder with the SAVE directory
        polarization (:py:class:`string`) : specifies the linear or circular polarization of the field
        Efield1 (:py:class:`list`) : direction of the first electric field
        Efield2 (:py:class:`list`) : direction of the second electric field. Useful for the circular polarization case
        removeTimeReversal (:py:class:`bool`) : if True remove the time reversal symmetry

    """
    from mppi import InputFiles as I
    l_fixsyms_file = os.path.join(run_dir,'l_fixsyms*')
    comm_str = 'rm %s'%l_fixsyms_file
    os.system(comm_str)

    fixSymm_inp = I.YamboInput('ypp -y -V all',folder=run_dir,filename='FixSymm.in')
    if removeTimeReversal: fixSymm_inp.removeTimeReversal()
    if polarization == 'circular':
        fixSymm_inp.set_ypp_extFields(Efield1=Efield1,Efield2=Efield2)
    elif polarization == 'linear':
        fixSymm_inp.set_ypp_extFields(Efield1=Efield1)
    else:
        print('Specify a correct polarization for the field')
    return fixSymm_inp

def build_pw_kpath(*kpoints,numstep=40):
    """
    Build a list of kpoints to be passed to the set_kpoints methods of the :class:`PwInput`
    for computing the band structure along a path.

    Example:
        >>> build_kpath(L,G,X,K,G,numstep=30)

    Args:
        kpoints : arguments that specify the high symmetry points along the k-path
        numstep (int): specifies the number of intermediate points used to build the path

    Returns:
        :py:class:`list` : list of kpoints as nedded by pw in the bands computation with the tpiba_b option

    """
    klist = []
    for k in kpoints[:-1]:
        klist.append(k+[numstep])
    klist.append(kpoints[-1]+[0])
    return klist
