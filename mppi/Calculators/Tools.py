"""
This module contains some functions that implement useful operations to deal
with the Yambo and QuantumESPRESSO calculators.
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

def make_p2y(yambo_dir = './', input_dir ='./', overwrite_if_found = False, p2y_command = 'p2y'):
    """
    Run p2y in the ``yambo_dir`` using the PW wave functions. The function creates the
    SAVE Folder using the command provided as the ``p2y_command`` parameter.

        cd $yambo_dir; $p2y_command -I $input_dir

    If the ``yambo_dir`` does not exists it is built by the function.
    If a SAVE folder is already found in the ``yambo_dir`` no operations are performed, unless the
    overwrite_if_found option is `True`.

    Args:
        yambo_dir (:py:class:`string`) : location of the yambo_dir, it can be a nested directory path and if the path
            is not found is created by the function using the :py:meth:`os.makedirs`
        input_dir (:py:class:`string`) : name of the folder with the PW wave functions
        overwrite_if_found (:py:class:`bool`) : if True delete the SAVE folder
        p2y_command (:py:class:`string`) : command for execution of the p2y program. Default is 'p2y'

    """
    if not os.path.isdir(yambo_dir):
        os.makedirs(yambo_dir)
        print('Create the folder path %s'%yambo_dir)
    if not os.path.isdir(input_dir):
        raise ValueError('The input directory', input_dir,
                         ' does not exists.')
    # Evaluate if the SAVE folder has to be removed if found
    save_dir = os.path.join(yambo_dir,'SAVE')
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
        comm_str = 'cd %s; %s'%(yambo_dir,p2y_command)
        comm_str += ' -I %s'%os.path.relpath(input_dir,start=yambo_dir)
        print('Executing command:', comm_str)
        os.system(comm_str)

def build_r_setup(yambo_dir = './', overwrite_if_found = False, yambo_command = 'yambo'):
    """
    Create the `r_setup` file by executing the `yambo_command` in the ``yambo_dir``. If an
    instance of the ``r_setup`` file is found do not perform any operation unless the
    `overwrite_if_found` option is `True`.

    Args:
        yambo_dir (:py:class:`string`) : location of the yambo_dir
        overwrite_if_found (:py:class:`bool`) : if True clean the yambo_dir
        yambo_command (:py:class:`string`) : command for generation the r_setup file. Default is 'yambo'

    """
    r_setup_file = os.path.join(yambo_dir,'r_setup')
    l_setup_file = os.path.join(yambo_dir,'l_setup')
    yambo_in_file =  os.path.join(yambo_dir,'yambo.in')
    if overwrite_if_found:
        print('Clean the yambo_dir %s to build a new r_setup file'%yambo_dir)
        for f in [r_setup_file,l_setup_file,yambo_in_file]:
            if os.path.isfile(f) :
                comm_str = 'rm %s'%f
                print('Executing command:', comm_str)
                os.system(comm_str)
    if not os.path.isfile(r_setup_file):
        print('Build the r_setup in the yambo_dir path %s'%yambo_dir)
        comm_str = 'cd %s; %s'%(yambo_dir,yambo_command)
        os.system(comm_str)

def init_yambo_dir(yambo_dir = '.', input_dir = './', overwrite_if_found = False,
        p2y_command = 'p2y', yambo_command = 'yambo') :
    """
    Create and initialize the run directory where Yambo computations can be performed. The function runs the
    `make_p2y` and `build_r_setup` function of the module.

    Args:
        yambo_dir (:py:class:`string`) : location of the yambo_dir, it can be a nested directory path and if the path
            is not found is created by the function using the :py:meth:`os.makedirs`
        input_dir (:py:class:`string`) : name of the folder with the PW wave functions
        make_link (:py:class:`bool`) : if True create a symbolic link of the SAVE folder, otherwise the SAVE folder is copied in
            the run_dir
        overwrite_if_found (:py:class:`bool`) : if True delete the SAVE folder in the run_dir and the r_setup and l_setup (if found)
            and build them again (also the 'yambo.in' input file is deleted to build the `r_setup` file)
        p2y_command (:py:class:`string`) : command for execution of the p2y program. Default is 'p2y'
        yambo_command (:py:class:`string`) : command for generation the r_setup file. Default is 'yambo'

    """
    make_p2y(yambo_dir=yambo_dir,input_dir=input_dir,overwrite_if_found=overwrite_if_found,p2y_command=p2y_command)
    build_r_setup(yambo_dir=yambo_dir,overwrite_if_found=overwrite_if_found,yambo_command=yambo_command)

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

def build_pw_klist(kpoints,kweight=None):
    """
    Build a list of kpoints with the structure [[k1_x,k1_y,k1_z,w1],...[kn_x,kn_y,kn_z,wn]], where wi is the
    weight of the i-th kpoint.

    Args:
        kpoints (:py:class:`array`) : array with the coordinates of the kpoints
        kweight (:py:class:`list`) : array with the weigth of each kpoint. If is None
            a uniform weight equal to 1 is attributed to each kpoint

    Returns:
        :py:class:`list` : list of kpoints with the structure [[k1_x,k1_y,k1_z,w1],...[kn_x,kn_y,kn_z,wn]]

    """
    if kweight is None: kweight = [1 for ind in range(len(kpoints))]
    klist = []
    for ind,k in enumerate(kpoints):
        klist.append(k.tolist()+[kweight[ind]])
    return klist
