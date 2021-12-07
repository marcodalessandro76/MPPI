"""
This module contains some functions that implement common operations useful to deal
with Yambo and QuantumESPRESSO and to perform data analysis.
The module can be loaded in the notebook in one of the following way

>>> from mppi.Utilities import Tools

>>> Tools.eval_FT

or as

>>> import mppi.Utilities.Tools

>>> Tools.eval_FT

or to load directly some elements

>>> from mppi.Utilities.Tools import eval_FT

>>> eval_FT

"""
import numpy as np
import os

def eval_FT(time,values,time_units = 'fs', verbose = True):
    """
    Compute the FT of the (real-valued) values array using the time array as x-values.
    The FT is provided only for positive energies since due to the real-valued of the
    values array the FT is assumed to be an even function.

    Args:
        time (:py:class:`np.array`) : x-values array (assumes a uniform sampling)
        values (:py:class:`np.array`) : values array (assumes real-value function)
        time_units (:py:class:`string`) : set the units to compute the energy array.
            Default is 'fs' and the other possible choice is 'ps'
        verbose (:py:class:`bool`) : defines the amount of information provided
            on terminal

    Return:
        :py:class:`dict`  : dictionary with the energy array (in eV), the real and imaginary part
            of the FT and the associated module

    """
    from mppi.Utilities import Constants as C
    h = C.Planck_ev_ps
    if time_units == 'fs':
        h *= 1e3
        if verbose : print('Time values expressed in fs. Energy array is provided in eV.')
    elif time_units == 'ps':
        if verbose : print('Time values expressed in ps. Energy array is provided in eV.')
    else:
        print('Unkwon time units. FT has not been computed.')
        return 0
    dt = time[1]-time[0]
    N = len(time)
    freqs = np.fft.fftfreq(N,d=dt)
    energy = h*freqs[0:int(N/2)]
    ft = np.fft.fft(values)[0:int(N/2)]
    ft_abs = np.sqrt(ft.real**2+ft.imag**2)
    if verbose :
        print('Maximum energy value =',energy[-1])
        print('Energy sampling step =',energy[1]-energy[0])
    return {'energy':energy,'ft_real':ft.real,'ft_im':ft.imag,'ft_abs':ft_abs}


def build_kpath(*kpoints,numstep=40):
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

def build_SAVE(source_dir, run_dir, command = 'p2y -a 2', make_link = True, overwrite_if_found = False):
    """
    Build the SAVE folder for a yambo computation.

    The function creates the SAVE folder in the source_dir using the command provided
    as the command parameter (the option -a 2 ensures that labelling of the
    high-symmetry kpoints is consistent in both QE and Yambo) and create a symbolic
    link (or a copy) of the SAVE folder in the run_dir. This procedure is performed only if the SAVE
    folder is not already found in the run_dir, unless the ``overwrite_if_found`` parameter is True.
    If the source_dir is not found an exception is raised.

    Args:
        source_dir (:py:class:`string`) : name of the folder with the source nscf QuantumESPRESSO computation
        run_dir (:py:class:`string`) : folder where the SAVE folder is linked or copied. The run_dir can be
            a nested directory path and if the path is not found is created by the function using the
            :py:meth:`os.makedirs`
        command (:py:class:`string`) : command for generation of the SAVE Folder. Default is 'p2y -a 2'
        make_link (:py:class:`bool`) : if True create a symbolic link
        overwrite_if_found (:py:class:`bool`) : if True delete the SAVE folder in the run_dir and the
            r_setup (if found) and build them again

    """
    SAVE_dir = os.path.join(run_dir,'SAVE')
    if not os.path.isdir(source_dir): # check if the source_dir exists
        raise ValueError('The source directory', source_dir,
                         ' does not exists.')
    if not os.path.isdir(run_dir):
        os.mkdirs(run_dir)
        print('Create folder %s'%run_dir)
    # Evaluate if the SAVE_dir folder has to be removed if found
    if os.path.isdir(SAVE_dir):
        if overwrite_if_found:
            print('clean the run_dir %s to build a new SAVE folder'%run_dir)
            comm_str = 'rm -r %s'%SAVE_dir
            print('Executing command:', comm_str)
            os.system(comm_str)
            r_setup_files = os.path.join(run_dir,'r_setup')
            comm_str = 'rm %s'%r_setup_files
            print('Executing command:', comm_str)
            os.system(comm_str)
        else:
            print('SAVE folder already present in %s. No operations performed.'%run_dir)
    # Actions performed if the SAVE_dir is not present (or if it has been removed)
    if not os.path.isdir(SAVE_dir):
        comm_str = 'cd %s; %s'%(source_dir,command)
        print('Executing command:', comm_str)
        os.system(comm_str)
        src = os.path.abspath(os.path.join(source_dir,'SAVE'))
        dest = os.path.abspath(os.path.join(run_dir,'SAVE'))
        if make_link: # copy (or create a symbolik link) of the SAVE folder in the run_dir
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

def make_FixSymm(run_dir, polarization= 'linear', Efield1 = [1.,0.,0.], Efield2 = [0.,1.,0.],
                removeTimeReversal = True, overwrite_if_found = False):
    """
    Perform the fixSymm procedure to remove the symmetries broken by the electric field.
    The procedure creates the FixSymm folder into run_dir and run yambo_rt into the FixSymm to generate the r_setup.
    If a SAVE folder is already present in the run_dir/FixSymm path no operations are performed,
    unless the ``overwrite_if_found`` parameter is True.

    Args:
        run_dir (:py:class:`string`) : folder with the SAVE directory
        polarization (:py:class:`string`) : specifies the linear or circular polarization of the field
        Efield1 (:py:class:`list`) : direction of the first electric field
        Efield2 (:py:class:`list`) : direction of the second electric field. Useful for the circular polarization case
        removeTimeReversal (:py:class:`bool`) : if True remove the time reversal symmetry
        overwrite_if_found (:py:class:`bool`) : if True delete the SAVE folder in the run_dir/FixSymm and the
            r_setup (if found) and build them again.

    Note:
        Although the function does not remove the content of the FixSymm  folder, when 'ypp -y' is executed this folder
        is erased. This fact must be considered if there are relevant data in the FixSymm


    """
    from mppi import InputFiles as I, Calculators as C
    fixsymm_dir = os.path.join(run_dir,'FixSymm')
    SAVE_dir = os.path.join(fixsymm_dir,'SAVE')
    # Evaluate if the SAVE_dir folder has to be removed if found
    if os.path.isdir(SAVE_dir):
        if overwrite_if_found:
            print('clean the FixSymm folder %s to build a new SAVE folder'%fixsymm_dir)
            comm_str = 'rm -r %s'%SAVE_dir
            print('Executing command:', comm_str)
            os.system(comm_str)
            l_fixsyms_file = os.path.join(run_dir,'l_fixsyms')
            comm_str = 'rm %s'%l_fixsyms_file
            print('Executing command:', comm_str)
            os.system(comm_str)
            r_setup_file = os.path.join(fixsymm_dir,'r_setup')
            comm_str = 'rm %s'%r_setup_file
            print('Executing command:', comm_str)
            os.system(comm_str)
            l_Fixsymm_file = os.path.join(fixsymm_dir,'l-FixSymm_fixsyms')
            comm_str = 'rm %s'%l_Fixsymm_file
            print('Executing command:', comm_str)
            os.system(comm_str)
            r_Fixsymm_file = os.path.join(fixsymm_dir,'r-FixSymm_fixsyms')
            comm_str = 'rm %s'%r_Fixsymm_file
            print('Executing command:', comm_str)
            os.system(comm_str)
        else:
            print('SAVE folder already present in %s. No operations performed.'%fixsymm_dir)
    if not os.path.isdir(SAVE_dir):
        print('Perform the fixSymm in the folder %s'%run_dir)
        fixSymm_inp = I.YamboInput('ypp -y',folder=run_dir,filename='FixSymm.in')
        if removeTimeReversal:
            fixSymm_inp.removeTimeReversal()
        if polarization == 'circular':
            fixSymm_inp.set_ypp_extFields(Efield1=Efield1,Efield2=Efield2)
        elif polarization == 'linear':
            fixSymm_inp.set_ypp_extFields(Efield1=Efield1,Efield2=[0.,0.,0.])
        else:
            print('Specify a correct polarization for the field')
        code = C.YamboCalculator(omp=1,mpi=1,executable='ypp',skip=False,verbose=False)
        code.run(input=fixSymm_inp,name='FixSymm',run_dir=run_dir)
        # build the real-time r_setup
        command = 'cd %s; OMP_NUM_THREADS=1 yambo_rt'%fixsymm_dir
        os.system(command)

def get_variable_from_db(ndb_file,var_name):
    """
    Extract the value of a variable from a ndb database

    Args:
        ndb_file (:py:class:`string`) : the name of the database
        var_name (:py:class:`string`) : name of the variable

    Return:
        :py:class:`numpy.ndarray`  : array with the values of the variable
    """
    from netCDF4 import Dataset as Ds
    db = Ds(ndb_file)
    var = np.array(db.variables[var_name])
    return var
