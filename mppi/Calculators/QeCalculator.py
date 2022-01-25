"""
This module manages calculations performed with QuantumESPRESSO.
The run of the computation is performed by the python subprocess package (direct scheduler)
or by the slurm scheduler.
"""

from .Runner import Runner
from mppi.Utilities import Tools
from mppi.Calculators.RunRules import build_slurm_header, mpi_command
import os

class QeCalculator(Runner):
    """
    Perform a QuantumESPRESSO calculation. The parameters used to define the parellelization
    strategy are defined in the runRules objet, provided in the constructor of this class

    Computations are managed by a scheduler that,
    in the actual implementation of the class, can be `direct` or `slurm`.

    Parameters:
       runRulues (:class:`RunRules`) : instance of the :class:`RunRules` class
       #omp (:py:class:`int`) : value of the OMP_NUM_THREADS variable
       #mpi (:py:class:`int`) : number of mpi processes
       #mpi_run (:py:class:`string`) : command for the execution of mpirun, e.g. 'mpirun -np' or 'mpiexec -np'
       executable (:py:class:`string`) : set the executable (pw.x, ph.x, ..) of the QuantumESPRESSO package
       #scheduler (:py:class:`string`) : choose the scheduler used to submit the job, actually the choices implemented are
            'direct', that runs the computation using the python subprocess package and 'slurm', that creates a slurm script
       skip (:py:class:`bool`) : if True evaluate if the computation can be skipped. This is done by checking if the log
            file of the run contains the string `job_done`, defined as a data member of this class
       clean_restart (:py:class:`bool`) : if True delete the folder $prefix.save before running the computation
       dry_run (:py:class:`bool`) : with this option enabled the calculator setup the calculations and write the script
            for submitting the job, but the computations are not run
       wait_end_run (:py:class:`bool`) : with this option disabled the run method does not wait the end of the run.
            This option may be useful for interacting with the code in particular in _asincronous_ computation managed
            by the slurm scheduler
       #sbatch_options (:py:class:`list`) : the elements of this list are strings used as options in the slurm script.
       #    For instance it is possible to specify the number of tasks per node as `--ntasks-per-node=16`
       activate_BeeOND (:py:class:`bool`) :  if True set I/O of the run in the BeeOND_dir created by the slurm scheduler.
            With this options enabled the ``out_dir`` of the run is set in the ``BeenOND_dir`` folder and the input wavefunction
            of the source folder (if needed) are copied in the ``BeeOND_dir``. At the end of the run the ``out_dir`` is moved
            in its original path. The value of the ``BeeOND_dir`` is written as a data member of the class and can be modified
            if needed
       verbose (:py:class:`bool`) : set the amount of information provided on terminal
       kwargs : other parameters that are stored in the _global_options dictionary

    Computations are performed in the folder specified by the ``run_dir`` parameter. The input and
    the log files are written in the run_dir. Instead, the $prefix.xml file and the $prefix.save
    folders are written in the ``out_dir`` path. The values of the prefix and out_dir variables
    are read from the input file. If the ``out_dir`` path is a relative path its root is located
    in the ``run_dir`` folder.

    Example:
     >>> rr = RunRules(scheduler='slurm',ntasks_per_node=4,memory='12GB')
     >>> code = calculator(rr,skip=True,clean_restart=True,verbose=True)
     >>> code.run(input = ..., run_dir = ...,name = ..., source_dir = ..., **kwargs)

     where the arguments of the run method are:

    Args:
        run_dir (:py:class:`string`) : the folder in which the simulation is performed
        input (:py:class:`string`) : instance of the :class:`PwInput` class
            that define the input object
        name (:py:class:`string`) : string with the name associated to the input file.
            Usually it is convenient to set the name equal to the prefix of the input object so
            the name of the input file and the prefix folder built by QuantumESPRESSO are the same
        source_dir (:py:class:`string`) : useful for the nscf computations. The source folder contains
            the wave-functions created by a scf calculation. If present the class copies this folder in the
            ``out_dir`` with the name $prefix.save
        kwargs : other parameters that are stored in the run_options dictionary

    """
    BeeOND_dir = '/mnt/${SLURM_JOB_USER}-jobid_${SLURM_JOB_ID}'
    job_done = 'JOB DONE.'

    def __init__(self,runRules,executable = 'pw.x', skip =  True, clean_restart = True,
                 dry_run = False, wait_end_run = True, activate_BeeOND = False,
                 verbose = True, **kwargs):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self,**runRules,executable=executable,skip=skip,clean_restart=clean_restart,
                        dry_run=dry_run,wait_end_run=wait_end_run,activate_BeeOND=activate_BeeOND,
                        verbose=verbose, **kwargs)
        print('Initialize a QuantumESPRESSO calculator')
    # def __init__(self,
    #              omp = os.environ.get('OMP_NUM_THREADS', 1), mpi = 2, mpi_run = 'mpirun -np',
    #              executable = 'pw.x', scheduler = 'direct', skip =  True, clean_restart = True,
    #              dry_run = False, wait_end_run = True, sbatch_options = [], activate_BeeOND = True,
    #               verbose = True, **kwargs):
    #     # Use the initialization from the Runner class (all options inside _global_options)
    #     Runner.__init__(self, omp=omp, mpi=mpi, mpi_run=mpi_run, executable=executable,
    #                     scheduler=scheduler,skip=skip, clean_restart=clean_restart,
    #                     dry_run=dry_run,wait_end_run=wait_end_run,sbatch_options=sbatch_options,
    #                     activate_BeeOND=activate_BeeOND,verbose=verbose, **kwargs)
    #     print('Initialize a QuantumESPRESSO calculator with scheduler %s'%self._global_options['scheduler'])

    def pre_processing(self):
        """
        Process local run dictionary to create the run directory and input file.
        If clean_restart is True the clean_run method is called before the run.
        Call the :meth:`copy_source_dir` that manages the source folder,
        if provided.

        """
        run_dir = self.run_options.get('run_dir', '.')
        input= self.run_options.get('input')
        name = self.run_options.get('name','default')
        clean_restart= self.run_options.get('clean_restart')
        verbose = self.run_options.get('verbose')
        self.is_to_run()
        is_to_run = self.run_options.get('is_to_run')

        # Create the run_dir and write the input file
        self._ensure_run_directory()
        if input is not None:
            input.write(os.path.join(run_dir,name)+'.in')
        else:
            print('input not provided')

        # Clean the folders before the run
        if is_to_run:
            if clean_restart:
                self.clean_run()
            else:
                if verbose: print('run performed starting from existing results')

        # Include the source .save folder (if provided) in the out_dir
        self.copy_source_dir()

        return {}

    def process_run(self):
        """
        Method associated to the running of the executable. The method runs the computation
        and wait the end of the computation before passing to the :meth:`post_processing` method.

        """
        is_to_run = self.run_options.get('is_to_run')

        if is_to_run:
            job = self.run_job()
            self.wait(job)
        return {}

    def post_processing(self):
        """
        Return the name, including the path, of the data-file-schema.xml file. If the file is absent the
        method displays a warning.

        Return:
            :py:class:`string` : name, including the path, of the xml data-file-schema file

        """
        input = self.run_options.get('input')
        prefix = input.get_prefix()
        out_dir = self._get_outdir_path()
        save_dir = os.path.join(out_dir,prefix)+'.save'
        result = os.path.join(save_dir,'data-file-schema.xml')
        if not os.path.isfile(result):
            print('Expected file %s not found'%result)
            print("""
                Check if wait_end_run is False or the dry_run option is active.
                Otherwise a possible error has occured during the computation""")
        return result

    def is_to_run(self):
        """
        The method evaluates if the computation has to be performed.  If ``skip`` is
        False the run is always performed, instead if ``skip`` is True the method
        checks if the log file exsists and contains the string `job_done`
        defined as a member of the class.
        The method adds the key `is_to_run` to ``the run_options`` of the class.

        """
        skip = self.run_options.get('skip')
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        logfile = os.path.join(run_dir,name)+'.log'
        verbose = self.run_options.get('verbose')

        if not skip:
            self.run_options['is_to_run'] = True
        else:
            if Tools.find_string_file(logfile,self.job_done) is not None:
                if verbose: print('Skip the run of',name)
                self.run_options['is_to_run'] = False
            else:
                self.run_options['is_to_run'] = True

    def run_job(self):
        """
        Run the computation. The operations performed depend on the scheduler adopted.
        If the dry_run option is enabled the run is not performed but the slurm script
        is written on disk.

        Return:
            The type of the object depends on the chosen scheduler. For scheduler `direct`
            job is an instance of Popen, while for `slurm` scheduler job is the name of the
            slurm script.

        """
        from subprocess import Popen
        run_dir = self.run_options.get('run_dir', '.')
        scheduler = self.run_options['scheduler']
        dry_run = self.run_options.get('dry_run')
        verbose = self.run_options.get('verbose')

        if scheduler == 'direct':
            # Set the OMP_NUM_THREADS variable in the environment
            os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])
            if not dry_run:
                comm_str = 'cd %s ; %s'%(run_dir,self.run_command())
                job = Popen(comm_str, shell = True)
            else:
                job = None
                if verbose: print('Dry_run option active. Computation not performed')
        elif scheduler == 'slurm':
            job = self.build_slurm_script()
            if not dry_run:
                slurm_submit = 'cd %s ; sbatch %s.sh' %(run_dir,job)
                if verbose: print('slurm submit: ',slurm_submit )
                slurm_run = Popen(slurm_submit, shell = True)
            else:
                if verbose: print('Dry_run option active. Script not submitted')
        else:
            print('scheduler unknown')
        return job

    def wait(self,job):
        """
        Wait the end of the job. If the dry_run option is enabled or wait_end_run is False
        the check is not performed.

        Args:
            jobs : The reference to the job to be executed. If the scheduler is `direct`
                jobs is an instance of Popen of the :py:class:subprocess package. If the
                scheduler is `slurm` jobs is a string with the name of the slurm script

        """
        import time
        dry_run = self.run_options.get('dry_run')
        wait_end_run = self.run_options.get('wait_end_run')
        name = self.run_options.get('name','default')
        verbose = self.run_options.get('verbose')
        delay = 1 # in second

        if wait_end_run is True and dry_run is False:
            message_written = False
            while not self.run_ended(job):
                if not message_written:
                    if verbose: print('computation %s is running...'%name)
                    message_written = True
                time.sleep(delay)
            if verbose: print('computation %s ended'%name)
        else:
            if verbose: print('The wait_end_run is False or the dry_run option is active. The calculator proceedes to the postprocessing')

    def build_slurm_script(self):
        """
        Create the slurm script associated to the run.

        Return:
            :py:class:`string`: string with the name of the slurm script

        """
        #omp = self.run_options.get('omp')
        #mpi = self.run_options.get('mpi')
        input = self.run_options.get('input')
        prefix = input.get_prefix()
        out_dir = input.get_outdir()
        name = self.run_options.get('name','default')
        run_dir = self.run_options.get('run_dir', '.')
        out_dir_path = self._get_outdir_path()
        save_dir = os.path.join(out_dir_path,prefix)+'.save'
        job = 'job_'+name

        #sbatch_options = self.run_options.get('sbatch_options')
        activate_BeeOND = self.run_options.get('activate_BeeOND')
        comm_str = self.run_command()

        lines = build_slurm_header(self.run_options)

        # lines = []
        # lines.append('#!/bin/bash')
        # lines.append('#SBATCH --ntasks=%s           ### Number of tasks (MPI processes)'%mpi)
        # lines.append('#SBATCH --cpus-per-task=%s    ### Number of threads per task (OMP threads)'%omp)
        # for option in sbatch_options: # add other SBATCH options
        #     lines.append('#SBATCH %s'%option)
        # lines.append('#SBATCH --output=%s.out'%job)
        # lines.append('')
        #
        # lines.append('export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK')
        # lines.append('export OUT_DIR=%s'%out_dir)
        # lines.append('export OUT_DIR_PATH=%s'%out_dir_path)
        # lines.append('export SAVE_DIR=%s'%save_dir)
        # lines.append('')
        #
        # lines.append('echo "Cluster name $SLURM_CLUSTER_NAME"')
        # lines.append('echo "Job name $SLURM_JOB_NAME "')
        # lines.append('echo "Job id $SLURM_JOB_ID"')
        # lines.append('echo "Job nodelist $SLURM_JOB_NODELIST"')
        # lines.append('echo "Number of nodes $SLURM_JOB_NUM_NODES"')
        # lines.append('echo "Number of mpi $SLURM_NTASKS"')
        # lines.append('echo "Number of threads per task $SLURM_CPUS_PER_TASK"')
        # lines.append('echo "OUT_DIR input parameter is $OUT_DIR"')
        # lines.append('echo "OUT_DIR path is $OUT_DIR_PATH"')
        # lines.append('echo "SAVE_DIR path is $SAVE_DIR"')
        # lines.append('echo " "')
        # lines.append('')

        #lines.append('#SBATCH --job-name=%s'%job)
        #lines.append('#SBATCH --output=%s.out'%job)

        lines.append('export OUT_DIR=%s'%out_dir)
        lines.append('export OUT_DIR_PATH=%s'%out_dir_path)
        lines.append('export SAVE_DIR=%s'%save_dir)
        lines.append('echo "OUT_DIR input parameter is $OUT_DIR"')
        lines.append('echo "OUT_DIR path is $OUT_DIR_PATH"')
        lines.append('echo "SAVE_DIR path is $SAVE_DIR"')
        lines.append('echo " "')
        lines.append('')

        if activate_BeeOND:
            lines.append('export BEEOND_DIR=%s'%self.BeeOND_dir)
            lines.append('echo "THe BeeOND option is activated. The I/O is performed in $BEEOND_DIR"')
            lines.append('echo "BEEOND_DIR path is $BEEOND_DIR"')
            lines.append('if [ ! -d $BEEOND_DIR ]; then')
            lines.append('echo "$BEEOND_DIR not found!"')
            lines.append('exit')
            lines.append('fi')
            lines.append('echo " "')
            lines.append('')
            lines.append('echo "Change the outdir key of the input from $OUT_DIR to $BEEOND_DIR"')
            lines.append('sed -i "/outdir/s:%s:%s:" %s.in'%(out_dir,self.BeeOND_dir,name))
            lines.append('echo " "')
            lines.append('')
            # If there is a save_dir (created by the pre_processing method) it is copied
            # in the BeeOND_dir
            if os.path.isdir(save_dir):
                lines.append('echo "found SAVE_DIR folder $SAVE_DIR. Copy the SAVE_DIR in the $BEEOND_DIR folder"')
                lines.append('echo "rsync -azv $SAVE_DIR $BEEOND_DIR"')
                lines.append('rsync -azv $SAVE_DIR $BEEOND_DIR')
                lines.append('echo " "')
                lines.append('')

        lines.append('echo "execute : %s"'%comm_str)
        lines.append(comm_str)
        lines.append('echo " "')
        lines.append('')

        if activate_BeeOND:
            lines.append('echo "Change the outdir key of the input to its original value $OUT_DIR"')
            lines.append('sed -i "/outdir/s:%s:%s:" %s.in'%(self.BeeOND_dir,out_dir,name))
            lines.append('echo "rsync -azv $BEEOND_DIR/ $OUT_DIR_PATH"')
            lines.append('rsync -azv $BEEOND_DIR/ $OUT_DIR_PATH')
            lines.append('echo " "')
            lines.append('')

        lines.append('echo "JOB_DONE"')
        f = open(os.path.join(run_dir,job+'.sh'),'w')
        f.write('\n'.join(lines))
        f.close()

        return job

    def run_command(self):
        """
        Define the run command used to run the computation.

        Return:
            :py:class:`string` : command that runs the computation

        """
        executable = self.run_options.get('executable')
        #mpi = self.run_options.get('mpi')
        #mpi_run = self.run_options.get('mpi_run')
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        verbose = self.run_options.get('verbose')

        mpi_run = mpi_command(self.run_options)
        command = mpi_run + ' ' + executable
        #command = mpi_run + ' ' + str(mpi) + ' ' + executable
        input_name = name + '.in'
        output_name = name + '.log'
        comm_str =  command + ' -inp %s > %s'%(input_name,output_name)
        if verbose: print('run command: %s' %comm_str)

        return comm_str

    def run_ended(self,job):
        """
        Check the status of the running job.

        Args:
            job : reference to the actual job. job is an istance of Popen for `direct` scheduler
                or a string for `slurm` scheduler

        Return:
            :py:class:`bool`: return True if the computation is ended and False if it is running

        """
        scheduler = self.run_options.get('scheduler')
        run_dir = self.run_options.get('run_dir', '.')

        if scheduler == 'direct':
            if job.poll() is not None: is_ended = True
            else: is_ended = False
        if scheduler == 'slurm':
            job_out = os.path.join(run_dir,job+'.out')
            if not os.path.isfile(job_out):
                is_ended = False
            else:
                with open(job_out, 'r') as f:
                    last_line = f.read().splitlines()[-1]
                    if last_line == 'JOB_DONE': is_ended = True
                    else: is_ended = False

        return is_ended

    def clean_run(self):
        """
        Clean the run before performing the computation. Delete the $name.log and
        the job_$name.out file, located in the `run_dir`, and the $prefix.xml file
        and the $prefix.save folder located in the `out_dir`. Finally, if the
        `out_dir` is empty it is deleted.

        """
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        input = self.run_options.get('input')
        verbose = self.run_options.get('verbose')
        prefix = input.get_prefix()
        out_dir = self._get_outdir_path()

        logfile = os.path.join(run_dir,name)+'.log'
        job_out = os.path.join(run_dir,'job_'+name+'.out')
        xmlfile = os.path.join(out_dir,prefix)+'.xml'
        save_dir = os.path.join(out_dir,prefix)+'.save'

        if os.path.isfile(logfile):
            if verbose: print('delete log file:',logfile)
            os.system('rm %s'%logfile)
        if os.path.isfile(job_out):
            if verbose: print('delete job_out script:',job_out)
            os.system('rm %s'%job_out)
        if os.path.isfile(xmlfile):
            if verbose: print('delete xml file:',xmlfile)
            os.system('rm %s'%xmlfile)
        if os.path.isdir(save_dir):
            if verbose: print('delete folder:',save_dir)
            os.system('rm -r %s'%save_dir)
        # delete the out_dir (if it is empty)
        if os.path.isdir(out_dir) and not os.listdir(out_dir):
            if verbose: print('delete the out_dir:',out_dir)
            os.system('rm -r %s'%out_dir)

    def _ensure_run_directory(self):
        """
        Create the run_dir, if it does not exists

        """
        run_dir = self.run_options.get('run_dir','.')
        verbose = self.run_options.get('verbose')

        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
            if verbose: print("create the run_dir folder : '%s'" %run_dir)

    def copy_source_dir(self):
        """
        Copy the source_dir (if provided) in the out_dir and atttibute to the copied folder
        the name $prefix.save. The operation is performed only if a folder with the target name
        of the source_dir is not found in the out_dir

        Args:
            source_dir: the name of the source_dir (tipically it is the .save folder
            of the scf calculation that contains the wave-functions of the ground state).

        """
        from shutil import copytree
        source_dir = self.run_options.get('source_dir',None)
        input = self.run_options.get('input')
        prefix = input.get_prefix()
        out_dir = self._get_outdir_path()
        verbose = self.run_options.get('verbose')

        if source_dir is not None:
            dest_dir = os.path.join(out_dir,prefix)+'.save'
            if not os.path.isdir(dest_dir):
                if verbose: print('copy source_dir %s in the %s'%(source_dir,dest_dir))
                copytree(source_dir,dest_dir)
            else:
                if verbose:
                    print('The folder %s already exists. Source_dir % s not copied'
                    %(dest_dir,source_dir))

    def _get_outdir_path(self):
        """
        Get the absolute out_dir path. The path is built using the ``outdir`` parameter
        of the input file. If ``outdir`` is provided as a relative address the path starts
        from the ``run_dir`` of the calculator.

        """
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options.get('input')
        out_dir = input.get_outdir()
        if os.path.isabs(out_dir):
            return out_dir
        else:
            return os.path.abspath(os.path.join(run_dir,out_dir))
