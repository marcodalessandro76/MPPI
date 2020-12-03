"""
This module manages calculations performed with QuantumESPRESSO.
Actually the run of the computation can be managed by the python subprocess package (direct scheduler)
or by the slurm scheduler.
"""

from .Runner import Runner
import os

class QeCalculator(Runner):
    """
    Perform a QuantumESPRESSO calculation. Computations are managed by a scheduler that,
    in the actual implementation of the class, can be `direct` or `slurm`.

    Parameters:
       omp (:py:class:`int`) : value of the OMP_NUM_THREADS variable
       mpi (:py:class:`int`) : number of mpi processes
       mpi_run (:py:class:`string`) : command for the execution of mpirun, e.g. 'mpirun -np' or 'mpiexec -np'
       executable (:py:class:`string`) : set the executable (pw.x, ph.x, ..) of the QuantumESPRESSO package
       scheduler (:py:class:`string`) : choose the scheduler used to submit the job, actually the choices implemented are
            'direct' that runs the computation using the python subprocess package and 'slurm' that creates a slurm script
       skip (:py:class:`bool`) : if True evaluate if the computation can be skipped. This is done by checking if the file
            $prefix.xml is present in the run_dir folder
       clean_restart (:py:class:`bool`) : if True delete the folder $prefix.save before running the computation
       verbose (:py:class:`bool`) : set the amount of information provided on terminal
       kwargs : other parameters that are stored in the _global_options dictionary

    Example:
     >>> code = calculator(omp=1,mpi=4,mpi_run='mpirun -np',skip=True,clean_restart=True,verbose=True,scheduler='direct')
     >>> code.run(input = ..., run_dir = ...,name = ..., source_dir = ..., **kwargs)

     where the arguments of the run method are:

    Args:
        run_dir (:py:class:`string`) : the folder in which the simulation is performed
        input (:py:class:`string`) : instance of the :class:`PwInput` class
            that define the input object
        name (:py:class:`string`) : string with the name associated to the input file.
            Usually it is convenient to set the name equal to the prefix of the input object so
            the name of the input file and the prefix folder built by QuantumESPRESSO are the same
        source_dir (:py:class:`string`) : location of the scf source folder for a nscf computation.
            If present the class copies this folder in the run_dir with the name $prefix.save
        kwargs : other parameters that are stored in the run_options dictionary

    The calculator looks for the following variables in the run_options dictionary. These options
    may be useful for interacting with the code in particular in _asincronous_ computation managed
    the slurm scheduler.

        `dry_run=True` with this option the calculator setup the calculations and write the script
        for submitting the job, but the computations are not run

        `wait_end_run=False` with this option the wait of the end of the run is suppressed

        `sbatch_options = [option1,option2,....]` allows the user to include further options in the slurm script

    """

    def __init__(self,
                 omp = os.environ.get('OMP_NUM_THREADS', 1), mpi = 2, mpi_run = 'mpirun -np',
                 executable = 'pw.x', scheduler = 'direct', skip =  True, clean_restart = True,
                 verbose = True, **kwargs):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=omp, mpi=mpi, mpi_run=mpi_run, executable=executable,
                        scheduler=scheduler,skip=skip, clean_restart=clean_restart,
                        verbose=verbose, **kwargs)
        print('Initialize a QuantumESPRESSO calculator with scheduler %s'%self._global_options['scheduler'])

    def pre_processing(self):
        """
        Process local run dictionary to create the run directory and input file.
        If clean_restart = True the run_dir is cleaned before the run.
        If the 'source_dir' key is passed to the run method copy the source folder
        in the run_dir with the name $prefix.

        """
        run_dir = self.run_options.get('run_dir', '.')
        input= self.run_options.get('input')
        name = self.run_options.get('name','default')
        skip = self.run_options.get('skip')
        clean_restart= self.run_options.get('clean_restart')
        verbose = self.run_options.get('verbose')

        # Create the run_dir and write the input file
        self._ensure_run_directory()
        if input is not None:
            input.write(os.path.join(run_dir,name)+'.in')
        else:
            print('input not provided')

        # Clean the run_dir
        if not skip:
            if clean_restart:
                self._clean_run_dir()
            else:
                if verbose: print('run performed starting from existing results')

        # Copy the source folder in the run_dir
        source_dir = self.run_options.get('source_dir')
        if source_dir is not None:
            self._copy_source_dir(source_dir)

        return {}

    def process_run(self):
        """
        Method associated to the running of the executable. The method runs the computation
        and wait the end of the computation before passing to the :meth:`post_processing` method.

        """
        to_run = self.is_to_run()
        if to_run:
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
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options['input']
        prefix = input['control']['prefix'].strip("'")
        prefix += '.save'
        result = os.path.join(run_dir,prefix,'data-file-schema.xml')
        if not os.path.isfile(result):
            print('Expected file %s not found'%result)
            print('Check if wait_end_run is False or the dry_run option is active. Otherwise a possible error has occured during the computation')
        return result

    def is_to_run(self):
        """
        The method evaluates if the computation can be skipped. This is done by
        checking if the file $prefix.xml is already present in the run_dir.

        Return:
            :py:class:`bool` : the boolean is True if the computation needs to be run
        """
        skip = self.run_options.get('skip')
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options.get('input')
        name = self.run_options.get('name','default')
        verbose = self.run_options.get('verbose')

        if not skip:
            return True
        else:
            prefix = input['control']['prefix'].strip("'")
            skipfile = os.path.join(run_dir,prefix)+'.xml'
            if os.path.isfile(skipfile):
                if verbose: print('Skip the run of',name)
                return False
            else:
                return True

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
        scheduler = self.run_options['scheduler']
        dry_run = self.run_options.get('dry_run',False)
        verbose = self.run_options.get('verbose')

        if scheduler == 'direct':
            # Set the OMP_NUM_THREADS variable in the environment
            os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])
            if not dry_run:
                comm_str = self.run_command()
                job = Popen(comm_str, shell = True)
            else:
                job = None
                if verbose: print('Dry_run option active. Computation not performed')
        elif scheduler == 'slurm':
            job = self.build_slurm_script()
            if not dry_run:
                run_dir = self.run_options.get('run_dir', '.')
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
        dry_run = self.run_options.get('dry_run',False)
        wait_end_run = self.run_options.get('wait_end_run',True)
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
        omp = self.run_options.get('omp')
        mpi = self.run_options.get('mpi')
        name = self.run_options.get('name','default')
        job = 'job_'+name
        run_dir = self.run_options.get('run_dir', '.')
        sbatch_options = self.run_options.get('sbatch_options', None)
        comm_str = self.run_command()

        lines = []
        lines.append('#!/bin/bash')
        lines.append('#SBATCH --ntasks=%s           ### Number of tasks (MPI processes)'%mpi)
        lines.append('#SBATCH --cpus-per-task=%s    ### Number of threads per task (OMP threads)'%omp)
        if sbatch_options is not None: # add other options if present in the run_options of the calculator
            for option in sbatch_options:
                lines.append('#SBATCH %s'%option)
        lines.append('#SBATCH --output=%s.out'%job)
        lines.append('')
        lines.append('export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK')
        lines.append('')
        lines.append('echo "Job id $SLURM_JOB_ID"')
        lines.append('echo "Number of mpi  $SLURM_NTASKS"')
        lines.append('echo "Number of threads per task $SLURM_CPUS_PER_TASK"')
        lines.append('')
        lines.append('echo "execute : %s"'%comm_str)
        lines.append(comm_str)
        lines.append('')
        lines.append('echo "JOB_DONE"')
        f = open(os.path.join(run_dir,job+'.sh'),'w')
        f.write('\n'.join(lines))
        f.close()

        return job

    def run_command(self):
        """
        Define the run command used to run the computation. The value of the command
        depends on the chosen scheduler.

        Return:
            :py:class:`string` : command that runs the computation

        """
        scheduler = self.run_options.get('scheduler')
        executable = self.run_options.get('executable')
        mpi = self.run_options.get('mpi')
        mpi_run = self.run_options.get('mpi_run')
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        verbose = self.run_options.get('verbose')

        if scheduler == 'direct':
            set_run_dir = 'cd %s; '%run_dir
            command = set_run_dir + mpi_run + ' ' + str(mpi) + ' ' + executable
        if scheduler == 'slurm':
            command = mpi_run + ' ' + str(mpi) + ' ' + executable
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

    def _ensure_run_directory(self):
        from mppi.Utilities import FutileUtils as f
        run_dir = self.run_options.get('run_dir', '.')
        # Restrict run_dir to a sub-directory
        if ("/" in run_dir or run_dir == ".."):
            raise ValueError(
                "run_dir '%s' must be a sub-directory"% run_dir)
        # Create the run_dir if not exist
        if f.ensure_dir(run_dir) and self.run_options['verbose']:
            print("Create the sub-directory '%s'" % run_dir)

    def _copy_source_dir(self,source_dir):
        """
        Copy the source_dir in the run_dir and atttibute to the copied folder
        the name $prefix, for all the inputs.

        Args:
            source_dir: the name of the source_dir including its relative path.
            A source_dir outer respect to the actual run_dir of the instance of
            QeCalculator can be used.

        """
        from shutil import copytree
        verbose = self.run_options.get('verbose')
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options.get('input')

        prefix = input['control']['prefix'].strip("'")
        dest_dir = os.path.join(run_dir,prefix)+'.save'
        if not os.path.isdir(dest_dir):
            if verbose: print('Copy source_dir %s in the %s'%(source_dir,dest_dir))
            copytree(source_dir,dest_dir)
        else:
            if verbose:
                print('The folder %s already exsists. Source folder % s not copied'
                %(dest_dir,source_dir))

    def _clean_run_dir(self):
        """
        Clean the run_dir before performing the computation. Delete the $name.log file,
        the $prefix.xml file, the job_$name.out file and the folder run_dir/prefix.save.

        """
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        input = self.run_options.get('input')
        verbose = self.run_options.get('verbose')

        logfile = os.path.join(run_dir,name)+'.log'
        prefix = input['control']['prefix'].strip("'")
        xmlfile = os.path.join(run_dir,prefix)+'.xml'
        job_out = os.path.join(run_dir,'job_'+name+'.out')
        outdir = os.path.join(run_dir,prefix)+'.save'

        if os.path.isfile(logfile):
            if verbose: print('delete log file:',logfile)
            os.system('rm %s'%logfile)
        if os.path.isfile(xmlfile):
            if verbose: print('delete xml file:',xmlfile)
            os.system('rm %s'%xmlfile)
        if os.path.isfile(job_out):
            if verbose: print('delete job_out script:',job_out)
            os.system('rm %s'%job_out)
        if os.path.isdir(outdir):
            if verbose: print('delete folder:',outdir)
            os.system('rm -r %s'%outdir)
