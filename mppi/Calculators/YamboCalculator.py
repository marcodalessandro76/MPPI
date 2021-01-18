"""
This module manages calculations performed with Yambo.
Actually the run of the computation can be managed by the python subprocess package (direct scheduler)
or by the slurm scheduler.
"""

from .Runner import Runner
import os

def get_output_files(outputPath):
    """
    Look for the names of the 'o-' file(s) produced by the execution of the code.

    Args:
        outputPath (:py:class:`string`) : folder with the 'o-*' files

    Return:
        :py:class:`list` : A list with the names, including the path, of the
        files 'o-*' produced by the run

    """

    output = []
    if os.path.isdir(outputPath):
        for file in os.listdir(outputPath):
            if 'o-' in file:
                output.append(os.path.join(outputPath,file))
    return output

def get_db_files(dbsPath):
    """
    Look for the files of tht type ndb.* in the dbsPath

    Return:
        :py:class:`dict`: A dictionary in which the keys are the extension of the
            databases found and the values are their names, including the path

    """

    dbs = {}
    if os.path.isdir(dbsPath):
        for file in os.listdir(dbsPath):
            if 'ndb' in file and 'fragment' not in file:
                key = file.split('.')[-1]
                dbs[key] = os.path.join(dbsPath,file)
    return dbs

def build_results_dict(run_dir, outputPath, dbsPath = None, verbose = True):
    """
    Return a dictionary with the names of the o- file(s), the `ns.db1` in the
    SAVE folder and the names of the ndb databases written by yambo in the
    jobname folder .

    Args:
        run_dir (:py:class:`string`) : `run_dir` folder of the calculation
        outputPath (:py:class:`string`) : folder with the 'o-' files
        dbsPath (:py:class:`string`) : folder with the ndb databases. If it is
            None the databases are sought in the outputPath
        verbose (:py:class:`bool`) : set the amount of information provided on terminal

    Return:
        :py:class:`dict` : the dictionary
            {'output' : [o-1,o-2,...],'dft':...,'dipoles':..., ....}

    """
    if dbsPath is None:
            dbsPath = outputPath
    results = dict(output=get_output_files(outputPath))
    if verbose and len(results['output']) == 0:
        print("""
        There are no o-* files.
        Maybe you have performed a ypp computation or wait_end_run and/or
        the dry_run option are active?
        Otherwise a possible error has occured during the computation
        """)
    # add the dft database from the SAVE folder
    dft = os.path.join(run_dir,'SAVE','ns.db1')
    if os.path.isfile(dft):
        results['dft'] = dft
    else:
        if verbose: print('ns.db1 database not found in SAVE folder')
    # add the dbs found in the jobname folder
    results.update(get_db_files(dbsPath))
    return results

class YamboCalculator(Runner):
    """
    Perform a Yambo calculation. Computations are managed by a scheduler that,
    in the actual implementation of the class, can be `direct` or `slurm`.

    Parameters:
       omp (:py:class:`int`) : value of the OMP_NUM_THREADS variable
       mpi (:py:class:`int`) : number of mpi processes
       mpi_run (:py:class:`string`) : command for the execution of mpirun, e.g. 'mpirun -np' or 'mpiexec -np'
       executable (:py:class:`string`) : set the executable (yambo, ypp, yambo_rt, ...) of the Yambo package
       scheduler (:py:class:`string`) : choose the scheduler used to submit the job, actually the choices implemented are
            'direct' that runs the computation using the python subprocess package and 'slurm' that creates a slurm script
       skip (:py:class:`bool`) : if True evaluate if the computation can be skipped. This is done by checking that the folder
            where yambo write the results contains at least one file 'o-*'
       clean_restart (:py:class:`bool`) : if True delete the folder with the output files and the database before running the computation
       verbose (:py:class:`bool`) : set the amount of information provided on terminal
       kwargs : other parameters that are stored in the _global_options dictionary

     Example:
     >>> code = YamboCalculator(omp=1,mpi=4,mpi_run='mpirun -np',executable='yambo',skip=True,verbose=True,scheduler='direct')
     >>> code.run(input = ..., run_dir = ...,name = ...,jobname = ..., **kwargs)

     where the arguments of the run method are:

    Args:
        run_dir (:py:class:`string`) : the folder in which the simulation is performed
        input (:py:class:`string`) : instance of the :class:`YamboInput` class
            that define the input objects
        name (:py:class:`string`) : string with the names associated to the input file (without extension).
            This string is used also as the radical of the folders in which results are written as well as a
            part of the name of the output files.
        jobname (:py:class:`string`) : string with the values of the jobname. If this variable is not specified
            the value of name is attributed to jobname by process_run.
        kwargs : other parameters that are stored in the run_options dictionary

    When the run method is called the class runs the command:
                executable_name -F name.in -J jobname -C name

    The calculator looks for the following variables in the run_options dictionary. These options
    may be useful for _asincronous_ computation managed the slurm scheduler.

        `dry_run=True` with this option the calculator setup the calculations and write the scrpt
        for submitting the jobs, but the computations are not run.

        `wait_end_run=False` with this option the wait of the end of the run is suppressed.

        `sbatch_options = [option1,option2,....]` allows the user to include further options in the slurm script

    """

    def __init__(self,
                 omp = os.environ.get('OMP_NUM_THREADS', 1), mpi = 2, mpi_run = 'mpirun -np',
                 executable = 'yambo', scheduler = 'direct', skip = True, clean_restart = True,
                 verbose = True, **kwargs):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=omp, mpi=mpi, mpi_run=mpi_run, executable=executable,
                        scheduler=scheduler, skip=skip, clean_restart=clean_restart,
                        verbose=verbose, **kwargs)
        print('Initialize a Yambo calculator with scheduler %s' %self._global_options['scheduler'])


    def pre_processing(self):
        """
        Process local run dictionary. Check that the run_dir exists and that it
        contains the SAVE folder.
        If clean_restart = True the run_dir is cleaned before the run.

        Note:
            If the run_dir and/or the SAVE folder do not exist an alert is
            written but the execution of the run method proceedes.

        """
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options.get('input')
        name = self.run_options.get('name')
        skip = self.run_options.get('skip')
        clean_restart= self.run_options.get('clean_restart')
        verbose = self.run_options.get('verbose')
        SAVE = os.path.join(run_dir,'SAVE')

        # check if the run_dir and SAVE folder exist and write the input
        if not os.path.isdir(run_dir):
            print('Run_dir %s does not exists'%run_dir)
        elif not os.path.isdir(SAVE):
            print('SAVE folder does not exists')
        else:
            if input is not None:
                input.write(run_dir,name+'.in')
            else:
                print('input not provided')

        # Clean the run dir
        if not skip:
            if clean_restart:
                self._clean_run_dir()
            else:
                if verbose: print('run performed starting from existing results')

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
        Return a dictionary with the names of the o- file(s), the `ns.db1` in the
        SAVE folder and the names of the ndb databases written by yambo in the
        jobname folder .

        Return:
            :py:class:`dict` : the dictionary
                {'output' : [o-1,o-2,...],'dft':...,'dipoles':..., ....}

        """
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name')
        outputPath = os.path.join(run_dir,name)
        jobname = self.run_options.get('jobname',name)
        dbsPath = os.path.join(run_dir,jobname)
        verbose = self.run_options.get('verbose')

        return build_results_dict(run_dir,outputPath,dbsPath=dbsPath,verbose=verbose)

    def is_to_run(self):
        """
        The method evaluates if the computation can be skipped. This is done by
        checking if the folder where yambo write the results contains at least one file 'o-*'.

        Return:
            :py:class:`bool` : the boolean is True if the computation needs to be run
        """
        skip = self.run_options.get('skip')
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name')
        outputPath = os.path.join(run_dir,name)
        verbose = self.run_options.get('verbose')

        if not skip:
            return True
        else:
            if len(get_output_files(outputPath)) > 0:
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
        name = self.run_options.get('name')
        jobname = self.run_options.get('jobname',name)
        verbose = self.run_options.get('verbose')

        if scheduler == 'direct':
            set_run_dir = 'cd %s; '%run_dir
            command = set_run_dir + mpi_run + ' ' + str(mpi) + ' ' + executable
        if scheduler == 'slurm':
            command = mpi_run + ' ' + str(mpi) + ' ' + executable

        input_name = name + '.in'
        comm_str =  command + ' -F %s -J %s -C %s'%(input_name,jobname,name)
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

    def _clean_run_dir(self):
        """
        Clean the run_dir before performing the computation. Delete the job_$name.out
        file, the out_dir and ndb_dir folders.

        """
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name')
        jobname = self.run_options.get('jobname',name)
        verbose = self.run_options.get('verbose')

        job_out = os.path.join(run_dir,'job_'+name+'.out')
        out_dir = os.path.join(run_dir,name)
        ndb_dir = os.path.join(run_dir,jobname)

        if os.path.isfile(job_out):
            if verbose: print('delete job_out script:',job_out)
            os.system('rm %s'%job_out)
        if os.path.isdir(out_dir):
            if verbose: print('delete folder:',out_dir)
            os.system('rm -r %s'%out_dir)
        if os.path.isdir(ndb_dir):
            if verbose: print('delete folder:',ndb_dir)
            os.system('rm -r %s'%ndb_dir)
