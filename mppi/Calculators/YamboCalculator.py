"""
This module manages calculations performed with Yambo.
The run of the computation is performed by the python subprocess package (direct scheduler)
or by the slurm scheduler.
"""

from .Runner import Runner
import os

def get_report(outputPath):
    """
    Look for the name of the 'r-*' file produced by the execution of the code.

    Args:
        outputPath (:py:class:`string`) : folder with the r-* and the 'o-*' files

    Return:
        :py:class:`string` : A string with the name, including the path, of the
        file 'r-*' produced by the run

    """

    report = ''
    if os.path.isdir(outputPath):
        for file in os.listdir(outputPath):
            if 'r-' in file:
                report = os.path.join(outputPath,file)
    return report

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

    Args:
        dbstPath (:py:class:`string`) : folder with the databases created by yambo

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
    Return a dictionary with the names of the o- file(s), the report 'r-' file,
    the `ns.db1` in the SAVE folder and the names of the ndb databases written by
    yambo in the jobname folder.

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
    results = dict(output=get_output_files(outputPath),report=get_report(outputPath))
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
       skip (:py:class:`bool`) : if True evaluate if the computation can be skipped. This is done by checking that the
            report file built by yambo exsists and contains the string `game_over` defined as a data member of this class
       clean_restart (:py:class:`bool`) : if True delete the folder with the output files and the database before running the computation
       dry_run (:py:class:`bool`) : with this option enabled the calculator setup the calculations and write the script
            for submitting the job, but the computations are not run
       wait_end_run (:py:class:`bool`) : with this option disabled the run method does not wait the end of the run.
            This option may be useful for interacting with the code in particular in _asincronous_ computation managed
            by the slurm scheduler
       sbatch_options (:py:class:`list`) : the elements of this list are strings used as options in the slurm script.
            For instance it is possible to specify the number of tasks per node as `--ntasks-per-node=16`
       activate_BeeOND (:py:class:`bool`) :  if True set I/O of the run in the BeeOND_dir created by the slurm scheduler.
            With this options enabled ....
       verbose (:py:class:`bool`) : set the amount of information provided on terminal
       kwargs : other parameters that are stored in the _global_options dictionary

     Computations are performed in the folder specified by the ``run_dir`` parameter. The ``name`` parameter is
     used as name of the yambo input and as the name of the folder where yambo writes the o- `output` files.
     The ``jobname`` parameter is the name of the folder where yambo writes the .ndb databases. If this parameter
     is not provided in the run method the assumption jobname=name is made by the calculator.

     Example:
     >>> code = YamboCalculator(omp=1,mpi=4,mpi_run='mpirun -np',executable='yambo',skip=True,verbose=True,scheduler='direct')
     >>> code.run(input = ..., run_dir = ...,name = ...,jobname = ..., **kwargs)

     When the run method is called the class runs the command:
        cd run_dir ; mpi_run mpi executable_name -F name.in -J jobname -C name

     where the arguments of the run method are:

    Args:
        run_dir (:py:class:`string`) : the folder in which the simulation is performed
        input (:py:class:`string`) : instance of the :class:`YamboInput` class
            that define the input objects
        name (:py:class:`string`) : string with the names associated to the input file (without extension).
            This string is used also as the name of the folder in which results are written as well as a
            part of the name of the output files.
        jobname (:py:class:`string`) : string with the values of the jobname. If this variable is not specified
            the value of name is attributed to jobname by process_run.
        kwargs : other parameters that are stored in the run_options dictionary

    """
    BeeOND_dir = '/mnt/${SLURM_JOB_USER}-jobid_${SLURM_JOB_ID}'
    game_over = 'Game Over & Game summary'
    time_profile = 'Time-Profile'

    def __init__(self,
                 omp = os.environ.get('OMP_NUM_THREADS', 1), mpi = 2, mpi_run = 'mpirun -np',
                 executable = 'yambo', scheduler = 'direct', skip = True, clean_restart = True,
                 dry_run = False, wait_end_run = True, sbatch_options = [], activate_BeeOND = True,
                 verbose = True, **kwargs):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=omp, mpi=mpi, mpi_run=mpi_run, executable=executable,
                        scheduler=scheduler, skip=skip, clean_restart=clean_restart,
                        dry_run=dry_run,wait_end_run=wait_end_run,sbatch_options=sbatch_options,
                        activate_BeeOND=activate_BeeOND,verbose=verbose, **kwargs)
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
        name = self.run_options.get('name','default')
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
                self.clean_run_dir()
            else:
                if verbose: print('run performed starting from existing results')

        return {}

    def process_run(self):
        """
        Method associated to the running of the executable. The method runs the computation
        and wait the end of the computation before passing to the :meth:`post_processing` method.

        """
        self.is_to_run()
        if self.run_options['is_to_run']:
            job = self.run_job()
            self.wait(job)
        return {}

    def post_processing(self):
        """
        Return a dictionary with the names of the o- file(s), the report 'r-' file,
        the `ns.db1` in the SAVE folder and the names of the ndb databases written by
        yambo in the jobname folder. The method performs a sanity check of the computation
        by checking that the `game_over` string is found in the report and shows the
        time needed to perform the simulation (if the `time_profile` string is found
        in the report).

        Return:
            :py:class:`dict` : the dictionary
                {'output' : [o-1,o-2,...],'report':...,'dft':...,'dipoles':..., ....}

        """
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        outputPath = os.path.join(run_dir,name)
        jobname = self.run_options.get('jobname',name)
        dbsPath = os.path.join(run_dir,jobname)
        verbose = self.run_options.get('verbose')

        results = build_results_dict(run_dir,outputPath,dbsPath=dbsPath,verbose=verbose)
        if verbose:
            report = results['report']
            if find_string_file(report,self.game_over) is None:
                print('Game_over string not found in report. Check the computation!')
            time_sim = find_string_file(report,self.time_profile)
            if  self.run_options['is_to_run'] is True and time_sim is not None:
                print('Run performed in %s'%time_sim.split()[-1])

        return results

    def is_to_run(self):
        """
        The method evaluates if the computation can be skipped. This is done by
        checking if the report file exsists and contains the string `game_over`
        defined as a member of the class.
        The method adds the key `is_to_run` to ``the run_options`` of the class

        """
        skip = self.run_options.get('skip')
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        outputPath = os.path.join(run_dir,name)
        report = get_report(outputPath)
        verbose = self.run_options.get('verbose')

        if not skip:
            self.run_options['is_to_run'] = True
        else:
            if find_string_file(report,self.game_over) is not None:
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
        dry_run = self.run_options.get('dry_run',False)
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
        run_dir = os.path.abspath(self.run_options.get('run_dir', '.'))
        name = self.run_options.get('name','default')
        jobname = self.run_options.get('jobname',name)
        job = 'job_'+name

        sbatch_options = self.run_options.get('sbatch_options')
        activate_BeeOND = self.run_options.get('activate_BeeOND')
        comm_str = self.run_command()

        lines = []
        lines.append('#!/bin/bash')
        lines.append('#SBATCH --ntasks=%s           ### Number of tasks (MPI processes)'%mpi)
        lines.append('#SBATCH --cpus-per-task=%s    ### Number of threads per task (OMP threads)'%omp)
        for option in sbatch_options: # add other SBATCH options
            lines.append('#SBATCH %s'%option)
        lines.append('#SBATCH --output=%s.out'%job)
        lines.append('')

        lines.append('export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK')
        lines.append('export RUN_DIR=%s'%run_dir)
        lines.append('export BEEOND_DIR=%s'%self.BeeOND_dir)
        lines.append('')

        lines.append('echo "Cluster name $SLURM_CLUSTER_NAME"')
        lines.append('echo "Job name $SLURM_JOB_NAME "')
        lines.append('echo "Job id $SLURM_JOB_ID"')
        lines.append('echo "Job nodelist $SLURM_JOB_NODELIST"')
        lines.append('echo "Number of nodes $SLURM_JOB_NUM_NODES"')
        lines.append('echo "Number of mpi $SLURM_NTASKS"')
        lines.append('echo "Number of threads per task $SLURM_CPUS_PER_TASK"')
        lines.append('echo "BEEOND_DIR path is $BEEOND_DIR"')
        lines.append('echo "RUN_DIR path is $RUN_DIR"')
        lines.append('echo " "')
        lines.append('')

        if activate_BeeOND:
            lines.append('echo "THe BeeOND option is activated. The I/O of the .ndb database is performed in $BEEOND_DIR"')
            lines.append('if [ ! -d $BEEOND_DIR ]; then')
            lines.append('echo "$BEEOND_DIR not found!"')
            lines.append('exit')
            lines.append('fi')
            lines.append('echo " "')
            lines.append('')
            lines.append('echo "Add the option -O $BEEOND_DIR to the run command of the calculator"')
            comm_str += ' -O %s'%self.BeeOND_dir
            lines.append('echo " "')
            lines.append('')

        lines.append('echo "execute : %s"'%comm_str)
        lines.append(comm_str)
        lines.append('echo " "')
        lines.append('')

        if activate_BeeOND:
            lines.append('echo "Copy the folder with the ndb database in the $RUN_DIR"')
            lines.append('echo "rsync -azv $BEEOND_DIR/ $RUN_DIR"')
            lines.append('rsync -azv $BEEOND_DIR/ $RUN_DIR')
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
        scheduler = self.run_options.get('scheduler')
        executable = self.run_options.get('executable')
        mpi = self.run_options.get('mpi')
        mpi_run = self.run_options.get('mpi_run')
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name')
        jobname = self.run_options.get('jobname',name)
        verbose = self.run_options.get('verbose')

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

    def clean_run_dir(self):
        """
        Clean the run_dir before performing the computation. Delete the job_$name.out
        file, the out_dir and ndb_dir folders.

        """
        run_dir = self.run_options.get('run_dir','.')
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
