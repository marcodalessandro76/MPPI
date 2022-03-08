"""
This module manages calculations performed with Yambo.
The run of the computation is performed by the python subprocess package (direct scheduler)
or by the slurm scheduler.
"""

from .Runner import Runner
from mppi.Utilities import Tools
from mppi.Calculators.RunRules import build_slurm_header, mpi_command
import os

def get_report(path):
    """
    Look for the name of the r-* file(s) produced by the execution of the code.

    Note that multiple instances of the report could be found, if the computations
    have been repeated more than once with the `clean_restart=False` option.

    Args:
        path (:py:class:`string`) : folder with the r-* and the o-* files

    Return:
        :py:class:`list` : A list with the name, including the path, of the
        file r-* produced by the run. If no report files are found a list with a
        single (empty) string is produced.

    """
    report = []
    if os.path.isdir(path):
        for file in os.listdir(path):
            if 'r-' in file:
                report.append(os.path.join(path,file))
    if len(report) == 0 : report.append(' ')
    return report

def get_output_files(path):
    """
    Look for the names of the o- file(s) produced by the execution of the code.

    Args:
        path (:py:class:`string`) : folder where the o-* files are stored

    Return:
        :py:class:`list` : A list with the names, including the path, of the
        files o-* produced by the run
    """

    output = []
    if os.path.isdir(path):
        for file in os.listdir(path):
            if 'o-' in file:
                output.append(os.path.join(path,file))
    return output

def get_db_files(path):
    """
    Look for the files of tht type ndb.* in the dbsPath.

    Args:
        dbstPath (:py:class:`list`) : list of folders in which the ndb databases created by yambo
            are sought

    Return:
        :py:class:`dict`: A dictionary in which the keys are the extension of the
            databases found and the values are their names, including the path
    """

    dbs = {}
    for dir in path:
        if os.path.isdir(dir):
            for file in os.listdir(dir):
                if 'ndb' in file and 'fragment' not in file:
                    key = file.split('.')[-1]
                    dbs[key] = os.path.join(dir,file)
    return dbs

def build_results_dict(run_dir, outputPath, dbsPath = None, verbose = True):
    """
    Return a dictionary with the names of the o- file(s), the report 'r-' file,
    the `ns.db1` in the SAVE folder and the names of the ndb database written by
    yambo in the jobname folder. ndb database are sought in the outputPath and in
    the dbsPath.

    Args:
        run_dir (:py:class:`string`) : `run_dir` folder of the calculation
        outputPath (:py:class:`string`) : folder with the 'o-' files
        dbsPath (:py:class:`list`) : list of folders with the ndb databases. If it is
            None the databases are sought in the outputPath
        verbose (:py:class:`bool`) : set the amount of information provided on terminal

    Return:
        :py:class:`dict` : the dictionary
            {'output' : [o-1,o-2,...],'dft':...,'dipoles':..., ....}

    """
    if dbsPath is None:
        dbsPath = [outputPath]
    else:
        dbsPath += [outputPath]
    results = dict(output=get_output_files(outputPath),report=get_report(outputPath))
    if verbose and len(results['output']) == 0:
        print("""
        There are no o-* files.
        Maybe you have performed a ypp computation or wait_end_run and/or
        the dry_run option are active.
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
    Perform a Yambo calculation. The parameters used to define the parellelization
    strategy are provided in the `runRules` object.

    Parameters:
       runRulues (:class:`RunRules`) : instance of the :class:`RunRules` class
       executable (:py:class:`string`) : set the executable (yambo, ypp, yambo_rt, ...) of the Yambo package
       skip (:py:class:`bool`) : if True evaluate if the computation can be skipped. This is done by checking that the
            report file built by yambo exsists and contains the string `game_over`, defined as a data member of this class
       clean_restart (:py:class:`bool`) : if True delete the folder with the output files and the database before running the computation
       dry_run (:py:class:`bool`) : with this option enabled the calculator setup the calculations and write the script
            for submitting the job, but the computations are not run
       wait_end_run (:py:class:`bool`) : with this option disabled the run method does not wait the end of the run.
            This option may be useful for interacting with the code in particular in _asincronous_ computation managed
            by the slurm scheduler
       activate_BeeOND (:py:class:`bool`) :  if True set I/O of the run in the BeeOND_dir created by the slurm scheduler.
            The value of the ``BeeOND_dir`` is written as a data member of the class and can be modified if needed
       verbose (:py:class:`bool`) : set the amount of information provided on terminal
       kwargs : other parameters that are stored in the _global_options dictionary

     Computations are performed in the folder specified by the ``run_dir`` parameter. The ``name`` parameter is
     used as name of the yambo input and as the name of the folder where yambo writes the o- `output` files.
     The ``jobname`` parameter is the name of the folder where yambo writes the .ndb databases. If this parameter
     is not provided in the run method the assumption jobname=name is made by the calculator.

     Example:
        >>> rr = RunRules(scheduler='slurm',ntasks_per_node=4,memory='124GB')
        >>> code = YamboCalculator(rr,executable='yambo',skip=True,verbose=True)
        >>> code.run(input = ..., run_dir = ...,name = ...,jobname = ..., **kwargs)

        When the run method is called the class runs the command:
            cd run_dir ; `mpirun command` executable_name -F name.in -J jobname -C name - O out_dir

     where the arguments of the run method are:

    Args:
        run_dir (:py:class:`string`) : the folder in which the simulation is performed
        input (:py:class:`string`) : instance of the :class:`YamboInput` class
            that define the input objects
        name (:py:class:`string`) : string with the names associated to the input file (without extension).
            This string is used also as the name of the folder in which results are written (argument of the -C option of yambo) as
            well as a part of the name of the output files
        jobname (:py:class:`list` or :py:class:`string`) : string (or list of strings) with the value(s) of the jobname folders
            (argument of the -J option of yambo). The first element is the folder name, where yambo writes the database.
            The other values (if provided) are the folders where yambo seeks for pre existing databases. All the elements of the
            list are assumed to be located in the  ``run_dir`` of the calculator. If this variable is not specified the value of
            name is attributed to jobname
        out_dir (:py:class:`string`) : position of the folder in which the $jobname folder is located. This parameter
            is automatically set by the calculator the value of ``BeeOND_dir`` if the option `activate_BeeOND` is enabled.
            Otherwise all the folders are written in the ``run_dir``
        kwargs : other parameters that are stored in the run_options dictionary

    """
    BeeOND_dir = '/mnt/${SLURM_JOB_USER}-jobid_${SLURM_JOB_ID}'
    game_over = 'Game Over & Game summary'
    time_profile = 'Time-Profile'

    def __init__(self,runRules, executable = 'yambo', skip = True, clean_restart = True,
                 dry_run = False, wait_end_run = True, activate_BeeOND = False,
                 verbose = True, **kwargs):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self,**runRules,executable=executable,skip=skip,clean_restart=clean_restart,
                        dry_run=dry_run,wait_end_run=wait_end_run,activate_BeeOND=activate_BeeOND,
                        verbose=verbose,**kwargs)
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
        clean_restart= self.run_options.get('clean_restart')
        verbose = self.run_options.get('verbose')
        SAVE = os.path.join(run_dir,'SAVE')
        self.is_to_run()
        is_to_run = self.run_options.get('is_to_run')

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

        # clean the run dir
        if is_to_run:
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
        is_to_run = self.run_options.get('is_to_run')

        if is_to_run:
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

        Note that the first element of the report list is used, so caution is needed if
        there is more than one report file in the report key of the results dictionary.

        Return:
            :py:class:`dict` : the dictionary
                {'output' : [o-1,o-2,...],'report':...,'dft':...,'dipoles':..., ....}

        """
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        outputPath = os.path.join(run_dir,name)
        jobname = self.run_options.get('jobname',name)
        verbose = self.run_options.get('verbose')

        if type(jobname) == str : dbsPath = [os.path.join(run_dir,jobname)]
        if type(jobname) == list :
            dbsPath = [os.path.join(run_dir,j) for j in jobname]

        results = build_results_dict(run_dir,outputPath,dbsPath=dbsPath,verbose=verbose)
        if verbose:
            report = results['report'][0]
            if Tools.find_string_file(report,self.game_over) is None:
                print('game_over string not found in report. Check the computation!')
            time_sim = Tools.find_string_file(report,self.time_profile)
            if  self.run_options['is_to_run'] is True and time_sim is not None:
                print('Run performed in %s'%time_sim.split()[-1])

        return results

    def is_to_run(self):
        """
        The method evaluates if the computation has to be performed. If ``skip`` is
        False the run is always performed, instead if ``skip`` is True the method
        checks if the report file exsists and contains the string `game_over`
        defined as a member of the class.
        The method adds the key `is_to_run` to ``the run_options`` of the class.

        Note that the first element of the report list is used, so caution is needed if
        there is more than one report file in the report key of the results dictionary.

        """
        skip = self.run_options.get('skip')
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        outputPath = os.path.join(run_dir,name)
        report = get_report(outputPath)[0]
        verbose = self.run_options.get('verbose')

        if not skip:
            self.run_options['is_to_run'] = True
        else:
            if Tools.find_string_file(report,self.game_over) is not None:
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
            os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp_num_threads'])
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
        run_dir = os.path.abspath(self.run_options.get('run_dir', '.'))
        name = self.run_options.get('name','default')
        jobname = self.run_options.get('jobname',name)
        activate_BeeOND = self.run_options.get('activate_BeeOND')
        job = 'job_'+name
        comm_str = self.run_command()

        if type(jobname) == str : dbs_dir = os.path.abspath(os.path.join(run_dir,jobname))
        if type(jobname) == list : dbs_dir = os.path.abspath(os.path.join(run_dir,jobname[0]))

        lines = build_slurm_header(self.run_options)
        if activate_BeeOND:
            lines.append('export BEEOND_DIR=%s'%self.BeeOND_dir)
            lines.append('export RUN_DIR=%s'%run_dir)
            lines.append('echo "The BeeOND option is activated. The I/O of the .ndb database is performed in $BEEOND_DIR"')
            lines.append('echo "RUN_DIR path is $RUN_DIR"')
            lines.append('echo "BEEOND_DIR path is $BEEOND_DIR"')
            lines.append('if [ ! -d $BEEOND_DIR ]; then')
            lines.append('echo "$BEEOND_DIR not found!"')
            lines.append('exit')
            lines.append('fi')
            lines.append('echo " "')
            lines.append('')
            # If the jobname foLder with pre-existing ndb is found, it is copied in the BeeOND_dir path
            if os.path.isdir(dbs_dir):
                lines.append('export DBS_DIR=%s'%dbs_dir)
                lines.append('echo "found DBS_DIR folder $DBS_DIR. Copy the DBS_DIR in the $BEEOND_DIR folder"')
                lines.append('echo "rsync -azv --no-o --no-g --omit-dir-times $DBS_DIR/ $BEEOND_DIR"')
                lines.append('rsync -azv --no-o --no-g --omit-dir-times $DBS_DIR/ $BEEOND_DIR')
                lines.append('echo " "')
                lines.append('')

            lines.append('echo "Add the option -O $BEEOND_DIR to the run command of the calculator"')
            comm_str += ' -O %s'%self.BeeOND_dir
            lines.append('echo " "')
            lines.append('')

        # We do not show the " character in the of the comm_str in the echo line
        lines.append('echo "execute : %s"'%(comm_str.replace('"','')))
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
        executable = self.run_options.get('executable')
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name')
        jobname = self.run_options.get('jobname',name)
        verbose = self.run_options.get('verbose')

        if type(jobname) == str : J = jobname
        if type(jobname) == list : J = '"'+','.join(jobname)+'"'

        mpi_run = mpi_command(self.run_options)
        command = mpi_run + ' ' + executable
        input_name = name + '.in'
        comm_str =  command + ' -F %s -J %s -C %s'%(input_name,J,name)
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
        file, the out_dir and folder and the folder with the databases. If several folders
        are provided in the jobname field, only the first one is deleted, since the others are
        assumed to contain databases obtained from other computations.

        """
        run_dir = self.run_options.get('run_dir','.')
        name = self.run_options.get('name')
        jobname = self.run_options.get('jobname',name)
        verbose = self.run_options.get('verbose')
        job_out = os.path.join(run_dir,'job_'+name+'.out')
        out_dir = os.path.join(run_dir,name)
        if type(jobname) == list:
            ndb_dir = os.path.join(run_dir,jobname[0])
        else:
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
