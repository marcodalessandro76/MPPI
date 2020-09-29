"""
This module manages parallel calculations with Yambo.
Both a scheduler like slurm or the python multiprocessing package can be used.
"""

from .Runner import Runner
import os

class YamboCalculator(Runner):
    """
    Manage (multiple) Yambo calculations performed in parallel. Computations
    are managed by a scheduler that, in the actual implementation of the class, can
    be `direct` or `slurm`.

    Parameters:
       omp (:py:class:`int`) : value of the OMP_NUM_THREADS variable
       mpi (:py:class:`int`) : number of mpi processes
       mpi_run (:py:class:`string`) : command for the execution of mpirun, e.g. 'mpirun -np' or 'mpiexec -np'
       executable (:py:class:`string`) : set the executable (yambo, ypp, yambo_rt, ...) of the Yambo package
       scheduler (:py:class:`string`) : choose the scheduler used to submit the job, actually the choices implemented are
            'direct' that runs the computation using the python multiprocessing package and 'slurm' that creates a slurm script
       multiTask  (:py:class:`bool`) : if true a single run_script is built and all the computations are performed in parallel,
            otherwise an independent script is built for each elements of inputs and the computations are performed sequentially
       skip (:py:class:`bool`) : if True evaluate if one (or many) computations can be skipped.
           This is done by checking that the folder where yambo write the results contains at least
           one file 'o-*', for each name in names
       verbose (:py:class:`bool`) : set the amount of information provided on terminal
       IO_time (int) : time step (in second) used by the wait method to check that the job is completed
       kwargs : other parameters that are stored in the _global_options dictionary
       clean_restart (:py:class:`bool`) : if True the delete the folder(s) with the output files and database before
            running the computation

    Example:
     >>> code = YamboCalculator(omp=1,mpi=4,mpi_run='mpirun -np',executable='yambo',skip=True,verbose=True,scheduler='direct')
     >>> code.run(inputs = ..., run_dir = ...,names = ...,jobnames = ...)

     where the arguments of the run method are:

    Args:
        run_dir (:py:class:`string`) : the folder in which the simulation is performed
        inputs (:py:class:`list`) : list with the instances of the :class:`YamboInput` class
            that define the input objects
        names (:py:class:`list`) : list with the names associated to the input files (without extension),
            given in the same order of the inputs list. These strings are used also as the radicals of the
            folders in which results are written as well as a part of the name of the output files.
        jobnames (:py:class:`list`) : list with the values of the jobname. If this variable is not specified
            the value of name is attributed to jobname by process_run.
        kwargs : other parameters that are stored in the run_options dictionary

    When the run method is called the class runs the command:
                executable_name -F name.in -J jobname -C name

    The calculator looks for the following variables in the run_options dictionary. These options
    may be useful for _asincronous_ computation managed the slurm scheduler.

        `dry_run=True` with this option the calculator setup the calculations and write the scrpt
        for submitting the jobs, but the computations are not run.

        `wait_end_run=False` with this option the wait of the end of the run is suppressed.

    """

    def __init__(self,
                 omp = os.environ.get('OMP_NUM_THREADS', 1), mpi = 2, mpi_run = 'mpirun -np',
                 executable = 'yambo', scheduler = 'direct', multiTask = True,
                 skip = True, verbose = True, IO_time = 5, clean_restart = True, **kwargs):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=omp, mpi=mpi, mpi_run=mpi_run, executable=executable,
                        scheduler=scheduler, multiTask=multiTask,
                        skip=skip, verbose=verbose, IO_time=IO_time,
                        clean_restart=clean_restart, **kwargs)
        if multiTask: task_str = 'parallel'
        else: task_str = 'serial'
        print('Initialize a %s Yambo calculator with scheduler %s' %
            (task_str,self._global_options['scheduler']))


    def pre_processing(self):
        """
        Process local run dictionary. Check that the run_dir exists and that it
        contains the SAVE folder. Check that the inputs objects have been provided
        in the run parameters and write the input files on disk.
        If skip = False clean the run_dir.

        Note:
            If the run_dir and/or the SAVE folder do not exist an alert is
            written but the execution of the run method proceedes.

        """
        run_dir = self.run_options.get('run_dir', '.')
        inputs = self.run_options.get('inputs')
        names = self.run_options.get('names')
        SAVE = os.path.join(run_dir,'SAVE')
        skip = self.run_options.get('skip')
        clean_restart= self.run_options.get('clean_restart')
        verbose = self.run_options.get('verbose')

        # check if the run_dir and SAVE folder exist and write the input
        if not os.path.isdir(run_dir):
            print('Run_dir %s does not exists'%run_dir)
        elif not os.path.isdir(SAVE):
            print('SAVE folder does not exists')
        else:
            if inputs is not None:
                for input,name in zip(inputs,names):
                    input.write(run_dir,name+'.in')
            else:
                print('input list not provided')

        # if skip = False clean the run_dir
        if not skip:
            if clean_restart:
                self._clean_run_dir()
            else:
                if verbose: print('run performed starting from existing results')

        return {}

    def process_run(self):
        """
        Method associated to the running of the executable. The method prepares the
        jobs script(s), then submit the jobs and wait the end of the computation before
        passing to the :meth:`post_processing` method. Computations are performed
        in parallel or serially accordingly to the value of the multiTask option.

        """
        multiTask = self.run_options.get('multiTask')
        dry_run = self.run_options.get('dry_run',False)
        wait_end_run = self.run_options.get('wait_end_run',True)

        to_run = self.select_to_run()
        jobs = self.build_run_script(to_run)

        if not dry_run:
            if multiTask:
                self.submit_job(jobs)
                if wait_end_run: self.wait(jobs,to_run)
            else:
                for index,job in zip(to_run,jobs):
                    self.submit_job([job])
                    if wait_end_run: self.wait([job],[index])

        return {}

    def post_processing(self):
        """
        Return a dictionary with the names of the o- file(s) and the name of the
        folder that contains the databases.
        The construction of the lists for the output key is managed by the :meth:_get_output.
        For the folders that contain the databases, the method return None if the corresponding
        folder does not exists.

        Return:
            :py:class:`dict` : the dictionary
                {'output' : [[o-1,o-2,...],[o-1,o-2,...],...], 'dbs' : [ndb_folder1,ndb_folder2,...]}
        """
        run_dir = self.run_options.get('run_dir', '.')
        names = self.run_options.get('names')
        jobnames = self.run_options.get('jobnames',names)

        results = {'output' : [], 'dbs' : []}
        for name,jobname in zip(names,jobnames):
            results['output'].append(self._get_output_files(name))
            db_folder = os.path.join(run_dir,jobname)
            if os.path.isdir(db_folder):
                results['dbs'].append(db_folder)
            else:
                results['dbs'].append(None)

        return results

    def select_to_run(self):
        """
        If the skip attribute of run_options is True the method evaluates which
        computations can be skipped. This is done by checking if the folder where
        yambo write the results  contains at least one file 'o-*'.

        Return:
            :py:class:`list` : list with numbers of the computations that have to
            be performed, in the same order provided in the run method
        """
        skip = self.run_options.get('skip')
        run_dir = self.run_options.get('run_dir', '.')
        names = self.run_options.get('names')
        verbose = self.run_options.get('verbose')

        if not skip:
            to_run = [index for index in range(len(names))]
            return to_run
        else:
            to_run = []
            for index,name in enumerate(names):
                if len(self._get_output_files(name)) > 0:
                    if verbose: print('Skip the run of',name)
                else:
                    to_run.append(index)
            return to_run

    def build_run_script(self,to_run):
        """
        Create the run script(s) that are executed by the :meth:`submit_job` method.
        The scripts depend on the scheduler adopted, and specific methods for
        `direct` and `slurm` scheduler are implemented.

        Args:
            to_run (:py:class:`string`) : list with the cardinal numbers of the runs
                to be performed

        Return:
            :py:class:`list` : list with jobs to run. The type of the object in the
            list depends on the chosen scheduler

        """
        scheduler = self.run_options['scheduler']

        jobs = None
        if scheduler == 'direct':
            jobs = self.direct_scheduler(to_run)
        elif scheduler == 'slurm':
            jobs = self.slurm_scheduler(to_run)
        else:
            print('scheduler unknown')
        return jobs

    def submit_job(self,jobs):
        """
        Submit the job.

        Args:
            jobs : The reference to the jobs to be executed. If the scheduler is `direct`
                jobs is a list with the instance of :py:class:multiprocessing. If the
                scheduler is `slurm` jobs is a list with the names of the slurm scripts
        """
        scheduler = self.run_options['scheduler']

        if scheduler == 'direct':
            # Set the OMP_NUM_THREADS variable in the environment
            os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])
            for run in jobs:
                run.start()
        if scheduler == 'slurm':
            run_dir = self.run_options.get('run_dir', '.')
            for job in jobs:
                slurm_submit = 'cd %s ; sbatch %s.sh' %(run_dir,job)
                print('slurm submit: ',slurm_submit )
                os.system(slurm_submit)

    def wait(self,jobs,to_run):
        """
        Wait the end of the jobs.

        Args:
            jobs : The reference to the jobs to be executed. If the scheduler is `direct`
                jobs is a list with the instance of :py:class:multiprocessing. If the
                scheduler is `slurm` jobs is a list with the names of the slurm scripts

            to_run (:py:class:`string`) : list with the cardinal numbers of the runs
                to be performed

        """
        verbose = self.run_options.get('verbose')
        IO_time = self.run_options.get('IO_time')
        import time

        while not all(self._jobs_terminated(jobs)):
            if verbose:
                s = ''
                for index,status in zip(to_run,self._jobs_terminated(jobs)):
                    s+='run'+str(index)+'_is_running: '+str(not status) + ' '
                print(s)
            time.sleep(IO_time)
        if verbose : print('Job completed')

    def direct_scheduler(self,to_run):
        """
        Define the list of Process (methods of multiprocessing) associated to the
        runs specified in the list to_run.

        Args:
            to_run (:py:class:`string`) : list with the cardinal numbers of the runs
                to be performed

        Return:
            :py:class:`list` : list of the :py:class:`multiprocessing` objects
            associated to the runs of the job

        """
        import multiprocessing
        def os_system_run(comm_str):
            os.system(comm_str)

        jobs = []
        for index in to_run:
            comm_str = self.run_command(index)
            p = multiprocessing.Process(target=os_system_run, args=(comm_str,))
            jobs.append(p)
        return jobs

    def slurm_scheduler(self,to_run):
        """
        Create the slurm script(s) associated to the runs specified in the list to_run.

        Args:
            to_run (:py:class:`string`) : list with the cardinal numbers of the runs
                to be performed

        Return:
            :py:class:`list`: list with the names of the slurm scripts associated to the
            computations that are not skipped

        """
        omp = self.run_options.get('omp')
        mpi = self.run_options.get('mpi')
        names = self.run_options.get('names')
        run_dir = self.run_options.get('run_dir', '.')
        sbatch_options = self.run_options.get('sbatch_options', None)

        lines_options = []
        lines_options.append('#!/bin/bash')
        lines_options.append('#SBATCH --ntasks=%s           ### Number of tasks (MPI processes)'%mpi)
        lines_options.append('#SBATCH --cpus-per-task=%s    ### Number of threads per task (OMP threads)'%omp)
        if sbatch_options is not None: # add other options if present in the run_options of the calculator
            for option in sbatch_options:
                lines_options.append('#SBATCH %s'%option)
        jobs = []
        for index in to_run:
            job_name = 'job_'+names[index]
            jobs.append(job_name)
            comm_str = self.run_command(index)
            lines_run = []
            lines_run.append('#SBATCH --output=%s.out'%job_name)
            lines_run.append('')
            lines_run.append('export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK')
            lines_run.append('')
            lines_run.append('echo "Job id $SLURM_JOB_ID"')
            lines_run.append('echo "Number of mpi  $SLURM_NTASKS"')
            lines_run.append('echo "Number of threads per task $SLURM_CPUS_PER_TASK"')
            lines_run.append('')
            lines_run.append('echo "execute : %s"'%comm_str)
            lines_run.append(comm_str)
            lines_run.append('')
            lines_run.append('echo "JOB_DONE"')
            f = open(os.path.join(run_dir,job_name+'.sh'),'w')
            f.write('\n'.join(lines_options+lines_run))
            f.close()

        return jobs

    def run_command(self,index):
        """
        Define the run command used to run the computation associated to the
        input file $names[index]. The value of the command depends on the
        chosen scheduler.

        Args:
            index (:py:class:`int`) : index of the computation to be performed

        Return:
            :py:class:`string` : command that runs the computation associated to
            the $names[index] input file
        """
        scheduler = self.run_options.get('scheduler')
        executable = self.run_options.get('executable')
        mpi = self.run_options.get('mpi')
        mpi_run = self.run_options.get('mpi_run')
        run_dir = self.run_options.get('run_dir', '.')
        names = self.run_options.get('names')
        jobnames = self.run_options.get('jobnames',names)
        verbose = self.run_options.get('verbose')

        if scheduler == 'direct':
            set_run_dir = 'cd %s; '%run_dir
            command = set_run_dir + mpi_run + ' ' + str(mpi) + ' ' + executable
        if scheduler == 'slurm':
            command = mpi_run + ' ' + str(mpi) + ' ' + executable

        input_name = names[index] + '.in'
        comm_str =  command + ' -F %s -J %s -C %s'%(input_name,jobnames[index],names[index])
        if verbose: print('run %s command: %s' %(index,comm_str))

        return comm_str

    def _jobs_terminated(self,jobs):
        """
        Check the status of the running jobs.

        Args:
            jobs (:py:class:`list`) : list with the reference to the running jobs

        Return:
            :py:class:`list`: list with the status of the jobs. The elements are True
                if the associated computation is terminated and False if it is running

        """
        scheduler = self.run_options.get('scheduler')
        run_dir = self.run_options.get('run_dir', '.')

        if scheduler == 'direct':
            jobs_terminated = [not job.is_alive() for job in jobs]
        if scheduler == 'slurm':
            jobs_terminated = []
            for job in jobs:
                job_out = os.path.join(run_dir,job+'.out')
                if not os.path.isfile(job_out):
                    jobs_terminated.append(False)
                else:
                    with open(job_out, 'r') as f:
                        last_line = f.read().splitlines()[-1]
                    if last_line == 'JOB_DONE': jobs_terminated.append(True)
                    else: jobs_terminated.append(False)

        return jobs_terminated

    def _get_output_files(self,name):
        """
        Look for the names of the 'o-' file(s) produced by the execution of the
        $name input file.

        Return:
            :py:class:`list`: A list with the names, including the path, of the
            files 'o-*' produced by the run
        """
        run_dir = self.run_options.get('run_dir', '.')
        out_dir = os.path.join(run_dir,name)

        output = []
        if os.path.isdir(out_dir):
            for file in os.listdir(out_dir):
                if 'o-' in file:
                    output.append(os.path.join(out_dir,file))
        return output

    def _clean_run_dir(self):
        """
        Clean the run_dir before performing the computation. Delete the out_dir
        and ndb_dir folders (if found).
        """
        run_dir = self.run_options.get('run_dir', '.')
        names = self.run_options.get('names')
        jobnames = self.run_options.get('jobnames',names)
        verbose = self.run_options.get('verbose')

        for name,jobname in zip(names,jobnames):
            out_dir = os.path.join(run_dir,name)
            ndb_dir = os.path.join(run_dir,jobname)

            if os.path.isdir(out_dir):
                if verbose: print('delete folder:',out_dir)
                os.system('rm -r %s'%out_dir)
            if os.path.isdir(ndb_dir):
                if verbose: print('delete folder:',ndb_dir)
                os.system('rm -r %s'%ndb_dir)
