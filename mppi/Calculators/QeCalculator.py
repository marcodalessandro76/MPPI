"""
This module manages parallel calculations with QuantumESPRESSO.
Both a scheduler like slurm or the python multiprocessing package can be used.
"""

from .Runner import Runner
import os

class QeCalculator(Runner):
    """
    Manage (multiple) QuantumESPRESSO calculations performed in parallel. Computations
    are managed by a scheduler that, in the actual implementation of the class, can
    be `direct` or `slurm`.

    Parameters:
       omp (:py:class:`int`) : value of the OMP_NUM_THREADS variable
       mpi (:py:class:`int`) : number of mpi processes
       mpi_run (:py:class:`string`) : command for the execution of mpirun, e.g. 'mpirun -np' or 'mpiexec -np'
       executable (:py:class:`string`) : set the executable (pw.x, ph.x, ..) of the QuantumESPRESSO package
       scheduler (:py:class:`string`) : choose the scheduler used to submit the job, actually the choices implemented are
            'direct' that runs the computation using the python multiprocessing package and 'slurm' that creates a slurm script
       multiTask  (:py:class:`bool`) : if true a single run_script is built and all the computations are performed in parallel,
            otherwise an independent script is built for each elements of inputs and the computations are performed sequentially
       skip (:py:class:`bool`) : if True evaluate if one (or many) computations can be skipped.
           This is done by checking if the file $name.xml is present in the prefix folder,
           for each name in names
       verbose (:py:class:`bool`) : set the amount of information provided on terminal
       IO_time (int) : time step (in second) used by the wait method to check that the job is completed
       kwargs : other parameters that are stored in the _global_options dictionary

    Example:
     >>> code = calculator(omp=1,mpi=4,mpi_run='mpirun -np',skip=True,verbose=True,scheduler='direct')
     >>> code.run(inputs = [...], run_dir = ...,names = [...], source_dir = ..., **kwargs)

     where the arguments of the run method are:

    Args:
        run_dir (:py:class:`string`) : the folder in which the simulation is performed
        inputs (:py:class:`list`) : list with the instances of the :class:`PwInput` class
            that define the input objects
        names (:py:class:`list`) : list with the names associated to the input files,
            given in the same order of the inputs list.
            Usually you can set the name equal to the prefix of the input object so
            the name of the input file and the prefix folder built by QuantumESPRESSO
            are equal
        source_dir (:py:class:`string`) : location of the scf source folder for a nscf computation.
            If present the class copies this folder in the run_dir with the name $prefix.save
        kwargs : other parameters that are stored in the run_options dictionary

    """

    def __init__(self,
                 omp = os.environ.get('OMP_NUM_THREADS', 1), mpi = 2, mpi_run = 'mpirun -np',
                 executable = 'pw.x', scheduler = 'direct', multiTask = True,
                 skip = False, verbose = True, IO_time = 5, **kwargs):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=omp, mpi=mpi, mpi_run=mpi_run, executable=executable,
                        scheduler=scheduler, multiTask=multiTask,
                        skip=skip, verbose=verbose, IO_time=IO_time, **kwargs)
        if multiTask: task_str = 'parallel'
        else: task_str = 'serial'
        print('Initialize a %s QuantumESPRESSO calculator with scheduler %s' %
            (task_str,self._global_options['scheduler']))

    def pre_processing(self):
        """
        Process local run dictionary to create the run directory and input files.
        If skip = False clean the run_dir.
        If the 'source_dir' key is passed to the run method copy the source folder
        in the run_dir with the name $prefix. This procedure is performed
        after the deletion run_dir/prefix.save since otherwise the copy of the
        source_dir is deleted.

        """
        run_dir = self.run_options.get('run_dir', '.')
        inputs = self.run_options.get('inputs')
        names = self.run_options.get('names')
        skip = self.run_options.get('skip')

        # Create the run_dir and write the input file
        self._ensure_run_directory()
        if inputs is not None:
            for input,name in zip(inputs,names):
                input.write(os.path.join(run_dir,name)+'.in')
        else:
            print('input list not provided')

        # if skip = False clean the run_dir
        if not skip:
            self._clean_run_dir()

        # Copy the source folder in the run_dir
        source_dir = self.run_options.get('source_dir')
        if source_dir is not None:
            self._copy_source_dir(source_dir)

        return {}

    def process_run(self):
        """
        Method associated to the running of the executable. The method prepares the
        jobs script(s), then submit the jobs and wait the end of the computation before
        passing to the :meth:`post_processing` method. Computations are performed
        in parallel or serially accordingly to the value of the multiTask option.

        """
        multiTask = self.run_options.get('multiTask')

        to_run = self.select_to_run()
        jobs = self.build_run_script(to_run)

        if multiTask:
            self.submit_job(jobs)
            self.wait(jobs,to_run)
        else:
            for index,job in zip(to_run,jobs):
                self.submit_job([job])
                self.wait([job],[index])

        return {}

    def post_processing(self):
        """
        Return a list with the names, including the path, of the data-file-schema.xml
        files for each element of inputs. If a file is absent the method returns
        None in the associated element of the list, making easy to understand if a specific
        computation has been correctly performed.

        Return:
            :py:class:`dict` : dictionary
                {'output' : []}
             where [] is a list with the names of the xml files (if the file exists) otherwise the associated element is set to None.

        """
        run_dir = self.run_options.get('run_dir', '.')
        inputs = self.run_options['inputs']
        results = {'output' : []}
        for input in inputs:
            prefix = input['control']['prefix'].strip("'")
            prefix += '.save'
            result = os.path.join(run_dir,prefix,'data-file-schema.xml')
            if os.path.isfile(result):
                results['output'].append(result)
            else:
                results['output'].append(None)

        return results

    def select_to_run(self):
        """
        If the skip attribute of run_options is True the method evaluates which
        computations can be skipped. This is done by  checking if the file
        $prefix.xml is already present in the run_dir.

        Return:
            :py:class:`list` : list with numbers of the computations that have to
            be performed, in the same order provided in the run method
        """
        skip = self.run_options.get('skip')
        run_dir = self.run_options.get('run_dir', '.')
        inputs = self.run_options.get('inputs')
        names = self.run_options.get('names')
        verbose = self.run_options.get('verbose')

        if not skip:
            to_run = [index for index in range(len(inputs))]
            return to_run
        else:
            to_run = []
            for index,input in enumerate(inputs):
                prefix = input['control']['prefix'].strip("'")
                skipfile = os.path.join(run_dir,prefix)+'.xml'
                if os.path.isfile(skipfile):
                    if verbose: print('Skip the run of',names[index])
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
        lines_options.append('#SBATCH --cpus_per_task=%s    ### Number of threads per task (OMP threads)'%omp)
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
        verbose = self.run_options.get('verbose')

        if scheduler == 'direct':
            set_run_dir = 'cd %s; '%run_dir
            command = set_run_dir + mpi_run + ' ' + str(mpi) + ' ' + executable
        if scheduler == 'slurm':
            command = mpi_run + ' ' + str(mpi) + ' ' + executable

        input_name = names[index] + '.in'
        output_name = names[index] + '.log'

        comm_str =  command + ' -inp %s > %s'%(input_name,output_name)
        if verbose: print('run %s command: %s' %(index,comm_str))

        return comm_str

    # def _is_terminated(self,jobs):
    #     """
    #     Check if all the runs of the jobs list have been performed.
    #     """
    #     scheduler = self.run_options['scheduler']
    #
    #     if scheduler == 'direct':
    #         terminated = all([not run.is_alive() for run in jobs])
    #         return terminated
    #     if scheduler == 'slurm':
    #         print('is_terminated method to be implemented for slurm scheduler')
    #         return True

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
        inputs = self.run_options.get('inputs')

        for input in inputs:
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
        the $prefix.xml file and the folder run_dir/prefix.save associated to all the
        inputs and names.
        """
        run_dir = self.run_options.get('run_dir', '.')
        names = self.run_options.get('names','default')
        inputs = self.run_options.get('inputs')
        verbose = self.run_options.get('verbose')

        for input,name in zip(inputs,names):
            logfile = os.path.join(run_dir,name)+'.log'
            prefix = input['control']['prefix'].strip("'")
            xmlfile = os.path.join(run_dir,prefix)+'.xml'
            outdir = os.path.join(run_dir,prefix)+'.save'

            if os.path.isfile(logfile):
                if verbose: print('delete log file:',logfile)
                os.system('rm %s'%logfile)
            if os.path.isfile(xmlfile):
                if verbose: print('delete xml file:',xmlfile)
                os.system('rm %s'%xmlfile)
            if os.path.isdir(outdir):
                if verbose: print('delete folder:',outdir)
                os.system('rm -r %s'%outdir)
