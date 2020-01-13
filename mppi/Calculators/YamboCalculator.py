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
       executable (:py:class:`string`) : set the executable (yambo, ypp, yambo_rt, ...) of the Yambo package
       multiTask  (:py:class:`bool`) : if true a single run_script is built and all the computations are performed in parallel,
            otherwise an independent script is built for each elements of inputs and the computations are performed sequentially
       mpi_run (:py:class:`string`) : command for the execution of mpirun, e.g. 'mpirun -np 4'
       cpus_per_task (:py:class:`int`) : set the `cpus_per_task` variable for the slurm script
       ntasks (:py:class:`int`) : `set the ntasks` variable for the slurm script
       skip (:py:class:`bool`) : if True evaluate if one (or many) computations can be skipped.
           This is done by checking that the folder where yambo write the results contains at least
           one file 'o-*', for each name in names
       verbose (:py:class:`bool`) : set the amount of information provided on terminal
       IO_time (int) : time step (in second) used by the wait method to check that the job is completed

    Example:
     >>> code = YamboCalculator(omp=1,mpi_run='mpirun -np 4',executable='yambo',skip=True,verbose=True,scheduler='direct')
     >>> code.run(input = ..., run_dir = ...,name = ...,jobname = ...)

     where the arguments of the run method are:

    Args:
        run_dir (:py:class:`string`) : the folder in which the simulation is performed
        inputs (:py:class:`list`) : list with the instances of the :class:`YamboInput` class
            that define the input objects
        names (:py:class:`list`) : list with the names associated to the input files (without extension),
            given in the same order of the inputs list. These strings are used also as the radicals of the
            folders in which results are written as well as a part of the name of the output files.
        jobnames (:py:class:`list`) : list with the values of the jobname. If it left to None the
            value of name is attributed to jobname by process_run.
        kwargs : other parameters that are stored in the run_options dictionary

    When the run method is called the class runs the command:
                executable_name -F name.in -J jobname -C name

    """

    def __init__(self,
                 omp = os.environ.get('OMP_NUM_THREADS', 1), executable = 'yambo',
                 multiTask = True, scheduler = 'direct',
                 mpi_run = 'mpirun -np 2', cpus_per_task = 4, ntasks = 3,
                 skip = False, verbose = True, IO_time = 5):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=omp, executable=executable,
                        multiTask=multiTask, scheduler=scheduler,
                        mpi_run=mpi_run, cpus_per_task=cpus_per_task, ntasks=ntasks,
                        skip=skip, verbose=verbose, IO_time=IO_time)
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
        verbose = self.run_options['verbose']

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
        skip = self.run_options['skip']
        if not skip:
            self._clean_run_dir()

        return {}

    def process_run(self):
        """
        Method associated to the running of the executable. The method prepare the
        job script in a way that depend on the chosen scheduler, then submit the job
        and wait the end of the computation before passing to :meth:`post_processing` method.

        Note:
            The wait method is actually not implemented for the slurm scheduler

        Returns:
           :py:class:`dict`: dictionary with names of the o- files and the name of the dbs folder, for each element of the inputs

        """
        # Set the OMP_NUM_THREADS variable in the environment
        os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])
        names = self.run_options.get('names')
        jobnames = self.run_options.get('jobnames',names)
        multiTask = self.run_options.get('multiTask')

        if multiTask:
            job = self.build_run_script(names,jobnames)
            self.submit_job(job)
            self.wait(job)
        else:
            for name,jobname in zip(names,jobnames):
                job = self.build_run_script([name],[jobname])
                self.submit_job(job)
                self.wait(job)

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
            results['output'].append(self._get_output(name))
            db_folder = os.path.join(run_dir,jobname)
            if os.path.isdir(db_folder):
                results['dbs'].append(db_folder)
            else:
                results['dbs'].append(None)

        return results

    def _get_output(self,name):
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

    def build_run_script(self,names,jobnames):
        """
        Create the run script that is executed by the :meth:`submit_job` method.
        The script depends on the scheduler adopted, and specific methods for
        `direct` and `slurm` scheduler are implemented.

        Args:
            names (:py:class:`list`) : list with names of the input files included in the
                script. If multiTask is False this list is composed by a single element
            jobnames (:py:class:`list`) : list the jobnames associated to the runs.
                If multiTask is False this list is composed by a single element
        """
        scheduler = self.run_options['scheduler']

        job = None
        if scheduler == 'direct':
            job = self.direct_runner(names,jobnames)
        elif scheduler == 'slurm':
            job = self.slurm_runner(names,jobnames)
        else:
            print('scheduler unknown')
        return job

    def submit_job(self,job):
        """
        Submit the job.

        Args:
            job : The reference to the job to be executed. If the scheduler is `direct`
                job is a list with the instance of :py:class:multiprocessing. If the
                scheduler is `slurm` job is the name of the slurm script
        """
        scheduler = self.run_options['scheduler']

        if scheduler == 'direct':
            for run in job:
                run.start()
        if scheduler == 'slurm':
            print('slurm submit to be implemented')

    def wait(self,job):
        """
        Wait the end of the job.

        Args:
            job : The reference to the job to be executed. If the scheduler is `direct`
                job is a list with the instance of :py:class:multiprocessing. If the
                scheduler is `slurm` job is the name of the slurm script

        Note:
            Actually implemented only for scheduler `direct`

        """
        verbose = self.run_options['verbose']
        IO_time = self.run_options['IO_time']
        import time

        while not self._is_terminated(job):
            if verbose:
                s = ''
                for ind,run in enumerate(job):
                    s+='run'+str(ind)+'_is_running:'+str(run.is_alive()) + '  '
                print(s)
            time.sleep(IO_time)
        if verbose : print('Job completed')

    def _is_terminated(self,job):
        """
        Check if all the runs of the jobs list have been performed.
        """
        scheduler = self.run_options['scheduler']

        if scheduler == 'direct':
            terminated = all([not run.is_alive() for run in job])
            return terminated
        if scheduler == 'slurm':
            print('is_terminated method to be implemented for slurm scheduler')
            return True

    def direct_runner(self,names,jobnames):
        """
        Define the list of Process (methods of multiprocessing) associated to the
        runs of the job. If a run can be skipped it is not included in the list of
        the job processes.

        Args:
            names (:py:class:`list`) : list with names of the input files included in the
                script. If multiTask is False this list is composed by a single element
            jobnames (:py:class:`list`) : list the jobnames associated to the runs.
                If multiTask is False this list is composed by a single element
        Return:
            :py:class:`list` : list of the :py:class:`multiprocessing` objects
            associated to the runs of the job

        """
        import multiprocessing
        def os_system_run(comm_str):
            os.system(comm_str)

        job = []
        for name,jobname in zip(names,jobnames):
            run_command = self.run_command(name,jobname)
            if run_command is not None:
                p = multiprocessing.Process(target=os_system_run, args=(run_command,))
                job.append(p)
        return job

    def slurm_runner(self,names,jobnames):
        """
        Create the slurm script associated to the runs of the job.

        Args:
            names (:py:class:`list`) : list with names of the input files included in the
                script. If multiTask is False this list is composed by a single element
            jobnames (:py:class:`list`) : list the jobnames associated to the runs.
                If multiTask is False this list is composed by a single element

        Return:
            :py:class:`string`: the name of the slurm script

        """
        cpus_per_task = self.run_options['cpus_per_task']
        ntasks = self.run_options['ntasks']
        run_dir = self.run_options.get('run_dir', '.')

        job = os.path.join(run_dir,'job_'+names[0])+'.sh'
        if os.path.isfile(job):
            print('delete slurm script %s'%job)
            os.system('rm %s'%job)

        lines_options = []
        #SBATCH options
        lines_options.append('#!/bin/bash')
        lines_options.append('#SBATCH --ntasks=%s'%ntasks)
        lines_options.append('#SBATCH --cpus_per_task = %s'%cpus_per_task)
        lines_options.append('')
        #add the srun of the taks
        lines_run = []
        for name,jobname in zip(names,jobnames):
            comm_str = self.run_command(name,jobname)
            if comm_str is not None:
                lines_run.append(comm_str)
        if len(lines_run) == 0:
            print('No tasks for srun. Slurm script not submitted')
            return None

        lines = lines_options+lines_run
        f = open(job,'w')
        f.write('\n'.join(lines))
        f.close()

        return job

    def run_command(self,name,jobname):
        """
        Define the run command used to run the computation associated to the
        input file $name.
        If the skip attribute of run_options is True the method evaluated if
        the computation can be skipped. This is done if the folder where yambo
        write the results contains at least one file 'o-*'.

        Return:
            :py:class:`string` : command that run the executable with the $name input file. If the computation can be skipped the method returns None

        """
        verbose = self.run_options['verbose']
        skip = self.run_options['skip']
        run_dir = self.run_options.get('run_dir', '.')
        can_skip = all([skip,len(self._get_output(name)) > 0])

        # check if the computation can be skipped and return the proper comm_str
        if can_skip:
            if verbose: print('Skip the computation for input',name)
            return None
        else:
            comm_str = self._get_comm_str(name,jobname)
            if verbose: print('Executing command:', self._get_comm_str(name,jobname))
            return comm_str

    def _get_comm_str(self,name,jobname):
        """
        Define the command string used to run the computation. The command string depends on the
        choice of the scheduler. For `direct` scheduler the string set the directory to the run_dir
        and use the mpi_run parameter of the calculator. For `slurm` scheduler the string set the
        srun command.

        Return:
            :py:class:`string` : command that run the executable with the $name input file

        """
        verbose = self.run_options['verbose']
        run_dir = self.run_options.get('run_dir', '.')

        scheduler = self.run_options['scheduler']
        if scheduler == 'direct':
            set_run_dir = 'cd %s; '%run_dir
            command = (set_run_dir + self._global_options['mpi_run'] + ' ' + self._global_options['executable']).strip()
        if scheduler == 'slurm':
            command = 'srun ' + self._global_options['executable']

        input_name = name + '.in'

        comm_str =  command + ' -F %s -J %s -C %s'%(input_name,jobname,name)
        return comm_str

    def _clean_run_dir(self):
        """
        Clean the run_dir before performing the computation. Delete the out_dir
        and ndb_dir folders (if found).
        """
        run_dir = self.run_options.get('run_dir', '.')
        names = self.run_options.get('names')
        jobnames = self.run_options.get('jobnames',names)
        verbose = self.run_options['verbose']

        for name,jobname in zip(names,jobnames):
            out_dir = os.path.join(run_dir,name)
            ndb_dir = os.path.join(run_dir,jobname)

            if os.path.isdir(out_dir):
                if verbose: print('delete folder:',out_dir)
                os.system('rm -r %s'%out_dir)
            if os.path.isdir(ndb_dir):
                if verbose: print('delete folder:',ndb_dir)
                os.system('rm -r %s'%ndb_dir)
