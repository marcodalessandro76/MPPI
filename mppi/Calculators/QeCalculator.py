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
       executable (:py:class:`string`) : set the executable (pw.x, ph.x, ..) of the QuantumESPRESSO package
       multiTask  (:py:class:`bool`) : if true a single run_script is built and all the computations are performed in parallel,
            otherwise an independent script is built for each elements of inputs and the computations are performed sequentially
       mpi_run (:py:class:`string`) : command for the execution of mpirun, e.g. 'mpirun -np 4'
       cpus_per_task (:py:class:`int`) : set the `cpus_per_task` variable for the slurm script
       ntasks (:py:class:`int`) : `set the ntasks` variable for the slurm script
       skip (:py:class:`bool`) : if True evaluate if one (or many) computations can be skipped.
           This is done by checking if the file $name.xml is present in the prefix folder,
           for each name in names
       verbose (:py:class:`bool`) : set the amount of information provided on terminal
       IO_time (int) : time step (in second) used by the wait method to check that the job is completed

    Example:
     >>> code = calculator(omp=1,mpi_run='mpirun -np 4',skip=True,verbose=True,scheduler='direct')
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
                 omp = os.environ.get('OMP_NUM_THREADS', 1), executable = 'pw.x',
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
        verbose = self.run_options['verbose']

        # Create the run_dir and write the input file
        self._ensure_run_directory()
        if inputs is not None:
            for input,name in zip(inputs,names):
                input.write(os.path.join(run_dir,name)+'.in')
        else:
            print('input list not provided')

        # if skip = False clean the run_dir
        skip = self.run_options['skip']
        if not skip:
            self._clean_run_dir()

        # Copy the source folder in the run_dir
        source_dir = self.run_options.get('source_dir')
        if source_dir is not None:
            self._copy_source_dir(source_dir)

        return {}

    def process_run(self):
        """
        Method associated to the running of the executable. The method prepare the
        job script in a way that depend on the chosen scheduler, then submit the job
        and wait the end of the computation before passing to the :meth:`post_processing` method.

        Note:
            The wait method is actually not implemented for the slurm scheduler

        Returns:
           :py:class:`dict`: The dictionary `results_file`
             values to be passed to the :meth:`post_processing` method

        """
        # Set the OMP_NUM_THREADS variable in the environment
        os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])
        names = self.run_options.get('names')
        multiTask = self.run_options.get('multiTask')

        if multiTask:
            job = self.build_run_script(names)
            self.submit_job(job)
            self.wait(job)
        else:
            for name in names:
                job = self.build_run_script([name])
                self.submit_job(job)
                self.wait(job)

        return {}

    def post_processing(self):
        """
        Return a list with the names, including the path, of the data-file-schema.xml
        files built by pw, for each element of inputs. If a file is absent the method returns
        None in the associated element of the list, making easy to understand if a specific
        computation has been correctly performed.

        Return:
            :py:class:`list` : list with the names of the xml files (if the file exists)
                otherwise return None.

        """
        run_dir = self.run_options.get('run_dir', '.')
        inputs = self.run_options['inputs']
        results = []
        for input in inputs:
            prefix = input['control']['prefix'].strip("'")
            prefix += '.save'
            result = os.path.join(run_dir,prefix,'data-file-schema.xml')
            if os.path.isfile(result):
                results.append(result)
            else:
                results.append(None)

        return results

    def build_run_script(self,names):
        """
        Create the run script that is executed by the :meth:`submit_job` method.
        The script depends on the scheduler adopted, and specific methods for
        `direct` and `slurm` scheduler are implemented.

        Args:
            names (:py:class:`list`) : list with names of the input files included in the
                script. If multiTask is False this list is composed by a single element

        """
        scheduler = self.run_options['scheduler']

        job = None
        if scheduler == 'direct':
            job = self.direct_runner(names)
        elif scheduler == 'slurm':
            job = self.slurm_runner(names)
        else:
            print('scheduler unknown')
        return job

    def submit_job(self,job):
        """
        Submit the job.
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

    def direct_runner(self,names):
        """
        Define the list of Process (methods of multiprocessing) associated to the
        runs of the job. If a run can be skipped it is not included in the list of
        the job processes.

        Args:
            names (:py:class:`list`) : list with names of the input files included in the
                script. If multiTask is False this list is composed by a single element

        Return:
            :py:class:`list` : list of the :py:class:`multiprocessing` objects
            associated to the runs of the job

        """
        import multiprocessing
        def os_system_run(comm_str):
            os.system(comm_str)

        job = []
        for name in names:
            run_command = self.run_command(name)
            if run_command is not None:
                p = multiprocessing.Process(target=os_system_run, args=(run_command,))
                job.append(p)
        return job

    def slurm_runner(self,names):
        """
        Create the slurm script associated to the runs of the job.

        Args:
            names (:py:class:`list`) : list with names of the input files included in the
                script. If multiTask is False this list is composed by a single element

        Return:
            :py:class:`string`: the name of the slurm script

        """
        cpus_per_task = self.run_options['cpus_per_task']
        ntasks = self.run_options['ntasks']
        run_dir = self.run_options.get('run_dir', '.')

        job_name = os.path.join(run_dir,'job_'+names[0])+'.sh'
        if os.path.isfile(job_name):
            print('delete slurm script %s'%job_name)
            os.system('rm %s'%job_name)

        lines_options = []
        #SBATCH options
        lines_options.append('#!/bin/bash')
        lines_options.append('#SBATCH --ntasks=%s'%ntasks)
        lines_options.append('#SBATCH --cpus_per_task = %s'%cpus_per_task)
        lines_options.append('')
        #add the srun of the taks
        lines_run = []
        for name in names:
            comm_str = self.run_command(name)
            if comm_str is not None:
                lines_run.append(comm_str)
        if len(lines_run) == 0:
            print('No tasks for srun. Slurm script not submitted')
            return None

        lines = lines_options+lines_run
        f = open(job_name,'w')
        f.write('\n'.join(lines))
        f.close()

        return job_name

    def run_command(self,name):
        """
        Define the run command used to run the computation associated to the
        input file $name.
        If the skip attribute of run_options is True the method evaluated if
        the computation can be skipped. This is done by  checking if the file
        $name.xml is already present in the path run_dir/prefix.save

        Return:
            :py:class:`string` : command that run the executable with the $name input file.
                If the computation can be skipped the method returns None

        """
        verbose = self.run_options['verbose']
        skip = self.run_options['skip']
        run_dir = self.run_options.get('run_dir', '.')
        skipfile = os.path.join(run_dir,name)+'.xml'

        # check if the computation can be skipped and return the proper comm_str
        can_skip = all([skip,os.path.isfile(skipfile)])
        if can_skip:
            if verbose: print('Skip the computation for input',name)
            return None
        else:
            comm_str = self._get_comm_str(name)
            if verbose: print('Executing command:', self._get_comm_str(name))
            return comm_str

    def _get_comm_str(self,name):
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
        output_name = name + '.log'

        comm_str =  command + ' -inp %s > %s'%(input_name,output_name)
        return comm_str

    # def _get_results(self):
    #     """
    #     Return a list with the name, including the path, of the data-file-schema.xml
    #     file built by pw.
    #     """
    #     run_dir = self.run_options.get('run_dir', '.')
    #     inputs = self.run_options['inputs']
    #     results = []
    #     for input in inputs:
    #         prefix = input['control']['prefix'].strip("'")
    #         prefix += '.save'
    #         results.append(os.path.join(run_dir,prefix,'data-file-schema.xml'))
    #     return results

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
        verbose = self.run_options['verbose']
        run_dir = self.run_options.get('run_dir', '.')
        inputs = self.run_options['inputs']

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
        the $name.xml file and the folder run_dir/prefix.save associated to all the
        inputs and names.
        """
        run_dir = self.run_options.get('run_dir', '.')
        names = self.run_options.get('names','default')
        inputs = self.run_options.get('inputs')
        verbose = self.run_options['verbose']

        for input,name in zip(inputs,names):
            logfile = os.path.join(run_dir,name)+'.log'
            xmlfile = os.path.join(run_dir,name)+'.xml'
            prefix = input['control']['prefix'].strip("'")
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
