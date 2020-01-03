"""
This module manages parallel calculations with QuantumESPRESSO.
Both a scheduler like slurm or the python multiprocessing package can be used.
"""

from .Runner import Runner
import os

class QeParallelCalculator(Runner):
    """
    Manage multiple QuantumESPRESSO calculations performed in parallel. Computations
    are managed by a scheduler that, in the actual implementation of the class, can
    be 'direct' or 'slurm'.

    Example:
     >>> code = calculator(omp=1,mpi_run='mpirun -np 4',skip=True,verbose=True)
     >>> code.run(inputs = ..., run_dir = ...,names = ..., source_dir = ...)

    Args:
        run_dir (str) : the folder in which the simulation is performed
        inputs (:py:class:`list`) : list with the instances of the PwInput class
        that define the input objects
        names (:py:class:`list`) : list with the names associated to the input files,
            given in the same order of the inputs list.
            Usually you can set the name equal to the prefix of the input object so
            the name of the input file and the prefix folder built by QuantumESPRESSO
            are equal
        source_dir (str) : location of the scf source folder for a nscf computation.
        If present the class copies this folder in the run_dir with the name $prefix.save.
        verbose (bool) : set the amount of information provided on terminal
        skip (bool) : if True evaluate if one (or many) computations can be skipped.
            This is done by checking if the file $name.xml is present in the prefix folder,
            for each name in names

    """

    def __init__(self,
                 omp = os.environ.get('OMP_NUM_THREADS', 1), mpi_run = 'mpirun -np 2',
                 isParallel = True, scheduler = 'direct', executable = 'pw.x',
                 skip = False, verbose = True):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=omp, mpi_run=mpi_run,
                        isParallel=isParallel, scheduler=scheduler,
                        executable=executable,
                        skip=skip, verbose=verbose)
        print('Initialize a parallel QuantumESPRESSO calculator with scheduler %s' %
            self._global_options['scheduler'])

    def pre_processing(self):
        """
        Process local run dictionary to create the run directory and input files.
        If skip = False clean the run_dir.
        If the 'source_dir' key is passed to the run method copy the source folder
        in the run_dir with the name $prefix. This procedure is performed
        after the deletion run_dir/prefix.save since otherwise the copy of the
        source_dir is deleted.

        Returns:
            :py:class:`dict`: dictionary containing the command to be passed to
            :meth:`process_run`

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
        Routine associated to the running of the executable.

        Returns:
           :py:class:`dict`: The dictionary `results_file`
             values to be passed to `post_processing` function
        """

        # Set the OMP_NUM_THREADS variable in the environment
        os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])

        # Prepare the run `script`
        job = self.build_run_script()
        # Submit the job
        self.submit_job(job)
        # Wait the end of the job
        self.wait(job)

        return {'results_file': self._get_results_file()}

    def post_processing(self, results_file):
        """
        Return results_name (if the file exists) otherwise return None.
        This check allows us to understand if the computation has been correctly
        performed.

        Return:
            ......

        """
        for ind,file in enumerate(results_file):
            if not os.path.isfile(file):
                result_file[ind] = None
        return results_file

    def build_run_script(self):
        """
        Create the run script
        """
        scheduler = self.run_options['scheduler']
        print('scheduler',scheduler)

        job = None
        if scheduler == 'direct':
            job = self.direct_runner()
        elif scheduler == 'slurm':
            print('to be implemented')
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
            print('to be implemented')

    def wait(self,job):
        """
        Wait the end of the job.
        """
        verbose = self.run_options['verbose']
        import time

        while not self._is_terminated(job):
            if verbose:
                s = ''
                for ind,run in enumerate(job):
                    s+='run'+str(ind)+'_is_running:'+str(run.is_alive()) + '  '
                print(s)
            time.sleep(5)
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
            print('is_direct to be implemented')
            return True

    def direct_runner(self):
        """
        Define the list of Process (methods of multiprocessing) associated to the
        runs of the job.

        Return:
            (:py:class:list) : list of the :py:class:multiprocessing objects associated
                to the runs of the job

        """
        import multiprocessing
        names = self.run_options.get('names')
        def os_system_run(comm_str):
            os.system(comm_str)

        job = []
        for name in names:
            p = multiprocessing.Process(target=os_system_run, args=(self.run_command(name),))
            job.append(p)
        return job

    def slurm_runner(self):
        """
        Create the slurm script associated to the runs of the job.
        To be implemented.

        Return:
            (:py:class:string) : name of the slurm script

        """
        job = 'slurm_job.sh'
        return job

    def run_command(self,name):
        """
        Define the run command used to the run the computation associated to the
        input file $name.
        If the skip attribute of run_options is True the method evaluated if
        the computation can be skipped. This is done by  checking if the file
        $name.xml is already present in the path run_dir/prefix.save

        Return:
            (string) : command that run the executable with the $name input file

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
            comm_str = 'cd ' + run_dir + '; ' + self._get_comm_str(name)
            if verbose: print('Executing command:', self._get_comm_str(name))
            return comm_str

    def _get_comm_str(self,name):
        """
        Define the command string used to run the computation.
        TODO : check the command if scheduler is slurm.

        """
        scheduler = self.run_options['scheduler']
        if scheduler == 'direct':
            command = (self._global_options['mpi_run'] + ' ' + self._global_options['executable']).strip()
        if scheduler == 'slurm':
            command = self._global_options['executable']

        input_name = name + '.in'
        output_name = name + '.log'

        comm_str =  command + ' -inp %s > %s'%(input_name,output_name)
        return comm_str

    def _get_results_file(self):
        """
        Return a list with the name, including the path, of the data-file-schema.xml
        file built by pw.
        """
        run_dir = self.run_options.get('run_dir', '.')
        inputs = self.run_options['inputs']
        results = []
        for input in inputs:
            prefix = input['control']['prefix'].strip("'")
            prefix += '.save'
            results.append(os.path.join(run_dir,prefix,'data-file-schema.xml'))
        return results

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
        inputs = self.run_options['input']

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
