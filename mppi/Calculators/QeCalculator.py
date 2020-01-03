"""
A class to perform a calculations using QuantumESPRESSO.
"""

from .Runner import Runner
import os

class QeCalculator(Runner):
    """
    Manage a single  QuantumESPRESSO calculation. Setup the number of omp and mpi,
    prepare the folder in which the computation is run and perform the computation.

    Note:
        If skip is False the class delete the name.log and name.xml files and the
        run_dir/prefix.save folder before execute the run.

    Example:
     >>> code = calculator(omp=1,mpi_run='mpirun -np 4',skip=True,verbose=True)
     >>> code.run(input = ..., run_dir = ...,name = ..., source_dir = ...)

    Args:
        run_dir (str) : the folder in which the simulation is performed
        input (PwInput) : the object the contain the instance of the input file
        name (str) : the name associated to the input file (without extension)
            Usually you can set the name equal to the prefix of the input object so
            the name of the input file and the prefix folder built by QuantumESPRESSO
            are equal
        source_dir (str) : location of the scf source folder for a nscf computation.
        If present the class copies this folder in the run_dir with the name $prefix.save.
        verbose (bool) : set the amount of information provided on terminal
        skip (bool) : if True evaluate if the computation can be skipped. This is done
            by checking if the file $name.xml is present in the prefix folder
    """

    def __init__(self,
                 omp=os.environ.get('OMP_NUM_THREADS', 1),
                 mpi_run='mpirun -np 4',executable='pw.x',
                 skip=False, verbose=True):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=omp, mpi_run=mpi_run, executable=executable,
                        skip=skip, verbose=verbose)
        self.command = (self._global_options['mpi_run'] + ' ' + executable).strip()
        print('Initialize a QuantumESPRESSO calculator with OMP_NUM_THREADS=%s and command %s' %
            (self._global_options['omp'], self.command))

    def pre_processing(self):
        """
        Process local run dictionary to create the run directory and input file.
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
        input = self.run_options.get('input')
        name = self.run_options.get('name','default')
        verbose = self.run_options['verbose']

        # Create the run_dir and write the input file
        self._ensure_run_directory()
        if input is not None:
            input.write(os.path.join(run_dir,name)+'.in')
        else:
            print('input not provided')

        # if skip = False clean the run_dir
        skip = self.run_options['skip']
        if not skip:
            self._clean_run_dir()

        # Copy the source folder in the run_dir
        source_dir = self.run_options.get('source_dir')
        if source_dir is not None:
            self._copy_source_dir(source_dir)

        return {'command': self._get_command()}

    def process_run(self,command):
        """
        Launch the code.

        Routine associated to the running of the executable.
        If the skip attribute of run_options is True the method evaluated if
        the computation can be skipped. This is done by  checking if the file
        $name.xml is already present in the path run_dir/prefix.save

        Args:
           command (str): the command as it is set by the ``pre_processing``
             method.

        Returns:
           :py:class:`dict`: The dictionary `results`
             values to be passed to `post_processing` function
        """
        verbose = self.run_options['verbose']
        skip = self.run_options['skip']
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        skipfile = os.path.join(run_dir,name)+'.xml'

        # Set the OMP_NUM_THREADS variable in the environment
        os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])

        # Set the run the command
        if verbose: print('Run directory', run_dir)
        comm_str = 'cd ' + run_dir + '; ' + command

        # check if the computation can be skipped and run
        can_skip = all([skip,os.path.isfile(skipfile)])
        if can_skip:
            if verbose: print('Skip the computation for input',name)
        else:
            if verbose: print('Executing command:', command)
            os.system(comm_str)

        return {'result_file': self._get_result_file()}

    def post_processing(self, command, result_file):
        """
        Return results_name (if the file exists) otherwise return None.
        This check allows us to understand if the computation has been correctly
        performed.
        """
        if os.path.isfile(result_file):
            return result_file
        else:
            return None

    def _get_command(self):
        name = self.run_options.get('name','default')
        input_name = name + '.in'
        output_name = name + '.log'

        comm_str =  self.command + ' -inp %s > %s'%(input_name,output_name)
        return comm_str

    def _get_result_file(self):
        """
        Return the name, including the path, of the data-file-schema.xml
        file built by pw.
        """
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options['input']
        prefix = input['control']['prefix'].strip("'")
        prefix += '.save'
        return os.path.join(run_dir,prefix,'data-file-schema.xml')

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
        the name $prefix

        Args:
            source_dir: the name of the source_dir including its relative path.
            A source_dir outer respect to the actual run_dir of the instance of
            QeCalculator can be used.
        """
        from shutil import copytree
        verbose = self.run_options['verbose']
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options['input']
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
        the $name.xml file and the folder run_dir/prefix.save.
        """
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        input = self.run_options.get('input')
        prefix = input['control']['prefix'].strip("'")
        verbose = self.run_options['verbose']

        logfile = os.path.join(run_dir,name)+'.log'
        xmlfile = os.path.join(run_dir,name)+'.xml'
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
