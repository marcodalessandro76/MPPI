"""
A class to perform a calculations using QuantumESPRESSO.
"""

from .Runner import Runner
#from mppi.Parsers import PwParser
import os

class QeCalculator(Runner):
    """
    Manage a single  QuantumESPRESSO calculation. Setup the number of omp and mpi.
    Prepare the folder in which the computation is run. Peform the computation
    and apply a post processing function to extract the results.

    Note:
        The argument self.run_dir has been removed. The value of the run_dir is
        extracted directly from the run_options. Evaluate if it can be better to
        introduce again and also to add the attribute self.prefix.

    Example:
     >>> code = calculator(omp=1,mpi_run='mpirun -np 4',skip=True)
     >>> code.run(input = ..., run_dir = ...,name = ...)

    Args:
        run_dir (str) : the folder in which the simulation is performed
        input : the object the contain the instance of the input file
        name (str) : the name associated to the input file (without extension).
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
        If the 'source_dir' key is passed to the run method copy the source folder
        in the run_dir with the name $prefix.

        Returns:
            :py:class:`dict`: dictionary containing the command to be passed to
            :meth:`process_run`

        """
        self._ensure_run_directory()

        # Create the input file
        run_dir = self.run_options.get('run_dir', '.')
        inp = self.run_options.get('input')
        name = self.run_options.get('name','default')
        if inp is not None:
            inp.write(os.path.join(run_dir,name)+'.in')
        else:
            print('input not provided')

        # Copy the source folder in the run_dir
        source_dir = self.run_options.get('source_dir')
        if source_dir is not None:
            self._copy_source_dir(source_dir)

        return {'command': self._get_command()}

    def process_run(self,command):
        """Launch the code.

        Routine associated to the running of the executable.
        If the skip attribute o run_options is True the method evaluated if
        the computation can be skipped. This is done by if the
        '$file self.run_options['name'].xml'
        is already present in the run_dir.

        The amount of information provided on terminal is set by the verbose key
        of the run.

        Args:
           command (str): the command as it is set by the ``pre_processing``
             method.

        Returns:
           :py:class:`dict`: The dictionary `results`
             values to be passed to `post_processing` function

        Note:
            It could happens that the '$file self.run_options['name'].xml' is present
            but the 'data-file-schema.xml' is not, so the computations is skipped but
            the parsing is not performed. Maybe is better to directly use the
            data-file-schema to establish if the run can be skipped.
        """

        verbose = self.run_options['verbose']
        skip = self.run_options['skip']
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','default')
        skipfile = os.path.join(run_dir,name)+'.xml'

        # Set the OMP_NUM_THREADS variable in the environment
        os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])

        # Run the command
        if verbose: print('Run directory', run_dir)
        comm_str = 'cd ' + run_dir + '; ' + command
        if skip:
            if os.path.isfile(skipfile):
                if verbose: print('Skip the computation for input ',name)
            else:
                if verbose: print('Executing command:', command)
                os.system(comm_str)
        else:
            if verbose: print('Executing command:', command)
            os.system(comm_str)

        return {'results_name': self._get_results()}

    def post_processing(self, command, results_name):
        """
        Parse the xml file that contains the results of the computation.

        Returns:
            PwParser: Instance of the PwParser class associated to the run
            which has been just performed. If the run failed for some reasons
            and the result_name file does not exists or it cannot be parsed it
            the attribute data of the PwParser object is set to None.
        """
        # version used if the parse is performed later....
        return {'xml_data': results_name}
        #results = PwParser(results_name,verbose=self.run_options['verbose'])
        #return results

    def _get_command(self):
        name = self.run_options.get('name','default')
        input = name + '.in'
        output = name + '.log'

        comm_str =  self.command + ' -inp %s > %s'%(input,output)
        return comm_str

    def _get_results(self):
        """
        Return the name, including the path, of the data-file-schema.xml
        file build by pw
        """
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options['input']
        prefix = input['control']['prefix'].strip("'")
        prefix += '.save'
        return os.path.join(run_dir,prefix,'data-file-schema.xml')

    def _ensure_run_directory(self):
        from mppi.Utilities import Futile_utils as f
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
