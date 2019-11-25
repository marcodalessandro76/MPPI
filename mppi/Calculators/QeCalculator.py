"""
A class to perform a calculations using QuantumESPRESSO.
"""

from .Runner import Runner
#from qepppy import qe
from mppi.Parsers import PwParser
import os

class QeCalculator(Runner):
    """
    Manage a single  QuantumESPRESSO calculation. Setup the number of omp and mpi.
    Prepare the folder in which the computation is run. Peform the computation
    and apply a post processing function to extract the results.

    """

    def __init__(self,
                 omp=os.environ.get('OMP_NUM_THREADS', '1'),
                 mpi_run='mpirun -np 4',executable='pw.x',
                 skip=False, verbose=True):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=str(omp), mpi_run=mpi_run, executable=executable,
                        skip=skip, verbose=verbose)
        self.command = (self._global_options['mpi_run'] + ' ' + executable).strip()
        print('Initialize a QuantumESPRESSO calculator with OMP_NUM_THREADS=%s and command %s' %
            (self._global_options['omp'], self.command))

    def pre_processing(self):
        """
        Process local run dictionary to create the run directory and input file.
        If the 'source_dir' key is passed to the run method copy the source folder in
        the run_dir with the name 'name'.

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
            print('copy source_dir')
            # to be implemented....

        return {'command': self._get_command()}

    def process_run(self,command):
        """Launch the code.

        Routine associated to the running of the executable. It evaluates if the
        computation can be skipped. For the QuantumESPRESSO computation we use the
        $prefix.xml file in the run_dir which coincides with the xml file that collect
        the results.

        The amount of information provided on terminal is set by the verbose key
        of the run.

        Arguments:
           command (str): the command as it is set by the ``pre_processing``
             method.

        Returns:
           :py:class:`dict`: The dictionary `results`
             values to be passed to `post_processing` function
        """

        verbose = self.run_options['verbose']
        skip = self.run_options['skip']
        run_dir = self.run_options['run_dir']
        name = self.run_options['name']
        skipfile = self._get_results_name()

        # Set the number of omp threads only if the variable is not present
        # in the environment
        if 'OMP_NUM_THREADS' not in os.environ:
            os.environ['OMP_NUM_THREADS'] = self.run_options['omp']

        # Run the command
        if verbose:
            if run_dir != '.':
                print('Run directory', run_dir)
        comm_str = 'cd ' + run_dir + '; ' + command
        if skip:
            if os.path.isfile(skipfile):
                if verbose: print('Skip the computation for input ',name)
            else:
                if verbose: print('Executing command: ', command)
                os.system(comm_str)
        else:
            if verbose: print('Executing command: ', command)
            os.system(comm_str)

        return {'results_name': self._get_results_name()}

    def post_processing(self, command, results_name):
        """
        Parse the xml file that contains the results of the computation.
        Check the existence and the result_name file and return an instance
        of Qe Parser?.

        Returns:
            PwParser: Instance of the PwParser class associated to the run
            which has been just performed. If the run failed for some reasons
            and the result_name file does not exists or it cannot be parsed it
            returns `None`.
        """

        verbose = self.run_options['verbose']

        run_dir = self.run_options['run_dir']
        input = self.run_options['input']
        prefix = input['control']['prefix'].strip("'")
        results = PwParser(prefix,path=run_dir,verbose=verbose)
        return results

        # if os.path.isfile(results_name):
        #     if verbose: print('parse file : ',command)
        #     #results = 'perform parsing'
        #     results= qe.pw_out(xml=results_name)
        # else:
        #     if verbose: print('ERROR: The file %s does not exists'%command)
        #     results = None
        # return results

    def _get_command(self):
        name = self.run_options.get('name','default')
        input = name + '.in'
        output = name + '.log'

        comm_str =  self.command + ' -inp %s > %s'%(input,output)
        return comm_str

    def _get_results_name(self):
        run_dir = self.run_options['run_dir']
        input = self.run_options['input']
        prefix = input['control']['prefix']
        # remove the first and last character that are in the qe input file format
        prefix = prefix[1:-1]
        results_name = prefix+'.xml'
        return os.path.join(run_dir,results_name)

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

        # the run_dir attribute has been removed...
        #self.run_dir = run_dir
        # """Run directory.
        # str: the directory where the inputfile has been copied to.
        #       Might be useful to associate to each of the calculation of a
        #       given run a different directory. Note that this is different than
        #       setting the ``outdir`` or the ``name`` arguments at it refers
        #       to the directory of the inputfile.
        # Note:
        #     This is not a global property of the calculator, as the same
        #     calculator instance can be used for various workflows.
        # """
