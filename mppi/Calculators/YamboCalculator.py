"""
A class to perform a calculations using Yambo.
"""

from .Runner import Runner
import os
#from .YamboParser import dict_parser, AttributeDict

class YamboCalculator(Runner):
    """
    Manage a single yambo calculation: setup the number of omp and mpi, check
    that the SAVE folder is present and perform the computation.
    Lastly perform a parsing of the results.

    In this cases, the calculator assumes that the run_dir exists and that
    it contains the SAVE folder built with the p2y postprocessing of a
    QuantumESPRESSO computation.

    Note:
        If skip is False the class delete the folder with the o-* files (if found)
        before the run. Is this choice correct or is it better to keep the folder?

    """
    def __init__(self,
                 omp=os.environ.get('OMP_NUM_THREADS', 1),
                 mpi_run='mpirun -np 4',executable='yambo',
                 skip=False, verbose=True):
        # Use the initialization from the Runner class (all options inside _global_options)
        Runner.__init__(self, omp=omp, mpi_run=mpi_run, executable=executable,
                        skip=skip, verbose=verbose)
        self.command = (self._global_options['mpi_run'] + ' ' + executable).strip()
        print('Initialize a Yambo calculator with OMP_NUM_THREADS=%s and command %s' %
            (self._global_options['omp'], self.command))

    def pre_processing(self):
        """
        Process local run dictionary. Check that the run_dir exists and that it
        contains the SAVE folder. Check that the input object has been provided
        in the run parameters and write it on disk.

        Note:
            If the run_dir and/or the SAVE folder do not exist an alert is
            written but the execution of the run method proceedes.

        Returns:
            :py:class:`dict`: dictionary containing the command to be passed to
            :meth:`process_run`

        """
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options.get('input')
        name = self.run_options.get('name','yambo')
        SAVE = os.path.join(run_dir,'SAVE')

        if not os.path.isdir(run_dir):
            print('Run_dir %s does not exists'%run_dir)
        elif not os.path.isdir(SAVE):
            print('SAVE folder does not exists')
        else:
            if not (input is None):
                input.write(run_dir,name+'.in')
            else :
                print('input not provided')

        return {'command': self._get_command()}

    def process_run(self,command):
        """Launch the code.

        Routine associated to the running of the executable.
        If the skip attribute o run_options is True the method evaluated if
        the computation can be skipped. This is done if the the folder where yambo
        write the results contains at least one file 'o-*'.

        The amount of information provided on terminal is set by the verbose key
        of the run.

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
        name = self.run_options.get('name','yambo')
        jobname = self.run_options.get('jobname',name)
        out_dir = os.path.join(run_dir,name)
        ndb_dir = os.path.join(run_dir,jobname)

        # Check if (at least) a file o-* is found in the out_dir
        outfile_found = False
        num_outfiles = len(self._get_output_names())
        if num_outfiles != 0: outfile_found = True

        # Set the OMP_NUM_THREADS variable in the environment
        os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])

        # Run the command
        if verbose: print('Run directory', run_dir)
        comm_str = 'cd ' + run_dir + '; ' + command
        if skip:
            if outfile_found:
                if verbose: print('Skip the computation for input ',name)
            else:
                if verbose: print('Executing command:', command)
                os.system(comm_str)
        else:
            # delete the out_dir and run_dir (if found) before running the code
            if os.path.isdir(out_dir):
                if verbose: print('delete folder:',out_dir)
                os.system('rm -r %s'%out_dir)
            if os.path.isdir(ndb_dir):
                if verbose: print('delete folder:',ndb_dir)
                os.system('rm -r %s'%ndb_dir)
            if verbose: print('Executing command:', command)
            os.system(comm_str)

        return {'output_names': self._get_output_names(),
                'ndb_names': self._get_ndb_names() }

    def post_processing(self, command, output_names, ndb_names):
        """
        In the actual implementation return a dictionary with the names
        of the o- file(s) and ndb database.
        """
        return {'output' : output_names, 'ndb' : ndb_names}

    def _get_command(self):
        name = self.run_options.get('name','yambo')
        jobname = self.run_options.get('jobname',name)
        input = name + '.in'
        comm_str =  self.command + ' -F %s -J %s -C %s'%(input,jobname,name)
        return comm_str

    def _get_output_names(self):
        """
        Look for the names of the 'o-' file(s) produced by yambo.

        Return:
            :py:class:`list`: A list with the names, including the path, of the
            files 'o-*' produced by the run
        """
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','yambo')
        out_dir = os.path.join(run_dir,name)

        output_names = []
        if os.path.isdir(out_dir):
            for file in os.listdir(out_dir):
                if 'o-' in file:
                    output_names.append(os.path.join(out_dir,file))

        return output_names

    def _get_ndb_names(self):
        """
        Look for the names of the ndb file(s) produced by yambo.

        Note:
            This method needs to be modified to extract only the database
            that we need.

        Return:
            :py:class:`list`: A list with the names, including the path, of the
            files ndb produced by the run

        """
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','yambo')
        jobname = self.run_options.get('jobname',name)
        ndb_dir = os.path.join(run_dir,jobname)

        ndb_names = []
        if os.path.isdir(ndb_dir):
            for file in os.listdir(ndb_dir):
                if 'ndb' in file:
                    ndb_names.append(os.path.join(ndb_dir,file))

        return ndb_names


    # def run(self,run_dir='run',input=None,name='test',jobname=None,post_processing=True):
    #     """
    #     Prepare the run, perform the computation and apply the post_processing
    #     function to extract the results. It returns the object built by YamboParser.
    #     Args:
    #         run_dir (str) : the folder in which the simulation is performed
    #         input : the object the contain the instance of the input file
    #         name (str) : the name associated to the input file (without extension).
    #         This string is used also as the radical of the folder in which results
    #         are written as well as a part of the name of the output file.
    #         jobname (str) : the value of the jobname. If it left to None the
    #         value of name is attributed to jobname by process_run.
    #     """
    #     self.pre_processing(run_dir=run_dir,input=input,name=name)
    #     self.process_run(run_dir=run_dir,name=name,jobname=jobname)
    #     results = None
    #     if post_processing:
    #         results = self.post_processing(run_dir=run_dir,name=name,jobname=jobname)
    #     return results
