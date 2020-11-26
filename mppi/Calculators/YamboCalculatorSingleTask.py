"""
This module is the old implementation of the YamboCalculator class and will be removed
after an extensive set of tests of the new class.
"""

from .Runner import Runner
import os

class YamboCalculatorSingleTask(Runner):
    """
    Manage a single yambo calculation.The class can be used to perform computations
    using all the executables of the Yambo package. It allows to setup the number of
    omp and the mpi_run, and check that the SAVE folder is present.

    Note:
        If skip is False the class deletes the folders with the o-* files and with
        the ndb database (if found) before the run.

    Example:
     >>> code = YamboCalculator(omp=1,mpi_run='mpirun -np 4',executable='yambo',skip=True,verbose=True)
     >>> code.run(input = ..., run_dir = ...,name = ...,jobname = ...)

    The parameters of the run method are:

    Args:
        run_dir (str) : the folder in which the simulation is performed
        input : the object the contain the instance of the input file
        name (str) : the name associated to the input file (without extension).
            This string is used also as the radical of the folder in which results
            are written as well as a part of the name of the output file.
        jobname (str) : the value of the jobname. If it left to None the
            value of name is attributed to jobname by process_run.
        verbose (bool) : set the amount of information provided on terminal
        skip (bool) : if True evaluate if the computation can be skipped. This is done
            by checking if the output folder contains at least one file o-*

    When the run method is called the class runs the command:
                executable_name -F name.in -J jobname -C name


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
        If skip = False clean the run_dir.

        Note:
            If the run_dir and/or the SAVE folder do not exist an alert is
            written but the execution of the run method proceedes.

        Returns:
            :py:class:`dict`: dictionary containing the command to be passed to
            :meth:`process_run`

        """
        verbose = self.run_options['verbose']
        run_dir = self.run_options.get('run_dir', '.')
        input = self.run_options.get('input')
        name = self.run_options.get('name','yambo')
        SAVE = os.path.join(run_dir,'SAVE')

        # check if the run_dir and SAVE folder exist and write the input
        if not os.path.isdir(run_dir):
            print('Run_dir %s does not exists'%run_dir)
        elif not os.path.isdir(SAVE):
            print('SAVE folder does not exists')
        else:
            if not (input is None):
                input.write(run_dir,name+'.in')
            else :
                print('input not provided')

        # if skip = False clean the run_dir
        skip = self.run_options['skip']
        if not skip:
            self._clean_run_dir()

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

        # Set the OMP_NUM_THREADS variable in the environment
        os.environ['OMP_NUM_THREADS'] = str(self.run_options['omp'])

        # Set the run command
        if verbose: print('Run directory', run_dir)
        comm_str = 'cd ' + run_dir + '; ' + command

        # check if the computation can be skipped and run
        can_skip = all([skip,len(self._get_output_names()) > 0])
        if can_skip:
            if verbose: print('Skip the computation for input',name)
        else:
            if verbose: print('Executing command:', command)
            os.system(comm_str)

        return {'output_names': self._get_output_names(),
                'ndb_folder': jobname}

    def post_processing(self, command, output_names, ndb_folder):
        """
        Return a dictionary with the names of the o- file(s) and the folder that
        contains the ndb database.
        """
        return {'output' : output_names, 'ndb_folder' : ndb_folder}

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

    # def _get_ndb_names(self):
    #     """
    #     Look for the names of the ndb file(s) produced by yambo.
    #
    #     Return:
    #         :py:class:`list`: A list with the names, including the path, of the
    #         files ndb produced by the run
    #     """
    #     run_dir = self.run_options.get('run_dir', '.')
    #     name = self.run_options.get('name','yambo')
    #     jobname = self.run_options.get('jobname',name)
    #     ndb_dir = os.path.join(run_dir,jobname)
    #
    #     ndb_names = []
    #     if os.path.isdir(ndb_dir):
    #         for file in os.listdir(ndb_dir):
    #             if 'ndb' in file:
    #                 ndb_names.append(os.path.join(ndb_dir,file))
    #     return ndb_names

    def _clean_run_dir(self):
        """
        Clean the run_dir before performing the computation. Delete the out_dir
        and ndb_dir (if found).
        """
        verbose = self.run_options['verbose']
        run_dir = self.run_options.get('run_dir', '.')
        name = self.run_options.get('name','yambo')
        jobname = self.run_options.get('jobname',name)
        out_dir = os.path.join(run_dir,name)
        ndb_dir = os.path.join(run_dir,jobname)

        if os.path.isdir(out_dir):
            if verbose: print('delete folder:',out_dir)
            os.system('rm -r %s'%out_dir)
        if os.path.isdir(ndb_dir):
            if verbose: print('delete folder:',ndb_dir)
            os.system('rm -r %s'%ndb_dir)