"""
This module defines a class to perform a calculations using QuantumESPRESSO.
The module is (deeply) inspired from the SystemCalculator class of BigDFT
"""

import os
from qepppy import qe

class QeCalculator():
    """
    Manage a single calculation. Setup the number of omp and mpi.
    Prepare the folder in which the computation is run. Peform the computation
    and apply a post processing function to extract the results

    """
    def __init__(self,omp=1,mpi_run='mpirun -np 4',executable='',skip=False, verbose=True):
        self.omp=omp
        self.mpi_run=mpi_run
        self.executable=executable
        """
        Choose the executable called by run
        """

        self.skip=skip
        """
        Set if skip the run if the file self.xml_file is found. Default is False
        """

        self.verbose=verbose
        """
        Set the verbosity option
        """

        self.xml_file = ''
        """
        The xml file (with its complete path) used both for post_processing and to
        establish if the run can be skipped.
        """

        self.command = ('OMP_NUM_THREADS='+ str(self.omp) + ' ' + self.mpi_run + ' ' + executable).strip()
        print('Initialize a qe calculator with command %s' %self.command)

    def pre_processing(self,**kwargs):
        """
        Check if the run_dir folder exists and write the input file.
        """
        run_dir=kwargs['run_dir']
        if not os.path.isdir(run_dir):
            print('Run_dir %s does not exists'%run_dir)
        else:
            input=kwargs['input']
            name=kwargs['name'] + '.in'
            if not (input is None):
                input.write(run_dir + '/' + name)
            else :
                print('input not provided')

    def process_run(self,**kwargs):
        """
        Set the proper dir and run the computation. If skip=True run the computations
        only if the .log file is not present in the folder
        The skip could be imposed on the .xml file that has the same name of the prefix
        """
        string = 'cd %s ; '%kwargs['run_dir']
        input=kwargs['name'] + '.in'
        output=kwargs['name'] + '.log'
        string +=  self.command + ' -inp %s > %s'%(input,output)
        if self.skip:
            if os.path.isfile(kwargs['run_dir']+'/'+output):
                if self.verbose : print('skip the computation for : '+output)
            else:
                if self.verbose : print('execute : '+string)
                os.system(string)
        else:
            if self.verbose : print('execute : '+string)
            os.system(string)

    def post_processing(self,**kwargs):
        """
        Use the pw_out() method os qepppy.qe to parse the xml file
        that contains the results of the computation
        """
        input = kwargs['input']
        results = None
        if 'prefix' in input.control:
            prefix = input.control['prefix'].strip("'")
            self.xml_file = kwargs['run_dir'] + '/' + prefix + '.save/data-file-schema.xml'
            if self.verbose : print('parse file : '+self.xml_file)
            results= qe.pw_out(xml=self.xml_file)
        else:
            print('.save folder not provided. Cannot read xml file')
        return results

    def run(self,run_dir='run',input=None,name='test',post_processing=True):
        """
        Prepare the run, perform the computation and apply the post_processing
        function to extract the results
        Args:
            run_dir (str) : the folder in which the simulation is performed
            input : the object the contain the instance of the input file
            name (str) : the name associated to the input file (without extension)
        """
        self.pre_processing(run_dir=run_dir,input=input,name=name)
        self.process_run(run_dir=run_dir,name=name)
        results = None
        if post_processing:
            results = self.post_processing(run_dir=run_dir,input=input)
        return results
