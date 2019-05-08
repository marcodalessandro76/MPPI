"""
A class to perform a calculations using Yambo.
The class is (deeply) inspired from the SystemCalculator class of BigDFT.
"""

import os
from .YamboParser import dict_parser, AttributeDict

class YamboCalculator():
    """
    Manage a single yambo calculation: setup the number of omp and mpi, apply a
    pre processing to check that the SAVE folder is present and peform the
    computation. Lastyl, it apply a post processing function to extract the
    results. The post processing is defined in the module YamboParser.

    The Class assumes that the run_dir where yambo is executed exists and that
    it contains the SAVE folder with the p2y postprocessing of a QE computation.
    The setup of this folder is managed at the level of Dataset.
    """
    def __init__(self,omp=1,mpi_run='mpirun -np 4',executable='yambo',\
                 suffix ='',skip=False, verbose=True):
        self.omp=omp
        self.mpi_run=mpi_run
        self.executable=executable
        """
        Choose the executable called by run.
        """

        self.skip=skip
        """
        Set if skip the run if the file self.output is found. Default is False.
        """

        self.verbose=verbose
        """
        Set the verbosity option.
        """

        self.output = None
        """
        The o-*.hf or o-*.qp file (with its complete path) used both for
        post_processing and to establish if the run can be skipped.
        """

        self.suffix = suffix
        """
        The suffix of the output file that is used both to extablish for
        post_processing and to establish if the run can be skipped.
        """

        self.command = ('OMP_NUM_THREADS='+ str(self.omp) + ' ' + self.mpi_run +\
                        ' ' + executable).strip()
        print('Initialize a Yambo calculator with command %s' %self.command)

    def pre_processing(self,**kwargs):
        """
        This method checks if the run_dir and the SAVE folder exists, and write
        the yambo input file in the run_dir. The construction of the run_dir and
        the setup of the SAVE folder is managed by the pre_processing features
        of Dataset.
        """
        run_dir=kwargs['run_dir']
        input=kwargs['input']
        name=kwargs['name'] + '.in'

        if not os.path.isdir(run_dir):
            print('Run_dir %s does not exists'%run_dir)
        else:
            if not os.path.isdir(run_dir+'/SAVE'):
                print('SAVE folder does not exists')
            else:
                if not (input is None):
                    input.write(run_dir + '/' + name)
                else :
                    print('input not provided')

    def _seek_output_file(self,run_dir,jobname):
        """
        Check if one the output file (with both extension .hf or .qp) are inside
        the jobname folder. Update the class member self.output
        """
        outfile_hf = run_dir+'/'+jobname+'/o-'+jobname+'.hf'
        outfile_qp = run_dir+'/'+jobname+'/o-'+jobname+'.qp'
        if os.path.isfile(outfile_hf):
            self.output = outfile_hf
        elif os.path.isfile(outfile_qp):
            self.output = outfile_qp
        else:
            self.output = None

    def process_run(self,**kwargs):
        """
        Set the proper dir and run the computation.

        If skip=True run the computations only if the output file (o-$name) is
        not present in the jobname folder. The skip search both for output with
        extension .hf and .qp. If skip is False delete the jobname folder
        (if found) before the run.
        """
        run_dir = kwargs['run_dir']
        jobname = kwargs['name']
        input = jobname+'.in'
        outfolder = run_dir+'/'+jobname # da togliere
        output = jobname+'/o-'+jobname+self.suffix

        string = 'cd %s ; '%run_dir
        string +=  self.command + ' -F %s -J %s -C %s'%(input,jobname,jobname)

        # seek the output file to establish if the run can be skipped
        self._seek_output_file(run_dir,jobname)

        if self.skip:
            if not (self.output is None):
                if self.verbose : print('skip the computation for : '+input)
            else:
                if self.verbose : print('execute : '+string)
                os.system(string)

        if not self.skip:
            if not (self.output is None):
                if self.verbose: print('delete folder : '+outfolder)
                os.system('rm -r %s'%outfolder)
            if self.verbose: print('execute : '+string)
            os.system(string)

        # set self.output for post_processing
        self._seek_output_file(run_dir,jobname)

    def post_processing(self):
        """
        Apply the post processing method of YamboParser to self.output
        """
        results = None
        if not (self.output is None):
            if self.verbose : print('parse file : '+self.output)
            results_dic = dict_parser(self.output)
            results = AttributeDict(**results_dic)
        return results

    def run(self,run_dir='run',input=None,name='test',post_processing=True):
        """
        Prepare the run, perform the computation and apply the post_processing
        function to extract the results. It returns the object built by YamboParser.
        Args:
            run_dir (str) : the folder in which the simulation is performed
            input : the object the contain the instance of the input file
            name (str) : the name associated to the input file (without extension).
            This string is used also as jobname and so it represent the folder in which
            results are written as well as a part of the name of the output file (o-$name.**)
        """
        self.pre_processing(run_dir=run_dir,input=input,name=name)
        self.process_run(run_dir=run_dir,name=name)
        results = None
        if post_processing:
            results = self.post_processing()
        return results
