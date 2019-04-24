"""
This module defines a class to perform a calculations using Yambo.
The module is (deeply) inspired from the SystemCalculator class of BigDFT
"""

import os
import YamboParser as yp

class YamboCalculator():
    """
    Manage a single calculation. Setup the number of omp and mpi.
    Prepare the folder in which the computation is run. Peform the computation
    and apply a post processing function to extract the results

    The Class assumes that the run_dir where yambo is executed exists and that
    it contains the SAVE folder with the p2y postprocessing of a QE computation.
    The setup of this folder is managed at the level of Dataset.
    """
    def __init__(self,omp=1,mpi_run='mpirun -np 4',executable='',skip=False, verbose=True):
        self.omp=omp
        self.mpi_run=mpi_run
        self.executable=executable
        self.skip=skip
        self.verbose=verbose
        #the o-*.hf or o-*.qp file (with its complete path)
        self.output = None
        self.command = ('OMP_NUM_THREADS='+ str(self.omp) + ' ' + self.mpi_run + ' ' + executable).strip()
        print('Initialize a Yambo calculator with command %s' %self.command)

    def pre_processing(self,**kwargs):
        """
        This method checks if the run_dir and the SAVE folder exists, and write the yambo input file
        in the run_dir. The setup of the SAVE folder is managed by the pre_processing features of Dataset.
        """
        run_dir=kwargs['run_dir']
        if not os.path.isdir(run_dir):
            print('Run_dir %s does not exists'%run_dir)
        else:
            if not os.path.isdir(run_dir+'/SAVE'):
                print('SAVE folder does not exists')
            else:
                input=kwargs['input']
                name=kwargs['name'] + '.in'
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
        if os.path.isfile(outfile_qp):
            self.output = outfile_qp

    def process_run(self,**kwargs):
        """
        Set the proper dir and run the computation.

        If skip=True run the computations only if the output file (o-$name) is not
        present in the jobname folder. The skip search both for output with extension
        .hf and .qp. If skip is False delete the jobname folder (if found) before the run.
        """
        string = 'cd %s ; '%kwargs['run_dir']
        run_dir = kwargs['run_dir']
        jobname = kwargs['name']
        input = jobname+'.in'
        outfolder = run_dir+'/'+jobname
        string +=  self.command + ' -F %s -J %s -C %s'%(input,jobname,jobname)

        self._seek_output_file(run_dir,jobname)

        if self.skip:
            if not (self.output is None):
                if self.verbose : print('skip the computation for : '+self.output)
            else:
                if self.verbose : print('execute : '+string)
                os.system(string)

        if not self.skip:
            if not (self.output is None):
                if self.verbose: print('delete folder : '+outfolder)
                os.system('rm -r %s'%outfolder)
            if self.verbose: print('execute : '+string)
            os.system(string)

        self._seek_output_file(run_dir,jobname)
        if self.verbose:print('output in : ',self.output)

    def post_processing(self):
        """
        """
        results = None
        if not (self.output is None):
            results = yp.YamboOut(self.output)
        return results

    def run(self,run_dir='run',input=None,name='test',post_processing=True):
        """
        Prepare the run, perform the computation and apply the post_processing
        function to extract the results
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
