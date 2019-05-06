"""
This module defines a class to perform and manage several calculations.
The module is (deeply) inspired from the Dataset class of BigDFT.

A major difference with respect to the Dataset class in BigDFT is given by the
usage of the pre_processing function. In the present case this function performs
different operation dependending of the type of simulation managed by the dataset
(scf, nscf, yambo,...). The pre_processing is called independently from the run
method, this is important in particular for Yambo. Indeed the pre_processing manage
the building of the run_dir with the SAVE folder, and the YamboIn needs this folder
to init its object so the pre_processing in this case _has to_ be called before
append_run.
"""

import os
from copy import deepcopy

# The __init__ file used to load this module allows us to call this methods also
# outside from the class. It can be used, for instance, for setting the prefix in
# the qe input file from the id dictionary that identify the computation
def name_from_id(id):
    """
    Hash the id into a run name
    Construct the name of the run from the id dictionary
    Args:
        id (dict): id associated to the run
    Returns:
       str: name of the run associated to the dictionary ``id``
    """
    keys=sorted(id.keys())
    name=''
    for k in keys:
        name += k+':'+str(id[k])+','
    return name.rstrip(',')

class Dataset():
    """
    Class to manage a set of calculations.

    This class contains the various instances of a set of calculations specified by the their calculator object.
    The different calculations are labelled by the values of some relevant parameters given as a dictionary.

    Args:
      label (str): The label of the dataset. It will be needed to identify the instance for example
          in plot titles or in the running directory.
      run_dir (str): path of the directory where the runs will be performed.

    """
    def __init__(self,label='dataset',run_dir='dataset_runs',pre_processing=None):
        """
        Set the dataset ready for appending new runs
        """
        self.label=label

        self.run_dir=run_dir

        self.runs=[]
        """List of the runs which have to be treated by the dataset these runs contain the input
        parameter to be passed to the various runners.
        """

        self.calculators=[]
        """
        Calculators which will be used by the run method, useful to gather the inputs in the case of a
        multiple run.
        """

        self.results={}
        """
        Set of the results of each of the runs. The set is not ordered as the runs may be executed
        asynchronously.
        """

        self.ids=[]
        """
        List of run ids, to be used in order to classify and fetch the results.
        """

        self.names=[]
        """
        List of run names, needed for distinguishing the logfiles and input files.
        """

        self.pre_processing = pre_processing
        """
        Specifies the type of pre_processing_function to be call before running the dataset.
        """

    def scf_pre_processing(self,**kwargs):
        """
        Define the pre_processing function for a scf dataset.
        The method build the run_dir folder if it does not exists.
        """
        if not os.path.isdir(self.run_dir):
            os.mkdir(self.run_dir)
            print('Create folder %s'%self.run_dir)

    def nscf_pre_processing(self,**kwargs):
        """
        Define the pre_processing function for a nscf dataset.
        The method build the run_dir folder if it does not exists. Then, for each id of
        the dataset copy the source_dir with name name_from_id(id).save.
        """
        source = kwargs['source_dir']
        if not os.path.isdir(self.run_dir):
            os.mkdir(self.run_dir)
            print('Create folder %s'%self.run_dir)

        for id in self.ids:
            dest = self.run_dir + '/' + name_from_id(id) + '.save'
            # check if the source folder exists
            if not os.path.isdir(source): print('source folder : ',source,'not found')
            else :
                # copy the source folder only if the dest is not present
                if not os.path.isdir(dest):
                    string = 'cp -r %s %s'%(source,dest)
                    print('execute : ',string)
                    os.system(string)
                else : print('nscf save folder already exsists. Source .save NOT COPIED')

    def yambo_pre_processing(self,**kwargs):
        """
        Define the pre_processing function for a Yambo dataset.
        The method build the run_dir folder if it does not exists.Then check if the SAVE
        folder exists, if not runs p2y in the source_dir and copy the SAVE folder in the run_dir.
        Lastly, executes yambo (without arguments) in the run_dir to build the r_setup
        """
        source = kwargs['source_dir']
        if not os.path.isdir(self.run_dir):
            os.mkdir(self.run_dir)
            print('Create folder %s'%self.run_dir)
        if not os.path.isdir(source):
            print('source folder : ',source,' not found')
        elif os.path.isdir(self.run_dir +'/SAVE'):
            print('SAVE folder already present in %s'%self.run_dir)
        else:
            # run p2 y
            string = 'cd %s;p2y'%source
            print('execute : ',string)
            os.system(string)
            # copy the SAVE folder
            savef = source + '/SAVE'
            string = 'cp -r %s %s'%(savef,self.run_dir)
            print('execute : ',string)
            os.system(string)
            # build the r_setup
            string = 'cd %s;OMP_NUM_THREADS=1 yambo'%self.run_dir
            print('execute : ',string)
            os.system(string)

    def pre_processing_function(self,**kwargs):
        """
        Choose a pre_processing function among the ones provided in the dataset.
        """
        if self.pre_processing == None: print('Specify a pre_processing for the dataset')
        if self.pre_processing == 'scf': self.scf_pre_processing(**kwargs)
        if self.pre_processing == 'nscf': self.nscf_pre_processing(**kwargs)
        if self.pre_processing == 'yambo': self.yambo_pre_processing(**kwargs)

    def append_run(self,id,calculator,input):
        """
        Add a run into the dataset.

        Append to the list of runs to be performed the corresponding calculator and the arguments which are
        associated to it.

        Args:
          id (dict): the id of the run, useful to identify the run in the dataset. It has to be a dictionary and,
          typically its elements are used to build the associated input.
          calculator : the calculator class to which the remaining keyword arguments will be passed at the input.

        Raises:
          ValueError: if the provided id is identical to another previously appended run.

        """
        name=name_from_id(id)
        if name in self.names:
            raise ValueError('The run id',name,' is already provided, modify the run id.')
        self.names.append(name)

        #get the number of this run
        irun=len(self.runs)
        #append it to the runs list
        self.runs.append(deepcopy(input))
        #append id and name
        self.ids.append(deepcopy(id))
        #search if the calculator already exists
        found = False
        for calc in self.calculators:
            if calc['calc'] == calculator:
                calc['runs'].append(irun)
                found=True
                break
        if not found:
            self.calculators.append({'calc': calculator, 'runs':[irun]})

    def process_run(self,post_processing):
        """
        Run the dataset, by performing explicit run of each of the item of the runs_list.
        If post_processing is True apply the post_processing function so that self.results
        contains the results of all the element of the dataset.
        """
        for c in self.calculators:
            calc=c['calc']
            for r in c['runs']:
                input=self.runs[r]
                name=self.names[r]
                self.results[r] = calc.run(input=input,name=name,run_dir=self.run_dir,post_processing=post_processing)

    def run(self,post_processing=True):
        """
        Execute process_run
        """
        #self.pre_processing_function()
        self.process_run(post_processing)

    def fetch_results(self,id=None,attribute=None):
        """
        Retrieve some attribute from some of the results.

        Selects out of the results the objects which have in their ``id``
        at least the dictionary specified as input. May return an attribute
        of each result if needed.

        Args:
           id (dict): dictionary of the retrieved id. Return a list of the runs that
               have the ``id`` argument inside the provided ``id`` in the order provided by :py:meth:`append_run`.
           attribute (str): if present, provide the attribute of each of the results instead of the result object

        Example:
           >>> study=Dataset()
           >>> study.append_run(id={'cr': 3},input=input1)
           >>> study.append_run(id={'cr': 4},input=input2)
           >>> #append other runs if needed
           >>> study.run()  #run the calculations
           >>> # returns a list of the energies of first and the third result in this example
           >>> data=study.fetch_results(id={'cr': 3},attribute='E_tot')
        """
        name='' if id is None else name_from_id(id)
        data=[]
        for irun,n in enumerate(self.names):
            if name not in n: continue
            r=self.results[irun]
            data.append(r if attribute is None else getattr(r,attribute))
        return data
