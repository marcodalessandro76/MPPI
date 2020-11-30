"""
This module defines the tools to perform and manage several calculations.
The usage of this module aims to simplify the approach to an ensemble calculations
using both QuantumESPRESSO and Yambo, and to deal with parallel executions of multiple
instances of the code.
"""

from mppi.Calculators.Runner import Runner

def name_from_id(id):
    """
    Convert the id into a run name. If id is a string, set name = id, if it is a
    dictionary build the name string of the run from the id dictionary.

    Args:
        id : id associated to the run
    Returns:
       name (:py:class:`str`): name of the run associated to the dictionary ``id``
    """
    if type(id) is str :
        name = id
    elif type(id) is dict :
        keys=sorted(id.keys())
        name=''
        for k in keys:
            name += k+'_'+str(id[k])+'-'
        name = name.rstrip('-')
    else :
        print('id type not recognized')
        name = None
    return name

class Dataset(Runner):
    """
    Class to perform a set of calculations and to manage the associated results.

    Parameters:
        label (:py:class:`str`): the label of the dataset, it can be useful for instance if more
            than one istance of the class is present
        run_dir (:py:class:`str`): path of the directory where the runs will be performed. This argument
            can be overwritten by including a run_dir keyword in the :meth:`append_run` method of the class.
            In this way the various elements of the dataset can be run in different folders
        num_tasks  (:py:class:`int`) : maximum number of computations performed in parallel by
            the :meth:`run` method of the class
        verbose (:py:class:`bool`) : set the amount of information provided on terminal
        **kwargs : all the parameters passed to the dataset and stored in its _global_options.
            Can be useful, for instance, in performing a post-processing of the results

    The class members are:

    Attributes:
        ids (:py:class:`list`) : list of run ids, to be used in order to identify and fetch the results
        runs (:py:class:`list`) : list of the runs which have to be treated by the dataset. The runs
            contain all the input parameters to be passed to the various runners.
        calculators (:py:class:`list`) : calculators which will be used by the run method
        results (:py:class:`dict`) : set of the results of each of the runs. The set is not ordered as the
            runs may be executed asynchronously.
        post_processing_function (function) : specify a postprocessing on the results provided by the runs
            of the dataset

    Example:
        >>> code = QeCalculator()
        >>> study=Dataset(label = .., run_dir = ..., **kwargs)
        >>> study.append_run(id={'ecut': 30, 'kpoints' : 4},input=...,runner=code,variable1=1)
        >>> study.append_run(id={'ecut': 40, 'kpoints' : 4},input=...,runner=code,variable2='periodic')
        >>> study.run()

    """

    def __init__(self, label = 'Dataset', run_dir = 'runs', num_tasks = 2, verbose = True, **kwargs):
        """
        Set the dataset ready for appending new runs
        """
        from copy import deepcopy
        kwargs = deepcopy(kwargs)
        Runner.__init__(self,label=label,run_dir=run_dir,num_tasks=num_tasks,verbose=verbose,**kwargs)
        print('Initialize a Dataset with %s parallel tasks'%self._global_options['num_tasks'])

        self.ids = []
        self.runs = []
        self.calculators = []
        self.results = {}
        self.post_processing_function = None

    def append_run(self, id, runner, **kwargs):
        """
        Add a run into the dataset.

        Append a run to the list of runs to be performed and associate to each appended
        item the corresponding runner instance.
        If the name of the input file is not provided, the method attribute it from the
        id of the run using the function name_from_id. If, for instance, a jobname
        has to be provided it can be passed as kwargs.

        Args:
            id : the id of the run, useful to identify the run in the dataset. It can be
                a dictionary or a string, as it may contain different keyword. For example a
                run can be classified as
                ``id = {'energy_cutoff': 60, 'kpoints': 6}``
            runner (:class:`Runner`) : the instance of  :class:`runner` class to which the
                keyword arguments will be passed at the input
            kwargs : these arguments contain the instance of the input and any other variable needed for
                 appended run. All these quantities are stored as an element of the runs list
                 and are passed to the calculator, together with the global options of the Dataset,
                 when the run method is called

        Raises:
          ValueError: if the provided id is identical to another previously
             appended run.

        """
        from copy import deepcopy
        if id in self.ids:
            raise ValueError('The run id', id,
                             ' is already provided, modify the run id.')
        irun = len(self.ids) # get the cardinal number of this run and append its id
        self.ids.append(id)
        # check if the calculator has been already used, otherwise add it to self.calculators
        # icalc identifies the position of the calculator in the self.calculators list
        calc_found = False
        for ind,calc in enumerate(self.calculators):
            if calc['calc'] == runner:
                calc['iruns'].append(irun)
                calc_found = True
                break
        if not calc_found:
            self.calculators.append({'calc': runner, 'iruns': [irun]})
        #add the kwargs and the global_options to the self.runs
        inp_to_append = deepcopy(self._global_options)
        inp_to_append.update(deepcopy(kwargs))
        #add the 'name' keyword if not provided by the user
        if not 'name' in inp_to_append:
            inp_to_append['name'] = name_from_id(id)
        self.runs.append(inp_to_append)

    def process_run(self):
        """
        Run the dataset by performing explicit run of each of the item of the
           runs list. If the list selection is provided in the call of the :py:meth:'run'
           the calculation is restricted to the elements of the list

        """
        selection = self.run_options.get('selection',None)
        self.run_the_calculations(selection)
        return {}

    def run_the_calculations(self, selection):
        """
        Method that manage the execution of the runs of the Dataset. The elements of the Dataset
        in the selection list are computed in parallel according to the limitation provided by the
        num_tasks attribute. The method uses the :py:class:`multiprocessing.Process` to manage the
        parallel runs.

        Args:
            selection (:py:class:`list`) : if not None only the runs in the list are computed.
                This parameter is used only when the method is called by the :meth:`fetch_results`
                method.

        """
        import multiprocessing, time
        delay = 1 # in seconds
        verbose = self.global_options().get('verbose')

        def calculator_run(runs,calculators,iruns,queue):
            for calc in calculators: #identify the calculator associated to the present run
                if irun in calc['iruns']:
                    break
            #run and append the dictionary with the result to the queue
            result = calc['calc'].run(**runs[irun])
            queue.put({irun : result})

        if selection is None:
            selection = [ind for ind in range(len(self.ids))]
        task_groups = self.build_taskgroups(selection)
        if verbose: print('Run the selection %s with the parallel task_groups %s \n'%(selection,task_groups))
        for task in task_groups:
            queue = multiprocessing.Queue()
            task_job = []
            task_alive = True
            if verbose: print('Run the task %s '%task)
            for irun in task:
                p = multiprocessing.Process(target=calculator_run, args=(self.runs,self.calculators,irun,queue,))
                task_job.append(p)
                p.start()
            while task_alive: #wait the end of the task
                task_alive = any([job.is_alive() for job in task_job])
                time.sleep(delay)
            while queue.qsize() != 0: #add the result to self.results
                self.results.update(queue.get())
            if verbose: print('Task %s ended \n'%task)

    def build_taskgroups(self,selection):
        """
        Identify the elements,among the runs provided as input, that can be executed in parallel.
        The number of parallel computations is specified by the self.num_tasks attribute of the class

        Args:
            selection (:py:class:`list`) : list with a selection of the runs of dataset.

        Return:
            (:py:class:`list`) : a list of list with the groups of parallel computation. The order
            of the runs respects the ones of the selection list

        """
        num_tasks = self.global_options().get('num_tasks')
        task_groups = [selection[i:i + num_tasks] for i in range(0,len(selection),num_tasks)]
        return task_groups

    def set_postprocessing_function(self, func):
        """
        Set the callback of run.
        Calls the function ``func`` after having performed the appended runs.

        Args:
           func (func): function that process the `inputs` `results` and
               returns the value of the `run` method of the dataset.
               The function is called as ``func(self)``.
        """
        self.post_processing_function = func

    def post_processing(self, **kwargs):
        """
        Calls the Dataset function with the results of the runs as arguments
        """
        if self.post_processing_function is not None:
            return self.post_processing_function(self)
        else:
            return self.results

    def fetch_results(self, id=None, attribute=None, run_if_not_present=True):
        """
        Retrieve the results that match some conditions that is specified through
        an `id` in the form of a string or a dictionary. Selects out of the results
        of the objects which have in their ``name`` keyword at least the `id` provided as input.

        Args:
           id : string or dictionary of the retrieved id.
           attribute (str): if present, provide the attribute of each of the
               results instead of the result object
           run_if_not_present (bool): If the run has not yet been performed
               in the dataset then perform it.

        Return:
            A list of the runs that match the condition in the order provided by
            :py:meth:`append_run` method.

        Example:
           >>> study=Dataset()
           >>> study.append_run(id={'ecut': 40, 'k' : 4}, input = ..., runner = )
           >>> study.append_run(id={'ecut': 40, 'k' : 6}, input = ..., runner = )
           >>> study.append_run(id={'ecut': 50, 'k' : 6}, input = ..., runner = )
           >>> #append other runs if needed
           >>> #set a post processing function that perform a parsing of the rsesults
           >>> #and contains 'energy' as an attribute of the results object
           >>> #run the calculations (optional if run_if_not_present=True)
           >>> study.run()
           >>> # returns a list of the energies of first and the second result
           >>> # in this example
           >>> data=study.fetch_results(id={'ecut': 40},attribute='energy')

        """

        #names = [name_from_id(id) for id in self.ids]
        names = [val['name'] for val in self.runs]
        id_name = name_from_id(id)
        fetch_indices = []
        selection_to_run = []
        for irun,name in enumerate(names):
            if id_name in name :
                fetch_indices.append(irun)
                if run_if_not_present and irun not in self.results:
                    selection_to_run.append(irun)
        if len(selection_to_run) > 0:
            self.run_the_calculations(selection=selection_to_run)

        data = []
        if self.post_processing_function is not None:
            for irun in fetch_indices:
                r = self.post_processing()[irun]
                data.append(r if attribute is None else getattr(r, attribute))
        else:
            print('Provide a post processing function able to parse the results')
        return data

    def seek_convergence(self, rtol=1.e-5, atol=1.e-8,
                        selection=None, **kwargs):
        """
        Search for the first result of the dataset which matches the provided
        tolerance parameter. The results are in dataset order
        (provided by the :py:meth:`append_run` method) if `selection` is not
        specified.
        Employs the numpy :py:meth:`allclose` method for comparison.

        Args:
          rtol (float): relative tolerance parameter
          atol (float): absolute tolerance parameter
          selection (list): list of the id of the runs in which to perform the
               convergence search. Each id should be unique in the dataset.
          **kwargs: arguments to be passed to the :py:meth:`fetch_results`
               method.

        Returns:
          id,result (tuple): the id of the last run which matches the
                convergence, together with the result, if convergence is
                reached.

        Raises:
           LookupError: if the parameter for convergence were not found.
               The dataset has to be enriched or the convergence parameters
               loosened.

        """
        from numpy import allclose
        to_get = self.ids if selection is None else selection

        id_ref = to_get[0]
        print('Fetching results for id "', id_ref, '"')
        ref = self.fetch_results(id=id_ref, **kwargs)
        ref = ref[0]
        for id in to_get[1:]:
            print('Fetching results for id "', id, '"')
            val = self.fetch_results(id=id, **kwargs)
            val = val[0]
            if allclose(ref, val, rtol=rtol, atol=atol):
                res = self.fetch_results(id=id_ref, **kwargs)
                label = self.get_global_option('label')
                print('Convergence reached in Dataset "' +
                      label+'" for id "', id_ref, '"')
                return (id_ref, res[0])
            ref = val
            id_ref = id
        raise LookupError('Convergence not reached, enlarge the dataset'
                          ' or change tolerance values')
