"""
This module contains some tools to perform parallel procedures using the python
multiprocessing package

"""
import multiprocessing as mp, time as tm
import numpy as np
from datetime import timedelta

def loop(func, pars, *args, ntasks = 4, verbose = True, **kwargs):
    """
    Perform a parallel loop over the values of the pars array and compute the
    values of the function func, using ntasks parallel processes

    Args:
        func (function) : a function that returns a value for each element of pars
        pars (:py:class:`array`) : array with the values iterate by the loop
        ntask (:py:class:`int`) : number of parallel tasks
        verbose (:py:class:`bool`) : determine the amount of information provided on terminal
        args, kwargs : arguments and keyword arguments passed to func

     """
    def func_loop(func,pars_subset,task,output,*args,**kwargs):
        """
        Evaluate the function func for all the values inside a single task.
        Add the dictionary with the results of the task to the queue of the multiprocess

        """
        results = []
        for p in pars_subset:
            results.append(func(p,*args,**kwargs))
        output.put({task:np.array(results)})

    pars_split = np.array_split(pars,ntasks)
    if verbose : print('Run a parallel loop with %s tasks...'%ntasks)
    t0 = tm.time()
    output = mp.Queue()
    tasks = [mp.Process(target=func_loop, args=(func,pars_split[task],task,output,*args,), kwargs=kwargs) for task in range(ntasks)]
    for p in tasks:
        p.start()
    results_dict = {}
    for p in tasks:
        results_dict.update(output.get())
    results = np.concatenate([results_dict[i] for i in range(ntasks)])
    if verbose :
        deltaTime = int(tm.time()-t0)
        dT_str = "{:0>8}".format(str(timedelta(seconds=deltaTime)))
        print('Loop executed in',dT_str)
    return results
