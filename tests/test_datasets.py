from mppi.Calculators.Runner import Runner
from mppi.Datasets import Dataset
from mppi.Datasets.Dataset import id_matches, name_from_id

def test_name_from_id():
    assert name_from_id('run') == 'run'
    assert name_from_id({'k':4,'ecut':40}) == 'ecut_40-k_4'
    assert name_from_id((1,'a')) == '1-a'

def test_id_matches():
    assert id_matches({'ecut':40,'k':4},'ecut_40-k_4',{'ecut':40})
    assert not id_matches({'ecut':40,'k':4},'ecut_40-k_4',{'ecut':4})
    assert not id_matches({'ecut':40},'ecut_40',{'ecut':40,'k':4})
    assert id_matches('LiF-delta','LiF-delta','delta')
    assert not id_matches('LiF-delta_pulse','LiF-delta_pulse','delta')

def build_dataset():
    study = Dataset(run_dir='runs',verbose=False)
    code = Runner()
    for ecut in [4,40,400]:
        for k in [2,4]:
            study.append_run(id={'ecut':ecut,'k':k},runner=code)
    # fill the results without running the calculations
    for irun,id in enumerate(study.ids):
        study.results[irun] = 10*id['ecut']+id['k']
    study.set_postprocessing_function(lambda dataset: dataset.results)
    return study

class SquareRunner(Runner):
    """Runner that returns the square of its x option, or fails if x is negative"""
    def post_processing(self):
        x = self.run_options['x']
        if x < 0: raise ValueError('negative x')
        return x**2

def test_dataset_run_multiprocessing():
    # on Windows and macOS the processes are started with spawn, so this also tests that the
    # objects passed to the processes can be pickled
    study = Dataset(run_dir='runs',num_tasks=2,verbose=False)
    code = SquareRunner()
    for x in [1,2,3,-1]:
        study.append_run(id={'x':x},runner=code,x=x)
    assert study.run() == {0:1,1:4,2:9,3:None}
    study.set_postprocessing_function(lambda dataset: dataset.results)
    assert study.fetch_results(id={'x':3}) == [9]

def test_pw_postprocessing_failed_run(ref_dir, tmp_path):
    import os
    from mppi.Datasets import PostProcessing as PP
    study = Dataset(verbose=False)
    study.results = {0:os.path.join(ref_dir,'random_grids','data-file-schema.xml'),
                     1:str(tmp_path/'missing.xml')}
    energy = PP.pw_get_energy(study)
    assert energy[1] is None and isinstance(energy[0],float)
    assert PP.pw_get_gap(study)[1] is None

def test_parallel_loop():
    import math
    import numpy as np
    from mppi.Utilities.Parallel import loop
    pars = np.arange(10.)
    assert np.allclose(loop(math.sqrt,pars,ntasks=3,verbose=False),np.sqrt(pars))

def test_fetch_results():
    study = build_dataset()
    assert study.fetch_results(id={'ecut':4},run_if_not_present=False) == [42,44]
    assert study.fetch_results(id={'ecut':40,'k':4},run_if_not_present=False) == [404]
    assert study.fetch_results(id={'k':2},run_if_not_present=False) == [42,402,4002]
