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

def test_fetch_results():
    study = build_dataset()
    assert study.fetch_results(id={'ecut':4},run_if_not_present=False) == [42,44]
    assert study.fetch_results(id={'ecut':40,'k':4},run_if_not_present=False) == [404]
    assert study.fetch_results(id={'k':2},run_if_not_present=False) == [42,402,4002]
