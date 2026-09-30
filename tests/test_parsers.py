import os
import numpy as np
import pytest
from mppi import Parsers as P
from mppi.Parsers import ParsersUtils

def test_pwparser(ref_dir):
    data = P.PwParser(os.path.join(ref_dir,'random_grids','data-file-schema.xml'),verbose=False)
    assert data.nkpoints == 100
    assert data.nbands == 8
    assert data.evals.shape == (100,8)
    assert np.isclose(data.weights.sum(),2.0)
    gap = data.get_gap(verbose=False)
    assert np.isclose(gap['gap'],0.465234,atol=1e-5)
    assert np.isclose(gap['gap'],gap['direct_gap'])

def test_pwparser_missing_file(tmp_path):
    data = P.PwParser(str(tmp_path/'missing.xml'),verbose=False)
    assert data.data is None

def test_pwparser_evals_scissor(ref_dir):
    data = P.PwParser(os.path.join(ref_dir,'random_grids','data-file-schema.xml'),verbose=False)
    evals = data.get_evals(set_gap=1.0,verbose=False)
    assert np.isclose(evals[:,data.nbands_full].min()-evals[:,data.nbands_full-1].max(),1.0)

def test_yambodftparser(ref_dir):
    data = P.YamboDftParser(os.path.join(ref_dir,'rt_results','ns.db1'),verbose=False)
    assert data.nkpoints == 222
    assert data.nbands == 100
    assert data.evals.shape == (222,100)

def test_dft_parsers_agree(ref_dir):
    folder = os.path.join(ref_dir,'dftParsers_results','WSe2_12x12x3_100bands')
    pw = P.PwParser(os.path.join(folder,'data-file-schema.xml'),verbose=False)
    yambo = P.YamboDftParser(os.path.join(folder,'ns.db1'),verbose=False)
    assert pw.nbands == yambo.nbands
    assert pw.nkpoints == yambo.nkpoints
    assert np.allclose(pw.get_lattice(),yambo.get_lattice(),atol=1e-4)
    assert np.isclose(pw.get_gap(verbose=False)['gap'],yambo.get_gap(verbose=False)['gap'],atol=1e-3)

def test_yambooutputparser_rt(ref_dir):
    data = P.YamboOutputParser.from_path(os.path.join(ref_dir,'rt_results'),verbose=False)
    assert set(data.keys()) == {'carriers','current','external_field','orbt_magnetization',
                                'polarization','spin_magnetization'}
    assert list(data['orbt_magnetization'].keys()) == ['time','Ml_x','Ml_y','Ml_z','Mi_x','Mi_y','Mi_z']
    assert list(data['polarization'].keys()) == ['time','Pol_x','Pol_y','Pol_z']

def test_yambooutputparser_hf_qp(ref_dir):
    hf = P.YamboOutputParser.from_file(os.path.join(ref_dir,'hf_results','o-hf_test1.hf'),verbose=False)
    assert list(hf['hf'].keys()) == ['kpoint','band','E0','Ehf','Vxc','Vnlxc']
    qp = P.YamboOutputParser.from_file(os.path.join(ref_dir,'qp_results','o-qp_test1.qp'),
                                       verbose=False,extendOut=False)
    assert list(qp['qp'].keys()) == ['kpoint','band','E0','EmE0','sce0']
    k, b = hf['hf']['kpoint'][0], hf['hf']['band'][0]
    assert hf.get_energy(k,b) == hf['hf']['Ehf'][0]

def test_rtcarriers(ref_dir):
    data = P.YamboRTCarriersParser(os.path.join(ref_dir,'rt_results','ndb.RT_carriers'),verbose=False)
    assert data.delta_f.shape == (16,1332)
    assert data.E_bare.shape == (1332,)

def test_get_variable_from_db(ref_dir):
    kpts = ParsersUtils.get_variable_from_db(os.path.join(ref_dir,'rt_results','ns.db1'),'K-POINTS')
    assert kpts.shape == (3,222)

NL_DELTA = os.path.join('nl_results','LiF-delta_pulse','ndb.Nonlinear')

def test_nldbparser(ref_dir):
    file = os.path.join(ref_dir,NL_DELTA)
    if not os.path.isfile(file): pytest.skip('nl_results reference data not available')
    data = P.YamboNLDBParser(file,verbose=False)
    assert data.n_runs == 1
    assert data.Efield[0]['name'] == 'DELTA'
    assert data.Polarization[0].shape == (3,len(data.IO_TIME_points))
