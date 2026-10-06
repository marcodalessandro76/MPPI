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
    mass, pseudo = list(data.atomic_species.values())[0]
    assert isinstance(mass,float)
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

def test_yambodftparser_expand_IBZ_kpoints(ref_dir):
    # WSe2: 12x12x3 Gamma-centered grid, 12 symmetries with inversion, no time reversal. The weights of the IBZ points
    # are compared with the ones of pw.x (normalized to one, same order of the k points)
    folder = os.path.join(ref_dir,'dftParsers_results','WSe2_12x12x3_100bands')
    data = P.YamboDftParser(os.path.join(folder,'ns.db1'),verbose=False)
    kbz = data.expand_IBZ_kpoints(verbose=False)
    assert kbz.shape == (432,3) and data.kpoints_bz.shape == (432,3)
    assert np.all((kbz >= 0) & (kbz < 1))
    grid = kbz*np.array([12,12,3])
    assert np.allclose(grid,np.round(grid),atol=1e-3)
    assert len({tuple(g) for g in np.round(grid).astype(int)}) == 432
    assert data.ibz_index.shape == (432,) and data.sym_index.max() < len(data.syms)
    assert np.isclose(data.weights.sum(),1.) and np.all(data.weights > 0)
    pw = P.PwParser(os.path.join(folder,'data-file-schema.xml'),verbose=False)
    w_pw = np.array(pw.weights).ravel()
    assert np.allclose(data.weights,w_pw/w_pw.sum())
    minus = data.get_minus_k_indexes()
    assert np.all(minus >= 0) and np.all(minus[minus] == np.arange(len(minus)))

def test_yambodftparser_expand_without_inversion(ref_dir):
    # 2 symmetries without inversion and no time reversal: the 222 IBZ points give the full 12x12x3 grid,
    # which is closed under k -> -k
    data = P.YamboDftParser(os.path.join(ref_dir,'rt_results','ns.db1'),verbose=False)
    assert not data.time_reversal
    assert not any(np.allclose(s,-np.eye(3)) for s in data.syms)
    kbz = data.expand_IBZ_kpoints(verbose=False)
    assert kbz.shape == (432,3)
    assert np.isclose(data.weights.sum(),1.)
    minus = data.get_minus_k_indexes()
    assert np.all(minus >= 0) and np.all(minus[minus] == np.arange(len(minus)))

def test_yambodftparser_minus_k_requires_expansion(ref_dir):
    data = P.YamboDftParser(os.path.join(ref_dir,'rt_results','ns.db1'),verbose=False)
    with pytest.raises(AttributeError):
        data.get_minus_k_indexes()

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
    assert data.n_runs == 1 and data.N_ext_fields == 1
    assert data.Efield[0]['name'] == 'DELTA'
    assert data.Polarization[0].shape == (3,len(data.IO_TIME_points))
    # plain numpy arrays (not netCDF masked arrays)
    assert not isinstance(data.Polarization[0],np.ma.MaskedArray)
    assert data.E_ext[0].shape == (3,len(data.IO_TIME_points)) and np.iscomplexobj(data.E_ext[0])
    from mppi.Utilities.Constants import FsToAu
    assert np.allclose(data.get_time(convert_to_fs=False),data.IO_TIME_points)
    assert np.allclose(data.get_time()*FsToAu,data.IO_TIME_points)

def test_nl_output_file(ref_dir):
    folder = os.path.join(ref_dir,'nl_results','LiF-delta_pulse')
    out = P.YamboOutputParser.from_file(os.path.join(folder,'o-lresponse-bands_3-6-delta.NL_pol_F1'),verbose=False)
    pol = out['NL_pol_F1']
    assert list(pol.keys()) == ['time','Pol_x','Pol_y','Pol_z','Dip_x','Dip_y','Dip_z']
    # the polarization of the o- file is the one of the database
    data = P.YamboNLDBParser(os.path.join(folder,'ndb.Nonlinear'),verbose=False)
    assert np.allclose(pol['time'],data.get_time())
    assert np.allclose(pol['Pol_x'],data.Polarization[0][0],rtol=1e-4,atol=1e-20)
