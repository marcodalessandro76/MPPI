import os
import numpy as np
from mppi.Utilities import Utils, Dos

def test_dict_merge():
    dest = {'a':{'x':1},'b':2}
    Utils.dict_merge(dest,{'a':{'y':3},'c':4})
    assert dest == {'a':{'x':1,'y':3},'b':2,'c':4}

def test_file_parser_skip(tmp_path):
    file = tmp_path/'data.txt'
    file.write_text('! comment\n1 2\n3 4\n')
    columns = Utils.file_parser(str(file),skip='!')
    assert np.allclose(columns,[[1,3],[2,4]])

def test_dos_default_broadening():
    dos = Dos(np.array([0.]),minVal=-1,maxVal=1,step=0.01,eta=0.05)
    x, histo = dos.dos[0]
    assert np.isclose(histo.max(),1/(np.pi*0.05),rtol=1e-3)
    assert np.isclose(np.trapezoid(histo,x),1.0,atol=0.05)

def test_dos_gaussian():
    dos = Dos(np.array([0.]),minVal=-1,maxVal=1,eta=0.05,broad_kind='gaussian')
    x, histo = dos.dos[0]
    assert np.isclose(histo.max(),1/(0.05*np.sqrt(2*np.pi)),rtol=1e-3)

def test_dos_from_pw_set_gap(ref_dir):
    from mppi import Parsers as P
    file = os.path.join(ref_dir,'random_grids','data-file-schema.xml')
    data = P.PwParser(file,verbose=False)
    nb = data.nbands_full
    # with set_gap the empty bands are shifted so that the gap has the requested value
    dos = Dos.from_Pw(file,set_gap=2.0,eta=0.01,step=0.001)
    x, histo = dos.dos[0]
    evals = data.get_evals(set_gap=2.0,verbose=False)
    assert np.isclose(evals[:,nb].min(),2.0)
    assert histo[np.argmin(abs(x-1.0))] < 0.05*histo.max()
