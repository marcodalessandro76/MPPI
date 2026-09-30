import os
import pytest
from mppi import InputFiles as I

def test_pwinput_parse(io_dir):
    inp = I.PwInput(os.path.join(io_dir,'si_scf.in'))
    assert inp.get_prefix() == 'si_scf'
    assert inp['system']['ecutwfc'] == 40
    assert inp['atomic_species']['Si'][1] == 'Si.pbe-mt_fhi.UPF'
    assert inp['atomic_positions']['type'] == 'crystal'
    assert inp['atomic_positions']['values'][1] == ['Si',[-0.125,-0.125,-0.125]]
    assert inp['kpoints'] == {'type':'automatic','values':([4.,4.,4.],[0.,0.,0.])}

def test_pwinput_roundtrip(io_dir, tmp_path):
    inp = I.PwInput(os.path.join(io_dir,'si_scf.in'))
    out = str(tmp_path/'si.in')
    inp.write(out)
    inp2 = I.PwInput(out)
    for key in inp.namelist + inp.cards:
        assert inp2[key] == inp[key]

PW_TEMPLATE = """&control
    calculation = 'relax'
/
&system
    ibrav = 2
    celldm(1) = 10.3
    nat = 2
    ntyp = 1
    ecutwfc = 30
/
&electrons
/
ATOMIC_SPECIES
  Si   28.086    Si.upf
ATOMIC_POSITIONS crystal
 Si 0.0 0.0 0.0 0 0 0
 Si 0.25 0.25 0.25
{kpoints}
"""

def test_pwinput_if_pos(tmp_path):
    file = tmp_path/'relax.in'
    file.write_text(PW_TEMPLATE.format(kpoints='K_POINTS gamma'))
    inp = I.PwInput(str(file))
    assert inp['atomic_positions']['values'][0] == ['Si',[0.,0.,0.],[0,0,0]]
    assert inp['atomic_positions']['values'][1] == ['Si',[0.25,0.25,0.25]]
    assert inp['kpoints'] == {'type':'gamma','values':[]}
    # the if_pos flags and the gamma card are preserved in the written input
    string = inp.convert_string()
    assert string.splitlines()[-3].split()[-3:] == ['0','0','0']
    assert string.splitlines()[-1].strip() == 'K_POINTS { gamma }'

def test_pwinput_kpoints_list(tmp_path):
    file = tmp_path/'bands.in'
    file.write_text(PW_TEMPLATE.format(kpoints='K_POINTS tpiba_b\n2\n0 0 0 10\n0 0 1 0'))
    inp = I.PwInput(str(file))
    assert inp['kpoints'] == {'type':'tpiba_b','values':[[0.,0.,0.,10.],[0.,0.,1.,0.]]}

def test_pwinput_wrong_kpoints(tmp_path):
    file = tmp_path/'bands.in'
    file.write_text(PW_TEMPLATE.format(kpoints='K_POINTS tpiba_b\n3\n0 0 0 10'))
    with pytest.raises(ValueError):
        I.PwInput(str(file))

def test_phinput_parse(tmp_path):
    inp = I.PhInput(prefix='si',outdir='out')
    inp.set_kpoints([[0.,0.,0.,1]])
    file = str(tmp_path/'ph.in')
    inp.write(file)
    inp2 = I.PhInput(file)
    assert inp2.get_prefix() == 'si'
    assert inp2.get_outdir() == 'out'
    assert inp2['inputph']['tr2_ph'] == 1e-12

YAMBO_INPUT = """HF_and_locXC
EXXRLvcs= 20 Ry
%QPkrange
 1| 1| 4| 5|
%
"""

def test_yamboinput_from_file(tmp_path):
    (tmp_path/'yambo.in').write_text(YAMBO_INPUT)
    inp = I.YamboInput(folder=str(tmp_path))
    assert inp['arguments'] == ['HF_and_locXC']
    assert inp['variables']['EXXRLvcs'] == [20.,'Ry']
    assert inp['variables']['QPkrange'] == [[1,1,4,5],'']
    inp.set_bandRange(3,6)
    assert inp['variables']['QPkrange'] == [[1,1,3,6],'']

def test_yamboinput_missing_file(tmp_path):
    with pytest.raises(IOError):
        I.YamboInput(folder=str(tmp_path),filename='missing.in')

@pytest.mark.requires_yambo
def test_yamboinput_generation(ref_dir, tmp_path):
    # to be validated on a machine with yambo: the SAVE folder contains only the ns.db1 database
    import shutil
    os.makedirs(tmp_path/'SAVE')
    shutil.copy(os.path.join(ref_dir,'dftParsers_results','WSe2_12x12x3_100bands','ns.db1'),tmp_path/'SAVE')
    inp = I.YamboInput('yambo -x -V rl',folder=str(tmp_path))
    assert 'HF_and_locXC' in inp['arguments']
