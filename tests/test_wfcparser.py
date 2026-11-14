"""
Tests of the WfcParser module on synthetic data (the tests on real wavefunctions are in the tutorial notebook, since
the outputs of the tutorials are not in git).
"""
import numpy as np
import pytest
from scipy.io import FortranFile
from mppi.Parsers.WfcParser import read_qe_wfc, WfcParser, _pseudo_type

A = 6.                                  # simple cubic lattice (bohr)
B = 2 * np.pi / A * np.eye(3)           # reciprocal lattice vectors as rows

def cubic_states(k=(0.1, 0.2, 0.3), nb=3, gmax=2, seed=0):
    # orthonormal random states on the G vectors with |m_i| <= gmax
    r = np.arange(-gmax, gmax + 1)
    mill = np.array(np.meshgrid(r, r, r, indexing='ij')).reshape(3, -1).T
    rng = np.random.default_rng(seed)
    c = rng.normal(size=(len(mill), nb)) + 1j * rng.normal(size=(len(mill), nb))
    c = np.linalg.qr(c)[0].T
    return {'ik': 1, 'k': np.asarray(k) @ B, 'b': B, 'mill': mill, 'bands': np.arange(1, nb + 1), 'c': c}

def write_qe_wfc(filename, k, mill, c, gamma_only=0, npol=1):
    with FortranFile(filename, 'w') as f:
        head = np.zeros(1, dtype='i4,3f8,i4,i4,f8')
        head[0] = (1, k, 1, gamma_only, 1.)
        f.write_record(head)
        f.write_record(np.array([len(mill), len(mill), npol, len(c)], dtype='i4'))
        f.write_record(B.astype('f8'))
        f.write_record(mill.astype('i4'))
        for cn in c:
            f.write_record(cn.astype(complex))

def test_read_qe_wfc(tmp_path):
    s = cubic_states()
    filename = str(tmp_path / 'wfc1.dat')
    write_qe_wfc(filename, s['k'], s['mill'], s['c'])
    w = read_qe_wfc(filename)
    assert np.allclose(w['k'], s['k']) and np.allclose(w['b'], B) and np.array_equal(w['mill'], s['mill'])
    assert np.allclose(w['c'], s['c']) and list(w['bands']) == [1, 2, 3]
    w = read_qe_wfc(filename, bands=[3, 1])
    assert np.allclose(w['c'], s['c'][[2, 0]]) and list(w['bands']) == [3, 1]
    with pytest.raises(ValueError):
        read_qe_wfc(filename, bands=[4])
    write_qe_wfc(filename, s['k'], s['mill'], s['c'], gamma_only=1)
    with pytest.raises(NotImplementedError):
        read_qe_wfc(filename)
    write_qe_wfc(filename, s['k'], s['mill'], s['c'], npol=2)
    with pytest.raises(NotImplementedError):
        read_qe_wfc(filename)

def test_overlap():
    s = cubic_states()
    assert np.allclose(WfcParser.overlap(s, s), np.eye(3))
    # the G vectors are matched through the Miller indices, whatever their order
    perm = np.random.default_rng(1).permutation(len(s['mill']))
    s2 = dict(s, mill=s['mill'][perm], c=s['c'][:, perm])
    assert np.allclose(WfcParser.overlap(s, s2), np.eye(3))
    # the same states written at k + b1
    s3 = dict(s, k=s['k'] + B[0], mill=s['mill'] - [1, 0, 0])
    assert np.allclose(WfcParser.overlap(s, s3), np.eye(3))
    with pytest.raises(ValueError):
        WfcParser.overlap(s, dict(s, k=s['k'] + 0.1 * B[0]))
    with pytest.raises(ValueError):
        WfcParser.overlap(s, dict(s, b=1.1 * B))

def test_transform():
    s = cubic_states()
    C4z = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    tau = np.array([0.3, -0.7, 1.1])
    for op in [dict(rotation=C4z), dict(translation=tau), dict(time_reversal=True),
               dict(rotation=C4z, translation=tau, time_reversal=True)]:
        t = WfcParser.transform(s, **op)
        assert np.allclose(WfcParser.overlap(t, t), np.eye(3))       # the operations are unitary
    # rotation: k' = Rk and the coefficient of G is moved to RG
    t = WfcParser.transform(s, rotation=C4z)
    assert np.allclose(t['k'], C4z @ s['k'])
    i = np.where((s['mill'] == [1, 0, 0]).all(1))[0][0]
    j = np.where((t['mill'] == [0, 1, 0]).all(1))[0][0]
    assert np.allclose(t['c'][:, j], s['c'][:, i])
    # translation: phase exp(-i(k+G).tau), and the inverse translation gives back the states
    t = WfcParser.transform(s, translation=tau)
    assert np.allclose(t['c'], s['c'] * np.exp(-1j * (s['k'] + s['mill'] @ B) @ tau))
    assert np.allclose(WfcParser.transform(t, translation=-tau)['c'], s['c'])
    # time reversal: k -> -k, c(G) -> c*(-G); applied twice it gives back the states
    t = WfcParser.transform(s, time_reversal=True)
    assert np.allclose(t['k'], -s['k']) and np.array_equal(t['mill'], -s['mill']) and np.allclose(t['c'], s['c'].conj())
    tt = WfcParser.transform(t, time_reversal=True)
    assert np.allclose(WfcParser.overlap(s, tt), np.eye(3))
    # a rotation that does not map the lattice onto itself
    c, sn = np.cos(np.pi / 6), np.sin(np.pi / 6)
    with pytest.raises(ValueError):
        WfcParser.transform(s, rotation=[[c, -sn, 0.], [sn, c, 0.], [0., 0., 1.]])

def test_pseudo_type(tmp_path):
    upf2 = tmp_path / 'X.UPF'
    upf2.write_text('<UPF version="2.0.1">\n<PP_HEADER\n   pseudo_type="US"\n   is_ultrasoft="T"/>\n')
    assert _pseudo_type(str(upf2)) == 'US'
    # UPF version 1: the type is read from the header, not from the other sections (e.g. the GIPAW data)
    upf1 = tmp_path / 'Y.UPF'
    upf1.write_text('<PP_INFO>\n</PP_INFO>\n<PP_HEADER>\n   0   Version Number\n   C   Element\n'
                    '   NC  Norm - Conserving pseudopotential\n</PP_HEADER>\n<PP_GIPAW_RECONSTRUCTION_DATA>\n')
    assert _pseudo_type(str(upf1)) == 'NC'
