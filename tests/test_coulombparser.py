"""
Tests of the YamboCoulombParser module on synthetic ndb.cutoff and ndb.RIM databases written in the tests with the
netCDF variables of Yambo.
"""
import os
import numpy as np
import pytest
from netCDF4 import Dataset
from scipy.special import exp1
from mppi.Parsers import YamboCoulombParser, YamboScreeningParser
from test_screeningparser import write_db1, write_screening, write_text, gvectors, A, C, QCRYS, FREQS

NQBZ = 16
ALAT = np.array([A, A, C])

def write_head(d, qcrys=QCRYS):
    for n, s in [('D_1', 1), ('D_2', 2), ('D_3', 3), ('D_4', 4), ('D_10', 10), ('D_100', 100), ('D_nq', len(qcrys))]:
        d.createDimension(n, s)
    d.createVariable('HEAD_QPT', 'f4', ('D_3', 'D_nq'))[:] = np.asarray(qcrys).T
    d.createVariable('HEAD_R_LATT', 'f4', ('D_4',))[:] = [len(qcrys), NQBZ, len(qcrys), NQBZ]

def write_cutoff(filename, bare_qpg, description='slab z', qcrys=QCRYS):
    """bare_qpg: array (nq, ng_cut)."""
    with Dataset(filename, 'w') as d:
        write_head(d, qcrys)
        d.createDimension('D_ng', bare_qpg.shape[1])
        write_text(d, 'GEOMETRY', description)
        write_text(d, 'DESCRIPTION', description)
        pars = np.zeros(10)
        pars[5] = bare_qpg.shape[1]
        d.createVariable('PARS', 'f4', ('D_10',))[:] = pars
        b = np.asarray(bare_qpg, dtype=complex).T
        d.createVariable('CUT_BARE_QPG', 'f4', ('D_ng', 'D_nq', 'D_2'))[:] = np.stack([b.real, b.imag], -1)

def write_rim(filename, rim_qpg, diagonal=True, rl_vol=None):
    """rim_qpg: array (nq, rim_ng, rim_ng), order [iq, G, G']."""
    rl_vol = (2 * np.pi)**3 / (A * A * C) if rl_vol is None else rl_vol
    with Dataset(filename, 'w') as d:
        write_head(d)
        d.createDimension('D_ng', rim_qpg.shape[1])
        write_text(d, 'CUTOFF', 'none')
        d.createVariable('PARS', 'f4', ('D_10',))[:] = [rim_qpg.shape[1], int(diagonal), 1000, rl_vol * 1.001, rl_vol,
                                                        0, 0, 0, 0, 0]
        d.createVariable('RIM_qpg', 'f4', ('D_ng', 'D_ng', 'D_nq'))[:] = np.transpose(rim_qpg, (2, 1, 0))

def slab_b(qcrys, iq, G, rim=True):
    """Analytic bare_qpg of the slab cutoff along z of Yambo (G in Miller indexes)."""
    lcut = C / 2
    qG = (qcrys[iq - 1] if iq > 1 else np.zeros(3)) * 2 * np.pi / ALAT + G * 2 * np.pi / ALAT
    norm = np.linalg.norm(qG, axis=1)
    g0 = np.all(G == 0, axis=1)
    if iq == 1:
        norm[g0] = 1e-5
    Qp = np.maximum(np.linalg.norm(qG[:, :2], axis=1), 1e-5)
    b = norm / np.sqrt(1 - np.exp(-Qp * lcut) * np.cos(qG[:, 2] * lcut))
    if iq == 1 and not rim:
        # average over a circle of area S: 2 pi/S * Ein(R lcut), Ein(x) = gamma + ln(x) + E1(x)
        S = (2 * np.pi)**3 / (A * A * C) / NQBZ / (2 * np.pi / C)
        x = np.sqrt(S / np.pi) * lcut
        b[g0] = 1 / np.sqrt(2 * np.pi / S * (np.euler_gamma + np.log(x) + exp1(x)))
    return b

@pytest.fixture
def save(tmp_path):
    write_db1(str(tmp_path / 'SAVE'))
    return str(tmp_path / 'SAVE')

def test_arguments(save, tmp_path):
    with pytest.raises(ValueError):
        YamboCoulombParser(save)
    with pytest.raises(IOError, match='ndb.cutoff'):
        YamboCoulombParser(save, cutoff=str(tmp_path))
    with pytest.raises(IOError, match='ndb.RIM'):
        YamboCoulombParser(save, rim=str(tmp_path / 'missing' / 'ndb.RIM'))

def test_cutoff(save, tmp_path):
    G = gvectors()
    b = 1 + np.random.default_rng(1).random((3, 30))                   # the first 30 G vectors of the SAVE
    write_cutoff(str(tmp_path / 'ndb.cutoff'), b)
    c = YamboCoulombParser(save, cutoff=str(tmp_path), verbose=False)
    assert c.geometry == 'slab z' and c.ng_cut == 30 and c.nq == 3 and c.nqbz == NQBZ and c.rim is None
    assert np.allclose(c.q_crystal, QCRYS, atol=1e-7) and np.array_equal(c.G_crystal, G)
    for iq in (1, 2, 3):
        assert np.allclose(c.bare_qpg(iq), b[iq - 1], rtol=1e-6)
        assert np.allclose(c.coulomb(iq, 10), 4 * np.pi / b[iq - 1, :10]**2, rtol=1e-6)
    with pytest.raises(ValueError):
        c.coulomb(1, 31)
    with pytest.raises(ValueError):
        c.coulomb_rim(1)

def test_no_cutoff(save, tmp_path):
    # with ndb.RIM only, the point values are computed without cutoff
    write_rim(str(tmp_path / 'ndb.RIM'), np.ones((3, 4, 4)))
    c = YamboCoulombParser(save, rim=str(tmp_path / 'ndb.RIM'), verbose=False)
    assert c.geometry == 'none' and c.ng_cut is None
    G = gvectors() * 2 * np.pi / ALAT
    q2 = QCRYS[1] * 2 * np.pi / A
    assert np.allclose(c.coulomb(2), 4 * np.pi / np.linalg.norm(q2 + G, axis=1)**2, rtol=1e-6)
    v1 = c.coulomb(1, q0_norm=2e-5)
    g0 = np.all(gvectors() == 0, axis=1)
    assert np.isclose(v1[g0][0], 4 * np.pi / 4e-10) and np.allclose(v1[~g0], 4 * np.pi / np.sum(G[~g0]**2, axis=1))

def test_rim(save, tmp_path):
    rng = np.random.default_rng(3)
    r = rng.random((3, 5, 5))                                           # not symmetric: checks the transposition
    write_rim(str(tmp_path / 'ndb.RIM'), r, diagonal=False)
    c = YamboCoulombParser(save, rim=str(tmp_path), verbose=False)
    assert c.rim_ng == 5 and not c.rim_is_diagonal and c.rim_points == 1000
    q_weight = (2 * np.pi)**3 / (A * A * C) / NQBZ / (2 * np.pi)**3
    for iq in (1, 2, 3):
        assert np.allclose(c.coulomb_rim(iq), 4 * np.pi * r[iq - 1] / (2 * q_weight), rtol=1e-6)
    write_rim(str(tmp_path / 'ndb.RIM'), r, diagonal=True)
    c = YamboCoulombParser(save, rim=str(tmp_path), verbose=False)
    assert np.allclose(c.coulomb_rim(2), 4 * np.pi * np.diag(r[1]) / (2 * q_weight), rtol=1e-6)

def test_rim_of_the_point_values(save, tmp_path):
    # a RIM_qpg equal to 2 q_weight / b^2 gives back the point values of ndb.cutoff
    G = gvectors()
    b = np.array([slab_b(QCRYS, iq, G) for iq in (1, 2, 3)])
    q_weight = (2 * np.pi)**3 / (A * A * C) / NQBZ / (2 * np.pi)**3
    r = np.array([np.diag(2 * q_weight / bi[:7]**2) for bi in b])
    write_cutoff(str(tmp_path / 'ndb.cutoff'), b)
    write_rim(str(tmp_path / 'ndb.RIM'), r)
    c = YamboCoulombParser(save, cutoff=str(tmp_path), rim=str(tmp_path), verbose=False)
    for iq in (1, 2, 3):
        assert np.allclose(c.coulomb_rim(iq), c.coulomb(iq, 7), rtol=1e-5)

def test_files_of_different_geometries(save, tmp_path):
    os.makedirs(str(tmp_path / 'a'))
    os.makedirs(str(tmp_path / 'b'))
    write_cutoff(str(tmp_path / 'a' / 'ndb.cutoff'), np.ones((3, 10)), qcrys=QCRYS * 0.5)
    write_rim(str(tmp_path / 'b' / 'ndb.RIM'), np.ones((3, 4, 4)))
    with pytest.raises(ValueError, match='do not agree'):
        YamboCoulombParser(save, cutoff=str(tmp_path / 'a'), rim=str(tmp_path / 'b'))

def test_check_slab(save, tmp_path):
    G = gvectors()
    for rim in (True, False):
        b = np.array([slab_b(QCRYS, iq, G, rim=rim) for iq in (1, 2, 3)])
        write_cutoff(str(tmp_path / 'ndb.cutoff'), b)
        c = YamboCoulombParser(save, cutoff=str(tmp_path), verbose=False)
        dev = c.check_slab(rim=rim)
        assert max(dev.values()) < 1e-6
        if not rim:
            assert c.check_slab(rim=True)['q=0,G=0'] > 0.1
    # values in a different order of the G vectors (not those of the SAVE)
    write_cutoff(str(tmp_path / 'ndb.cutoff'), b[:, ::-1])
    c = YamboCoulombParser(save, cutoff=str(tmp_path), verbose=False)
    assert c.check_slab(rim=False)['q'] > 0.1
    write_cutoff(str(tmp_path / 'ndb.cutoff'), b, description='box xyz')
    with pytest.raises(ValueError):
        YamboCoulombParser(save, cutoff=str(tmp_path), verbose=False).check_slab()

def test_planar_slab_coulomb(tmp_path):
    # screening database with X = 0 and the slab-cut v: W = v^1/2 (1 + X) v^1/2 = v, and in the mixed representation
    # v(Q; z, z') = (2 pi/Q) e^{-Q|z - z'|} for |z - z'| < c/2 (the representation is periodic in z)
    N = 600
    G = np.array([[0, 0, m] for m in range(-N, N + 1)])
    write_db1(str(tmp_path / 'SAVE'), G=G)
    zeros = {iq: np.zeros((len(FREQS), len(G), len(G))) for iq in (1, 2, 3)}
    write_screening(str(tmp_path / 'scr'), zeros, G=G, cutoff='slab z')
    with np.errstate(divide='ignore', invalid='ignore'):
        b = np.array([slab_b(QCRYS, iq, G) for iq in (1, 2, 3)])
    write_cutoff(str(tmp_path / 'ndb.cutoff'), b)
    scr = YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), verbose=False)
    coul = YamboCoulombParser(str(tmp_path / 'SAVE'), cutoff=str(tmp_path), verbose=False)
    iq = 2
    s = np.sqrt(coul.coulomb(iq, scr.ng))
    W = s[:, None] * scr.eps_inv(iq)[0] * s[None, :]
    z = np.linspace(-8.5, 9, 15)
    P = scr.planar(W, iq, (0, 0), z, 1.0).real[:, 0]
    Q = np.linalg.norm(scr.inplane_Q(iq)[:2])
    assert np.allclose(P, 2 * np.pi / Q * np.exp(-Q * np.abs(z - 1.0)), rtol=2e-3)
