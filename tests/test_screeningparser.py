"""
Tests of the YamboScreeningParser module on synthetic databases written in the tests with the netCDF variables of
Yambo (the tests on a real database are in the tutorial notebook, since the outputs of the tutorials are not in git).
"""
import os
import numpy as np
import pytest
from netCDF4 import Dataset
from mppi.Parsers.YamboScreeningParser import YamboScreeningParser

A, C = 5., 20.                                  # tetragonal cell (bohr): a1, a2 in the xy plane, a3 along z
LATTICE = np.diag([A, A, C])                    # rows
QCRYS = np.array([[0, 0, 0], [0.25, 0, 0], [0.25, 0.25, 0]])
Q0 = np.array([1e-5, 0, 0])
FREQS = np.array([0, 1j])                       # Ha

def gvectors(m12=1, m3=2):
    r12, r3 = np.arange(-m12, m12 + 1), np.arange(-m3, m3 + 1)
    return np.array(np.meshgrid(r12, r12, r3, indexing='ij')).reshape(3, -1).T

def write_text(d, name, text):
    v = d.createVariable(name, 'S1', ('D_1', 'D_100'))
    v[0, :] = np.array(list(text.ljust(100)), dtype='S1')

def write_db1(folder, lattice=LATTICE):
    os.makedirs(folder, exist_ok=True)
    with Dataset(os.path.join(folder, 'ns.db1'), 'w') as d:
        for n, s in [('D_3', 3), ('D_1', 1), ('D_nsym', 2), ('D_20', 20)]:
            d.createDimension(n, s)
        sym = d.createVariable('SYMMETRY', 'f4', ('D_nsym', 'D_3', 'D_3'))
        sym[:] = np.array([np.eye(3), np.diag([-1., -1., 1.])])
        d.createVariable('LATTICE_VECTORS', 'f8', ('D_3', 'D_3'))[:] = lattice.T
        d.createVariable('LATTICE_PARAMETER', 'f8', ('D_3',))[:] = np.diag(lattice)
        d.createVariable('EIGENVALUES', 'f4', ('D_1', 'D_1', 'D_3'))[:] = 0.
        d.createVariable('K-POINTS', 'f4', ('D_3', 'D_1'))[:] = 0.
        dims = np.zeros(20)
        dims[5], dims[6], dims[11], dims[14], dims[9] = 3, 1, 1, 2, 1
        d.createVariable('DIMENSIONS', 'f4', ('D_20',))[:] = dims

def write_screening(folder, X, G=None, kind='pp', cutoff='none', kernel='HARTREE', q0=Q0):
    """X: dict iq -> array (nfreq, ng, ng) [w, G, G']. G in Miller indexes (= yambo units for this cell)."""
    G = gvectors() if G is None else G
    os.makedirs(folder, exist_ok=True)
    with Dataset(os.path.join(folder, 'ndb.' + kind), 'w') as d:
        for n, s in [('D_1', 1), ('D_3', 3), ('D_5', 5), ('D_10', 10), ('D_100', 100), ('D_ng', len(G)),
                     ('D_nq', len(QCRYS))]:
            d.createDimension(n, s)
        d.createVariable('X_PARS_1', 'f4', ('D_5',))[:] = [len(G), 1, 8, 0, 0]
        p3 = np.zeros(10)
        p3[3:6] = q0
        d.createVariable('X_PARS_3', 'f4', ('D_10',))[:] = p3
        d.createVariable('X_RL_vecs', 'f4', ('D_3', 'D_ng'))[:] = G.T
        d.createVariable('HEAD_QPT', 'f4', ('D_3', 'D_nq'))[:] = QCRYS.T
        d.createVariable('FRAGMENTED', 'f4', ('D_1',))[:] = 1
        write_text(d, 'X_TDDFT_KERNEL', kernel)
        write_text(d, 'CUTOFF', cutoff)
    for iq in range(1, len(QCRYS) + 1):
        with Dataset(os.path.join(folder, 'ndb.%s_fragment_%d' % (kind, iq)), 'w') as d:
            for n, s in [('D_1', 1), ('D_2', 2), ('D_nw', len(FREQS)), ('D_ng', len(G))]:
                d.createDimension(n, s)
            d.createVariable('FREQ_PARS_sec_iq%d' % iq, 'f4', ('D_1',))[:] = 0
            d.createVariable('FREQ_sec_iq%d' % iq, 'f4', ('D_nw', 'D_2'))[:] = np.stack([FREQS.real, FREQS.imag], -1)
            x = np.transpose(X[iq], (0, 2, 1))          # Fortran X(G, G', w) seen by netCDF as [w, G', G]
            d.createVariable('X_Q_%d' % iq, 'f4', ('D_nw', 'D_ng', 'D_ng', 'D_2'))[:] = np.stack([x.real, x.imag], -1)

def write_cutoff(filename, bare_qpg, description='slab z'):
    """bare_qpg: array (nq, ng)."""
    with Dataset(filename, 'w') as d:
        for n, s in [('D_1', 1), ('D_2', 2), ('D_3', 3), ('D_100', 100), ('D_ng', bare_qpg.shape[1]),
                     ('D_nq', len(QCRYS))]:
            d.createDimension(n, s)
        d.createVariable('HEAD_QPT', 'f4', ('D_3', 'D_nq'))[:] = QCRYS.T
        write_text(d, 'DESCRIPTION', description)
        b = bare_qpg.T
        d.createVariable('CUT_BARE_QPG', 'f4', ('D_ng', 'D_nq', 'D_2'))[:] = np.stack([b.real, b.imag], -1)

def random_X(ng, seed=0):
    rng = np.random.default_rng(seed)
    X = {}
    for iq in range(1, len(QCRYS) + 1):
        M = rng.normal(size=(len(FREQS), ng, ng)) + 1j * rng.normal(size=(len(FREQS), ng, ng))
        X[iq] = (0.05 * M).astype(np.complex64).astype(complex)       # exactly representable in single precision
    return X

@pytest.fixture
def db(tmp_path):
    write_db1(str(tmp_path / 'SAVE'))
    X = random_X(len(gvectors()))
    write_screening(str(tmp_path / 'scr'), X)
    return YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), verbose=False), X, tmp_path

def test_read(db):
    s, X, _ = db
    assert s.kind == 'pp' and s.ng == 45 and s.nq == 3 and s.bands == (1, 8) and s.kernel == 'HARTREE'
    assert s.cutoff_geometry == 'none' and s.nfreqs == 2
    assert np.allclose(s.frequencies, FREQS * 27.211386, rtol=1e-6)
    assert np.array_equal(s.G_crystal, gvectors())
    assert np.allclose(s.G, gvectors() * 2 * np.pi / np.array([A, A, C]))
    assert np.allclose(s.q_crystal, QCRYS) and np.allclose(s.q0, Q0)
    for iq in (1, 2, 3):
        assert np.array_equal(s.get_X(iq), X[iq])           # order [w, G, G'] (not symmetric: checks the transposition)
    with pytest.raises(ValueError):
        s.get_X(4)
    assert s.check_frequencies()

def test_kind_detection(tmp_path):
    write_db1(str(tmp_path / 'SAVE'))
    X = random_X(45)
    with pytest.raises(IOError):
        YamboScreeningParser(str(tmp_path), str(tmp_path / 'SAVE'), verbose=False)
    write_screening(str(tmp_path / 'scr'), X, kind='em1s')
    assert YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), verbose=False).kind == 'em1s'
    write_screening(str(tmp_path / 'scr'), X, kind='pp')
    with pytest.raises(IOError):
        YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), verbose=False)

def test_derived_quantities(db):
    s, X, _ = db
    I = np.eye(s.ng)
    for iq in (1, 2):
        assert np.allclose(s.eps_inv(iq), I + X[iq])
        chi0 = s.chi0(iq)
        assert np.allclose(np.linalg.solve(I - chi0, chi0), X[iq])          # Dyson: X = (1 - chi0)^-1 chi0
        v = s.coulomb(iq)
        assert np.allclose(s.W(iq), np.diag(v) + v[:, None] * s.chi(iq) * v[None, :])
    q2 = s.q[1]
    assert np.allclose(s.coulomb(2), 4 * np.pi / np.linalg.norm(q2 + s.G, axis=1)**2)
    v1 = s.coulomb(1)
    g0 = np.all(s.G_crystal == 0, axis=1)
    assert np.isclose(v1[g0][0], 4 * np.pi / 1e-10, rtol=1e-5)
    assert np.allclose(v1[~g0], 4 * np.pi / np.linalg.norm(s.G[~g0], axis=1)**2)
    # a matrix given by the user is used instead of the one of the database
    assert np.allclose(s.eps_inv(2, X=np.zeros((s.ng, s.ng))), I)

def test_chi0_needs_rpa(tmp_path):
    write_db1(str(tmp_path / 'SAVE'))
    write_screening(str(tmp_path / 'scr'), random_X(45), kernel='ALDA')
    s = YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), verbose=False)
    with pytest.raises(ValueError):
        s.chi0(1)

def test_cutoff(tmp_path, capsys):
    write_db1(str(tmp_path / 'SAVE'))
    G = gvectors()
    write_screening(str(tmp_path / 'scr'), random_X(len(G)), G=G, cutoff='slab z')
    s = YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), verbose=False)
    assert 'Warning' in capsys.readouterr().out
    with pytest.raises(ValueError):
        s.coulomb(2)
    b = 1 + np.random.default_rng(1).random((3, len(G) + 5))           # ndb.cutoff can have more G vectors
    write_cutoff(str(tmp_path / 'ndb.cutoff'), b)
    s = YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), cutoff=str(tmp_path), verbose=False)
    for iq in (1, 2, 3):
        assert np.allclose(s.coulomb(iq), 4 * np.pi / b[iq - 1, :len(G)]**2, rtol=1e-6)

def test_transform_translation(db):
    s, X, _ = db
    tau = np.array([0.3, -1.1, 2.5])
    Xt = s.transform(X[2], 2, translation=tau)
    ph = np.exp(-1j * s.G @ tau)
    assert np.allclose(Xt, ph[:, None] * X[2] * ph.conj()[None, :])
    assert np.allclose(s.transform(X[2], 2, translation=LATTICE[0] - 2 * LATTICE[2]), X[2])   # lattice vector

def test_transform_rotations(db):
    s, X, _ = db
    idx = lambda vectors: np.array([s._gkey[tuple(m)] for m in vectors])
    mirror_z = np.diag([1., 1., -1.])                # in-plane q: R^-1 q = q, direct branch
    m = idx(s.G_crystal * [1, 1, -1])
    for iq in (1, 2, 3):
        assert np.allclose(s.transform(X[iq], iq, rotation=mirror_z), X[iq][:, m][:, :, m])
    C2z = np.diag([-1., -1., 1.])                    # R^-1 q = -q (also q0 at iq = 1): time-reversal branch
    m = idx(-(s.G_crystal * [-1, -1, 1]))
    for iq in (1, 2, 3):
        Xt = s.transform(X[iq], iq, rotation=C2z)
        assert np.allclose(Xt, np.transpose(X[iq][:, m][:, :, m], (0, 2, 1)))
        assert np.allclose(s.transform(Xt, iq, rotation=C2z), X[iq])        # C2z twice is the identity
    mirror_y = np.diag([1., -1., 1.])                # q0 along x: direct at iq = 1, q of iq = 3 not mapped onto +-q
    m = idx(s.G_crystal * [1, -1, 1])
    assert np.allclose(s.transform(X[1], 1, rotation=mirror_y), X[1][:, m][:, :, m])
    with pytest.raises(NotImplementedError):
        s.transform(X[3], 3, rotation=mirror_y)
    with pytest.raises(ValueError):                  # 45 degrees about x (q0 unchanged): not a lattice operation
        r = np.sqrt(0.5)
        s.transform(X[1], 1, rotation=[[1, 0, 0], [0, r, -r], [0, r, r]])
    with pytest.raises(ValueError):                  # not orthogonal
        s.transform(X[1], 1, rotation=2 * np.eye(3))

def test_transform_checks(tmp_path):
    write_db1(str(tmp_path / 'SAVE'))
    G = gvectors()
    G = G[G[:, 2] < 2]                               # G set not closed under the mirror z
    write_screening(str(tmp_path / 'scr'), random_X(len(G)), G=G)
    s = YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), verbose=False)
    with pytest.raises(ValueError, match='not in the G list'):
        s.transform(None, 2, rotation=np.diag([1., 1., -1.]))
    # a cut Coulomb interaction that is not invariant under the rotation
    G = gvectors()
    write_screening(str(tmp_path / 'scr2'), random_X(len(G)), G=G, cutoff='slab z')
    b = 1 + np.random.default_rng(2).random((3, len(G)))
    write_cutoff(str(tmp_path / 'ndb.cutoff'), b)
    s = YamboScreeningParser(str(tmp_path / 'scr2'), str(tmp_path / 'SAVE'), cutoff=str(tmp_path / 'ndb.cutoff'),
                             verbose=False)
    with pytest.raises(ValueError, match='invariant'):
        s.transform(None, 2, rotation=np.diag([1., 1., -1.]))

def test_write(db):
    s, X, tmp = db
    dst = str(tmp / 'new')
    Xnew = {2: 2 * X[2]}
    r = s.write(dst, Xnew)
    assert r == {'written': [2], 'copied': [1, 3], 'unchanged': []}
    n = YamboScreeningParser(dst, str(tmp / 'SAVE'), verbose=False)
    assert np.array_equal(n.get_X(2), 2 * X[2]) and np.array_equal(n.get_X(1), X[1])
    assert s.write(dst, Xnew) == {'written': [], 'copied': [], 'unchanged': [1, 2, 3]}
    assert s.write(dst, {2: 2 * X[2] + 1e-9})['unchanged'] == [1, 2, 3]       # below the tolerance
    assert s.write(dst, lambda iq: 3 * X[iq] if iq == 3 else None)['written'] == [3]
    with pytest.raises(ValueError):
        s.write(str(tmp / 'scr'), Xnew)                                         # the template folder
    with pytest.raises(ValueError):
        s.write(dst, {1: X[1][0]})                                              # wrong shape
    open(os.path.join(dst, 'ndb.RIM_W'), 'w').close()
    with pytest.raises(ValueError, match='RIM_W'):
        s.write(dst, Xnew)
    s.write(dst, Xnew, overwrite=True)
    assert not os.path.isfile(os.path.join(dst, 'ndb.RIM_W'))

def test_write_different_header(db):
    s, X, tmp = db
    other = str(tmp / 'other')
    write_screening(other, X, kind='pp', kernel='ALDA')
    with pytest.raises(ValueError, match='different header'):
        s.write(other, {})
    assert s.write(other, {}, overwrite=True)['copied'] == [1, 2, 3]
    assert YamboScreeningParser(other, str(tmp / 'SAVE'), verbose=False).kernel == 'HARTREE'

def test_inplane_blocks_and_planar(db):
    s, X, _ = db
    blocks = s.inplane_blocks()
    assert len(blocks) == 9 and all(len(v) == 5 for v in blocks.values())
    assert np.all(np.diff(s.G[blocks[(1, 0)], 2]) > 0)
    assert np.allclose(s.inplane_Q(2, (1, 0)), s.q[1] + 2 * np.pi / A * np.array([1, 0, 0]))
    # one element G = G' = (0, 0, 1): M(Q; z, z') = e^{i Gz (z - z')} / c
    M = np.zeros((2, s.ng, s.ng), dtype=complex)
    i = s._gkey[(0, 0, 1)]
    M[:, i, i] = [1, 2]
    z, z2 = np.linspace(-C / 2, C / 2, 7), np.array([0.3, 1.])
    P = s.planar(M, 2, (0, 0), z, z2)
    Gz = 2 * np.pi / C
    assert P.shape == (2, 7, 2)
    assert np.allclose(P[1], 2 * np.exp(1j * Gz * (z[:, None] - z2[None, :])) / C)
    assert np.allclose(s.planar(M, 2, (1, 0), z, z2), 0)
    with pytest.raises(ValueError):
        s.planar(M, 2, (5, 0), z, z2)

def test_planar_slab_coulomb(tmp_path):
    # slab-cut Coulomb interaction: v(Q; z, z') = (2 pi/Q) e^{-Q|z - z'|} for |z - z'| < c/2
    write_db1(str(tmp_path / 'SAVE'))
    N = 600
    G = np.array([[0, 0, m] for m in range(-N, N + 1)])
    zeros = {iq: np.zeros((len(FREQS), len(G), len(G))) for iq in (1, 2, 3)}
    write_screening(str(tmp_path / 'scr'), zeros, G=G, cutoff='slab z')
    Gc = G * 2 * np.pi / np.array([A, A, C])
    b = np.zeros((3, len(G)))
    for iq in (1, 2, 3):
        qG = QCRYS[iq - 1] * 2 * np.pi / A + Gc
        Qp = np.linalg.norm(qG[:, :2], axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            b[iq - 1] = np.linalg.norm(qG, axis=1) / np.sqrt(1 - np.exp(-Qp * C / 2) * np.cos(qG[:, 2] * C / 2))
    b[0, N] = 1.                                     # q = 0, G = 0 (not used)
    write_cutoff(str(tmp_path / 'ndb.cutoff'), b)
    s = YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), cutoff=str(tmp_path / 'ndb.cutoff'),
                             verbose=False)
    iq = 2
    W = s.W(iq)[0]                                   # X = 0: W = v
    z = np.linspace(-8.5, 9, 15)                     # |z - z'| < c/2 = 10 (the representation is periodic in z)
    P = s.planar(W, iq, (0, 0), z, 1.0).real[:, 0]
    Q = np.linalg.norm(s.inplane_Q(iq)[:2])
    exact = 2 * np.pi / Q * np.exp(-Q * np.abs(z - 1.0))
    assert np.allclose(P, exact, rtol=2e-3)

def test_wrong_save(tmp_path):
    lattice = np.array([[A, 0, 0], [0, A, 0], [1., 0, C]])
    write_db1(str(tmp_path / 'SAVE'), lattice=lattice)
    write_screening(str(tmp_path / 'scr'), random_X(45))
    with pytest.raises(ValueError):                  # G not on the reciprocal lattice of this SAVE
        YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), verbose=False)

def test_not_2d(tmp_path):
    lattice = np.array([[A, 0, 0], [0, A, 0], [1., 0, C]])
    write_db1(str(tmp_path / 'SAVE'), lattice=lattice)
    recip = 2 * np.pi * np.linalg.inv(lattice).T
    G_iku = gvectors() @ recip * np.diag(lattice) / (2 * np.pi)
    write_screening(str(tmp_path / 'scr'), random_X(45), G=G_iku)
    s = YamboScreeningParser(str(tmp_path / 'scr'), str(tmp_path / 'SAVE'), verbose=False)
    assert np.array_equal(s.G_crystal, gvectors())
    with pytest.raises(ValueError, match='2D'):
        s.inplane_blocks()
