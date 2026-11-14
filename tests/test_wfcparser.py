"""
Tests of the WfcParser module on synthetic data (the tests on real wavefunctions are in the tutorial notebook, since
the outputs of the tutorials are not in git).
"""
import numpy as np
import pytest
from scipy.io import FortranFile
from scipy.integrate import trapezoid
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

# tools for 2D systems, on a hexagonal cell with the third vector along z
LAT2D = np.array([[4.6, 0., 0.], [2.3, 2.3 * np.sqrt(3), 0.], [0., 0., 20.]])
B2D = 2 * np.pi * np.linalg.inv(LAT2D).T

def states_2d(coefs, k=(0., 0., 0.), gmax=2, gzmax=4):
    # states with the given coefficients {(m1, m2, m3): c} (one dict per band) on the G vectors with |m1|, |m2| <= gmax
    # and |m3| <= gzmax
    r, rz = np.arange(-gmax, gmax + 1), np.arange(-gzmax, gzmax + 1)
    mill = np.array(np.meshgrid(r, r, rz, indexing='ij')).reshape(3, -1).T
    index = {tuple(m): i for i, m in enumerate(mill)}
    c = np.zeros((len(coefs), len(mill)), dtype=complex)
    for ib, cb in enumerate(coefs):
        for m, v in cb.items():
            c[ib, index[m]] = v
    return {'ik': 1, 'k': np.asarray(k) @ B2D, 'b': B2D, 'mill': mill, 'bands': np.arange(1, len(coefs) + 1), 'c': c}

def random_states_2d(nb=2, seed=0):
    s = states_2d([{}] * nb)
    rng = np.random.default_rng(seed)
    c = rng.normal(size=s['c'].shape) + 1j * rng.normal(size=s['c'].shape)
    return dict(s, c=c / np.linalg.norm(c, axis=1)[:, None], k=np.array([0.1, 0.2, 0.]) @ B2D)

def test_z_profile():
    c, h = LAT2D[2, 2], 1 / np.sqrt(2)
    s = states_2d([{(0, 0, 0): 1.}, {(0, 0, 0): h, (0, 0, 1): h}, {(0, 0, 0): h, (1, 0, 1): h}])
    z, p = WfcParser.z_profile(s, nz=64)
    assert z.min() == -c / 2 and z.max() < c / 2 and np.all(np.diff(z) > 0)
    assert np.allclose(p[0], 1 / c)                                    # one plane wave: flat profile
    assert np.allclose(p[1], (1 + np.cos(2 * np.pi * z / c)) / c)       # two G_z with the same G_par interfere
    assert np.allclose(p[2], 1 / c)                                    # different G_par do not interfere
    z, p = WfcParser.z_profile(random_states_2d())
    assert np.allclose(p.sum(axis=1) * (z[1] - z[0]), 1.)
    with pytest.raises(ValueError):
        WfcParser.z_profile(s, nz=5)

def test_z_weight():
    c, h = LAT2D[2, 2], 1 / np.sqrt(2)
    g = 2 * np.pi / c
    s = states_2d([{(0, 0, 0): h, (0, 0, 1): h}])
    z1, z2 = -3., 5.
    exact = (z2 - z1) / c + (np.sin(g * z2) - np.sin(g * z1)) / (g * c)
    assert np.isclose(WfcParser.z_weight(s, z1, z2)[0], exact)
    s = random_states_2d()
    assert np.allclose(WfcParser.z_weight(s, -c / 2, c / 2), 1.)
    # the same as the integral of the profile, also for a window that crosses the boundary of the cell
    z, p = WfcParser.z_profile(s, nz=4000)
    zz, pp = np.concatenate([z - c, z, z + c]), np.tile(p, 3)
    for z1, z2 in [(-2., 3.), (8., 14.)]:
        sel = (zz >= z1) & (zz <= z2)
        assert np.allclose(WfcParser.z_weight(s, z1, z2), trapezoid(pp[:, sel], zz[sel]), atol=2e-3)
    with pytest.raises(ValueError):
        WfcParser.z_weight(s, 1., -1.)

def test_to_real_space_and_inplane_map():
    s = random_states_2d()
    psi = WfcParser.to_real_space(s, grid=(6, 7, 12))
    volume = abs(np.linalg.det(LAT2D))
    dV = volume / psi[0].size
    assert np.allclose((np.abs(psi)**2).sum(axis=(1, 2, 3)) * dV, 1.)
    # value at a point of the grid vs the direct sum over the G vectors
    i, j, l = 2, 5, 9
    r = i / 6 * LAT2D[0] + j / 7 * LAT2D[1] + l / 12 * LAT2D[2]
    direct = s['c'] @ np.exp(1j * (s['k'] + s['mill'] @ B2D) @ r) / np.sqrt(volume)
    assert np.allclose(psi[:, i, j, l], direct)
    with pytest.raises(ValueError):
        WfcParser.to_real_space(s, grid=(3, 7, 12))
    # None sets the minimum number of points along a lattice vector
    assert WfcParser.to_real_space(s, grid=(8, None, None)).shape == (2, 8, 5, 9)
    # in-plane maps: over the whole cell they integrate to one, over a window to the weight in the window
    m = WfcParser.inplane_map(s, grid=(6, 7, 12))
    dA = volume / LAT2D[2, 2] / (6 * 7)
    assert np.allclose(m.sum(axis=(1, 2)) * dA, 1.)
    m = WfcParser.inplane_map(s, z_window=(-4., 4.), grid=(6, 7, 400))
    assert np.allclose(m.sum(axis=(1, 2)) * dA, WfcParser.z_weight(s, -4., 4.), atol=1e-2)

def test_layer_cell_check():
    s = random_states_2d()
    tilted = LAT2D.copy()
    tilted[2] = [1., 0., 20.]
    with pytest.raises(ValueError):
        WfcParser.z_profile(dict(s, b=2 * np.pi * np.linalg.inv(tilted).T))

def test_plots_2d():
    import matplotlib
    matplotlib.use('Agg')
    wfc = WfcParser.__new__(WfcParser)          # only the atoms are needed by the plots
    wfc.atoms = [['C', np.zeros(3)], ['C', np.array([2.3, 1.33, 0.])]]
    s = random_states_2d()
    ax = wfc.plot_z_profile(s, bands=[2])
    assert len(ax.lines) == 2                     # the profile and one plane of atoms
    ax = wfc.plot_inplane_map(s, 1, z_window=(-2., 2.), grid=(6, 7, 12), ncells=(2, 2))
    assert ax.get_title() == 'band 1'
    with pytest.raises(ValueError):
        wfc.plot_z_profile(s, bands=[5])
