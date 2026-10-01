import os
import numpy as np
from mppi.Utilities import BandStructure, Constants, LatticeUtils

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ypp')

def test_fcc_high_sym_points():
    # the cartesian (2pi/alat) and crystal dictionaries describe the same points (L are two symmetry equivalent points)
    a = np.array([[-0.5,0,0.5],[0,0.5,0.5],[-0.5,0.5,0]]) # pw ibrav=2 lattice vectors in units of alat
    b = LatticeUtils.get_reciprocal_lattice(a,1.,rescale=True)
    for point in ['G','X','W','K','U']:
        cart = LatticeUtils.convert_to_cartesian(b,np.array(Constants.high_sym_fcc_crystal[point]))
        assert np.allclose(cart,Constants.high_sym_fcc[point])
    L = LatticeUtils.convert_to_cartesian(b,np.array(Constants.high_sym_fcc_crystal['L']))
    assert np.allclose(abs(L),Constants.high_sym_fcc['L'])
    # W lies on the square face of the Brillouin zone centered in X: W-X is orthogonal to G-X
    X,W = np.array(Constants.high_sym_fcc['X']),np.array(Constants.high_sym_fcc['W'])
    assert np.isclose(np.dot(W-X,X),0)

def test_bands_from_ypp_file():
    file = os.path.join(DATA,'o-bands.bands_interpolated')
    hs = {'L':[0.5,0.5,0.5],'G':[0.,0.,0.],'X':[0.,0.5,0.5]} # crystal coordinates (cooOut = rlu)
    bands = BandStructure.from_Ypp_file(file,high_sym_points=hs)
    assert bands.bands.shape[0] == 8 and bands.kpoints.shape[1] == 3
    assert np.all(np.diff(bands.kpath) >= 0)                 # the abscissa of ypp is used
    labels,positions = bands.get_high_sym_positions()
    assert labels == ['L',r'$\Gamma$','X']
    assert positions == sorted(positions)

def test_kpath_cartesian():
    kpoints = [[0,0,0],[1,0,0],[1,1,0]]
    bands = BandStructure(kpoints=kpoints,bands=[[0,1,2]])
    assert np.allclose(bands.kpath,[0,1,2])
