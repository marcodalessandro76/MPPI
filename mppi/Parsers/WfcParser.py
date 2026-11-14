r"""
Class to read the Kohn-Sham wavefunctions of a QuantumESPRESSO computation (``<prefix>.save`` folder) or of the
``SAVE`` folder of Yambo, and to perform some basic operations on them (scalar products, symmetry transformations).

The wavefunctions are expanded in plane waves,

.. math::
    \psi_{n\mathbf{k}}(\mathbf{r}) = \frac{1}{\sqrt{\Omega}}\sum_\mathbf{G} c_n(\mathbf{G})\,e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}

with :math:`\sum_\mathbf{G}|c_n(\mathbf{G})|^2 = 1` (norm-conserving pseudopotentials), and the G vectors are identified by
their Miller indices, :math:`\mathbf{G} = m_1\mathbf{b}_1+m_2\mathbf{b}_2+m_3\mathbf{b}_3`. Only computations without spin
(no lsda, no spinors) are supported. For 2D systems (first two lattice vectors in the xy plane, third one along z) the
class computes the planar average of the density of the states along z, their weight in a window of z and their
density integrated along z in the plane. The module can be loaded in the notebook as follows

>>> from mppi.Parsers import WfcParser

>>> wfc = WfcParser.from_qe('run_dir/prefix.save')
>>> states = wfc.get_states(1, bands=[1, 2, 3, 4])

"""
import os
import re
import numpy as np
from mppi.Utilities.Constants import HaToeV

def _band_indexes(bands, nbands):
    """
    Convert the band indexes (1-based, None for all the bands) into an array of 0-based indexes and check them.
    """
    if bands is None:
        return np.arange(nbands)
    idx = np.atleast_1d(np.asarray(bands, dtype=int)) - 1
    if idx.min() < 0 or idx.max() >= nbands:
        raise ValueError('The band indexes must be between 1 and %d' % nbands)
    return idx

def read_qe_wfc(filename, bands=None):
    """
    Read the plane-wave coefficients of a k point from a ``wfc<ik>.dat`` file of QuantumESPRESSO (Fortran unformatted
    format, computations without spin).

    Args:
        filename (:py:class:`string`) : name of the file
        bands (:py:class:`list` or None) : bands to read (1-based). Default is None (all the bands)

    Returns:
        :py:class:`dict` : with the keys 'k' (cartesian, bohr^-1), 'b' (reciprocal lattice vectors as rows, bohr^-1),
            'mill' (Miller indices of the G vectors, shape (ng, 3)), 'bands' (1-based) and 'c' (coefficients, shape
            (len(bands), ng))
    """
    from scipy.io import FortranFile
    with FortranFile(filename, 'r') as f:
        head = f.read_record('i4,3f8,i4,i4,f8')            # ik, k, ispin, gamma_only, scale factor
        if head['f3'][0]:
            raise NotImplementedError('Computations with gamma_only = .true. are not supported')
        ngw, npw, npol, nbnd = f.read_ints('i4')
        if npol != 1:
            raise NotImplementedError('Computations with spinors (npol = %d) are not supported' % npol)
        b = f.read_reals('f8').reshape(3, 3)
        mill = f.read_ints('i4').reshape(npw, 3)
        idx = _band_indexes(bands, nbnd)
        c = np.zeros((len(idx), npw), dtype=complex)
        for ib in range(idx.max() + 1):
            record = f.read_reals('f8')
            for j in np.where(idx == ib)[0]:
                c[j] = record.view(complex)
    return {'k': head['f1'][0], 'b': b, 'mill': mill, 'bands': idx + 1, 'c': c}

def _pseudo_type(upf):
    """
    Type of a pseudopotential (NC, SL, US, PAW, ...) read from the header of its UPF file, None if not found.
    """
    with open(upf, errors='ignore') as f:
        text = f.read()
    start = text.find('<PP_HEADER')
    if start < 0:
        return None
    header = text[start:start + 5000]
    match = re.search(r'pseudo_type\s*=\s*"\s*(\w+)', header)       # UPF version 2
    if not match:
        match = re.search(r'^\s*(NC|SL|US|PAW)\s', header, re.M)     # UPF version 1
    return match.group(1).upper() if match else None

def _gvector_keys(mill1, mill2):
    """
    One integer for each G vector of two sets of Miller indices, equal for the same G vector.
    """
    allm = np.vstack([mill1, mill2]).astype(np.int64)
    lo = allm.min(axis=0)
    span = allm.max(axis=0) - lo + 1
    key = lambda m: ((m[:, 0] - lo[0]) * span[1] + (m[:, 1] - lo[1])) * span[2] + (m[:, 2] - lo[2])
    return key(mill1.astype(np.int64)), key(mill2.astype(np.int64))

class WfcParser():
    """
    Class to read the Kohn-Sham wavefunctions from the output of QuantumESPRESSO (:py:meth:`from_qe`) or from the
    ``SAVE`` folder of Yambo (:py:meth:`from_yambo`). The global data are read when the object is built, the
    plane-wave coefficients only when they are requested, for the chosen k point and bands (:py:meth:`get_states`).
    The k points and the bands are numbered from 1, as in QuantumESPRESSO and Yambo.

    The states of one k point are a dictionary with the keys 'ik', 'k' (cartesian, bohr^-1), 'b' (reciprocal lattice
    vectors as rows, bohr^-1), 'mill' (Miller indices of the G vectors), 'bands' and 'c' (coefficients, shape
    (number of bands, number of G vectors)). The static methods :py:meth:`overlap` and :py:meth:`transform` act on
    these dictionaries, so they can be used also with the states of different objects (e.g. of two computations).

    Args:
        save_dir (:py:class:`string`) : the ``<prefix>.save`` folder of QuantumESPRESSO or the ``SAVE`` folder of Yambo
        source (:py:class:`string`) : 'qe' or 'yambo'
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        save_dir, source : the folder and the kind of the data
        lattice (:py:class:`numpy.ndarray`) : lattice vectors as rows (bohr)
        reciprocal_lattice (:py:class:`numpy.ndarray`) : reciprocal lattice vectors as rows (bohr^-1), with
            :math:`\\mathbf{a}_i\\cdot\\mathbf{b}_j = 2\\pi\\delta_{ij}`
        volume (:py:class:`float`) : volume of the cell (bohr^3)
        atoms (:py:class:`list`) : list with the name and the cartesian position (bohr) of each atom
        kpoints (:py:class:`numpy.ndarray`) : k points in cartesian coordinates (bohr^-1)
        kpoints_crystal (:py:class:`numpy.ndarray`) : k points in crystal coordinates (units of the b_i)
        evals (:py:class:`numpy.ndarray`) : Kohn-Sham energies (eV), shape (nkpoints, nbands)
        nkpoints, nbands (:py:class:`int`) : number of k points and of bands
        num_electrons (:py:class:`float`) : number of electrons
    """

    def __init__(self, save_dir, source='qe', verbose=True):
        self.save_dir = save_dir
        self.source = source
        self.verbose = verbose
        if source == 'qe':
            self._read_qe_metadata()
        elif source == 'yambo':
            self._read_yambo_metadata()
        else:
            raise ValueError("The source must be 'qe' or 'yambo'")
        self.volume = abs(np.linalg.det(self.lattice))
        self.reciprocal_lattice = 2 * np.pi * np.linalg.inv(self.lattice).T
        self.kpoints_crystal = self.kpoints @ self.lattice.T / (2 * np.pi)
        self.nkpoints, self.nbands = self.evals.shape
        if verbose:
            self.get_info()

    @classmethod
    def from_qe(cls, save_dir, verbose=True):
        """
        Build the object from the ``<prefix>.save`` folder of a QuantumESPRESSO computation (metadata from the
        ``data-file-schema.xml`` file, coefficients from the ``wfc<ik>.dat`` files).

        Args:
            save_dir (:py:class:`string`) : the ``<prefix>.save`` folder
            verbose (:py:class:`boolean`) : define the amount of information provided on terminal
        """
        return cls(save_dir, source='qe', verbose=verbose)

    @classmethod
    def from_yambo(cls, save_dir, verbose=True):
        """
        Build the object from the ``SAVE`` folder of Yambo (metadata from the ``ns.db1`` database, coefficients from
        the ``ns.wf`` database and its fragments). Yambo stores the coefficients in single precision.

        Args:
            save_dir (:py:class:`string`) : the ``SAVE`` folder
            verbose (:py:class:`boolean`) : define the amount of information provided on terminal
        """
        return cls(save_dir, source='yambo', verbose=verbose)

    def _read_qe_metadata(self):
        """
        Read the global data from the ``data-file-schema.xml`` file of QuantumESPRESSO.
        """
        from mppi.Parsers.PwParser import PwParser
        xml = os.path.join(self.save_dir, 'data-file-schema.xml')
        if not os.path.isfile(xml):
            raise IOError('File %s not found' % xml)
        pw = PwParser(xml, verbose=False)
        if pw.lsda or pw.noncolin:
            raise NotImplementedError('Computations with spin (lsda or noncolinear) are not supported')
        if pw.data.find('output/basis_set/gamma_only').text.strip() == 'true':
            raise NotImplementedError('Computations with gamma_only = .true. are not supported')
        self.lattice = pw.lattice
        self.atoms = [[name, np.array(pos)] for name, pos in pw.atomic_positions]
        self.kpoints = pw.kpoints * 2 * np.pi / pw.alat
        self.evals = pw.evals * HaToeV
        self.num_electrons = pw.num_electrons
        for name, (mass, pseudo) in pw.atomic_species.items():
            upf = os.path.join(self.save_dir, pseudo)
            kind = _pseudo_type(upf) if os.path.isfile(upf) else None
            if kind not in ('NC', 'SL'):
                print('Warning: the pseudopotential of %s (%s) is of type %s. The scalar products of the wavefunctions '
                      'are correct only for norm-conserving pseudopotentials' % (name, pseudo, kind))

    def _read_yambo_metadata(self):
        """
        Read the global data from the ``ns.db1`` database of Yambo, through the :py:class:`YamboDftParser` class, and
        the G vectors of the wavefunctions of each k point.
        """
        from netCDF4 import Dataset
        from mppi.Parsers.YamboDftParser import YamboDftParser
        db1 = os.path.join(self.save_dir, 'ns.db1')
        if not os.path.isfile(db1):
            raise IOError('File %s not found' % db1)
        dft = YamboDftParser(db1, verbose=False)
        if dft.spin != 1:
            raise NotImplementedError('Computations with spin are not supported')
        self.lattice = dft.lattice
        self.atoms = [[name, np.array(pos)] for name, pos in dft.atomic_positions]
        self.kpoints = dft.get_kpoints(use_scalar_alat=True) * 2 * np.pi / dft.alat[0]
        self.evals = dft.evals * HaToeV
        self.num_electrons = float(dft.num_electrons)
        with Dataset(db1) as d:
            gvectors = np.array(d.variables['G-VECTORS'][:]).T * 2 * np.pi / dft.alat
            self._wfc_ng = np.array(d.variables['WFC_NG'][:], dtype=int)
            self._wfc_grid = np.array(d.variables['WFC_GRID'][:], dtype=int)
        mill = gvectors @ self.lattice.T / (2 * np.pi)
        self._mill = np.rint(mill).astype(int)
        if np.abs(mill - self._mill).max() > 1e-3:
            raise ValueError('The G vectors of %s are not integer combinations of the reciprocal lattice vectors' % db1)

    def get_info(self):
        """
        Provide information on the data read
        """
        print('Wavefunctions of %s (%s)' % (self.save_dir, 'QuantumESPRESSO' if self.source == 'qe' else 'Yambo'))
        print('number of k points', self.nkpoints, '- number of bands', self.nbands, '- number of electrons',
              self.num_electrons)
        print('atoms:', ', '.join('%s %s' % (name, np.round(pos, 4)) for name, pos in self.atoms), '(bohr)')

    def get_states(self, ik, bands=None):
        """
        Read the plane-wave coefficients of the chosen bands at the k point ik.

        Args:
            ik (:py:class:`int`) : index of the k point (1-based)
            bands (:py:class:`list` or None) : bands to read (1-based). Default is None (all the bands)

        Returns:
            :py:class:`dict` : the states, with the keys 'ik', 'k', 'b', 'mill', 'bands' and 'c' (see the
            documentation of the class)
        """
        if ik < 1 or ik > self.nkpoints:
            raise ValueError('The index of the k point must be between 1 and %d' % self.nkpoints)
        if self.source == 'qe':
            filename = os.path.join(self.save_dir, 'wfc%d.dat' % ik)
            if not os.path.isfile(filename):
                if os.path.isfile(os.path.join(self.save_dir, 'wfc%d.hdf5' % ik)):
                    raise NotImplementedError('The wavefunctions in the hdf5 format are not supported')
                raise IOError('File %s not found' % filename)
            states = read_qe_wfc(filename, bands)
            states['k'] = self.kpoints[ik - 1]
        else:
            states = self._read_yambo_states(ik, bands)
        states['b'] = self.reciprocal_lattice
        states['ik'] = ik
        return states

    def _read_yambo_states(self, ik, bands):
        """
        Read the coefficients of the k point ik from the ``ns.wf`` database of Yambo and from its fragments.
        """
        from netCDF4 import Dataset
        idx = _band_indexes(bands, self.nbands)
        ng = self._wfc_ng[ik - 1]
        mill = self._mill[self._wfc_grid[ik - 1, :ng] - 1]
        c = np.zeros((len(idx), ng), dtype=complex)
        with Dataset(os.path.join(self.save_dir, 'ns.wf')) as d:
            fragmented = int(d.variables['FRAGMENTED'][0])
            ngroups = int(d.variables['BAND_GROUPS'][1]) if 'BAND_GROUPS' in d.variables else 1
        first = 0  # first band (0-based) of the band group
        for ig in range(1, ngroups + 1):
            if first > idx.max():
                break
            name = 'WF_COMPONENTS_@_SP_POL1_K%d_BAND_GRP_%d' % (ik, ig)
            filename = os.path.join(self.save_dir, 'ns.wf_fragments_%d_%d' % (ik, ig) if fragmented else 'ns.wf')
            with Dataset(filename) as d:
                var = d.variables[name]
                nb = var.shape[0]
                sel = np.where((idx >= first) & (idx < first + nb))[0]
                if len(sel):
                    data = np.array(var[idx[sel] - first, 0, :ng, :], dtype=float)
                    c[sel] = data[..., 0] + 1j * data[..., 1]
            first += nb
        return {'k': self.kpoints[ik - 1], 'mill': mill, 'bands': idx + 1, 'c': c}

    @staticmethod
    def overlap(states_1, states_2):
        r"""
        Compute the scalar products :math:`\langle\phi_m|\psi_n\rangle = \sum_\mathbf{G}c^{(1)*}_m(\mathbf{G})c^{(2)}_n(\mathbf{G})`
        between two sets of states at the same k point (also from different computations, with the same cell). The G
        vectors are matched through their Miller indices; if the two k points differ by a reciprocal lattice vector
        :math:`\mathbf{G}_0` the Miller indices of the second set are shifted accordingly (the states at
        :math:`\mathbf{k}+\mathbf{G}_0` are the ones at :math:`\mathbf{k}`). The result is correct only for norm-conserving
        pseudopotentials.

        Args:
            states_1, states_2 (:py:class:`dict`) : states returned by :py:meth:`get_states` (or by :py:meth:`transform`)

        Returns:
            :py:class:`numpy.ndarray` : the matrix of the scalar products, with shape (bands of states_1, bands of states_2)
        """
        b1, b2 = np.asarray(states_1['b']), np.asarray(states_2['b'])
        if not np.allclose(b1, b2, rtol=1e-5, atol=1e-5 * np.abs(b1).max()):
            raise ValueError('The two sets of states have different reciprocal lattices')
        dk = (np.asarray(states_2['k']) - np.asarray(states_1['k'])) @ np.linalg.inv(b1)   # in units of the b_i
        g0 = np.rint(dk).astype(int)
        if np.abs(dk - g0).max() > 1e-4:
            raise ValueError('The two sets of states are not at the same k point')
        key1, key2 = _gvector_keys(np.asarray(states_1['mill']), np.asarray(states_2['mill']) + g0)
        _, i1, i2 = np.intersect1d(key1, key2, return_indices=True)
        return np.asarray(states_1['c'])[:, i1].conj() @ np.asarray(states_2['c'])[:, i2].T

    @staticmethod
    def transform(states, rotation=None, translation=None, time_reversal=False):
        r"""
        Apply a symmetry operation to a set of states: the result is :math:`\psi'(\mathbf{r}) = \psi(R^{-1}(\mathbf{r}-\boldsymbol{\tau}))`,
        or its complex conjugate if time_reversal is True. The rotation maps the k point and the G vectors as
        :math:`\mathbf{k}' = R\mathbf{k}` and :math:`\mathbf{G}' = R\mathbf{G}` (it must map the reciprocal lattice onto itself),
        the time reversal as :math:`\mathbf{k}' = -\mathbf{k}`, :math:`\mathbf{G}' = -\mathbf{G}` with the conjugate coefficients,
        and the translation multiplies the coefficients by :math:`e^{-i(\mathbf{k}'+\mathbf{G}')\cdot\boldsymbol{\tau}}`.

        Args:
            states (:py:class:`dict`) : states returned by :py:meth:`get_states`
            rotation (:py:class:`numpy.ndarray` or None) : orthogonal 3x3 matrix R in cartesian coordinates
            translation (:py:class:`numpy.ndarray` or None) : cartesian translation :math:`\boldsymbol{\tau}` (bohr)
            time_reversal (:py:class:`boolean`) : if True the complex conjugate of the states is taken

        Returns:
            :py:class:`dict` : the transformed states, with the same keys of the input
        """
        k, mill, c = np.asarray(states['k'], dtype=float), np.asarray(states['mill']), np.asarray(states['c'])
        b = np.asarray(states['b'])
        if rotation is not None:
            R = np.asarray(rotation, dtype=float)
            if not np.allclose(R @ R.T, np.eye(3), atol=1e-6):
                raise ValueError('The rotation matrix must be orthogonal')
            m = (mill @ b) @ R.T @ np.linalg.inv(b)
            mill = np.rint(m).astype(int)
            if np.abs(m - mill).max() > 1e-4:
                raise ValueError('The rotation does not map the reciprocal lattice onto itself')
            k = R @ k
        if time_reversal:
            k, mill, c = -k, -mill, c.conj()
        if translation is not None:
            tau = np.asarray(translation, dtype=float)
            c = c * np.exp(-1j * ((k + mill @ b) @ tau))[None, :]
        return dict(states, k=k, mill=mill, c=c)

    # tools for 2D systems: the first two lattice vectors are in the xy plane, the third one is along z

    @staticmethod
    def _layer_cell(states):
        """
        Lattice vectors (rows, bohr) of the states and length c of the cell along z. Check that the first two lattice
        vectors are in the xy plane and that the third one is along z (with a positive z component).
        """
        lattice = 2 * np.pi * np.linalg.inv(np.asarray(states['b'], dtype=float)).T
        tol = 1e-6 * np.abs(lattice).max()
        if np.abs(lattice[:2, 2]).max() > tol or np.abs(lattice[2, :2]).max() > tol or lattice[2, 2] <= 0:
            raise ValueError('The methods for 2D systems need the first two lattice vectors in the xy plane and the '
                             'third one along z')
        return lattice, lattice[2, 2]

    @staticmethod
    def _select_bands(states, bands):
        """
        States restricted to the chosen bands (1-based labels of the key 'bands'), None for all the bands.
        """
        if bands is None:
            return states
        labels = list(np.asarray(states['bands']))
        try:
            idx = [labels.index(b) for b in np.atleast_1d(bands)]
        except ValueError:
            raise ValueError('The bands %s are not in the states (bands %s)' % (bands, labels))
        return dict(states, bands=np.asarray(states['bands'])[idx], c=np.asarray(states['c'])[idx])

    @staticmethod
    def z_profile(states, nz=None):
        r"""
        Compute the planar average of the density of each state (2D systems),

        .. math::
            p_n(z) = \int_A|\psi_n|^2\,dx\,dy = \frac{1}{c}\sum_{\mathbf{G}_\parallel}\Big|\sum_{G_z}
            c_n(\mathbf{G}_\parallel,G_z)\,e^{iG_zz}\Big|^2 , \qquad \int_0^c p_n(z)\,dz = 1

        where A is the area of the cell and c its length along z. The plane waves with different in-plane components
        do not interfere, so only one Fourier transform along z is needed for each :math:`\mathbf{G}_\parallel`.

        Args:
            states (:py:class:`dict`) : states returned by :py:meth:`get_states`
            nz (:py:class:`int` or None) : number of points of the uniform grid along z. Default is None, that is
                four times the minimum number needed to represent the states (at least 200)

        Returns:
            :py:class:`tuple` : the array z (bohr) of the grid, in the interval [-c/2,c/2), and the array of the
            profiles (bohr^-1), with shape (bands, nz)
        """
        lattice, c = WfcParser._layer_cell(states)
        mill, coef = np.asarray(states['mill']), np.asarray(states['c'])
        nmin = 2 * np.abs(mill[:, 2]).max() + 1
        if nz is None:
            nz = max(4 * nmin, 200)
        if nz < nmin:
            raise ValueError('The number of points along z must be at least %d' % nmin)
        pairs, ipair = np.unique(mill[:, :2], axis=0, return_inverse=True)
        ipair = np.ravel(ipair)
        p = np.zeros((len(coef), nz))
        for ib, cb in enumerate(coef):
            A = np.zeros((len(pairs), nz), dtype=complex)
            A[ipair, mill[:, 2] % nz] = cb
            p[ib] = (np.abs(np.fft.ifft(A, axis=1) * nz)**2).sum(axis=0) / c
        z = np.arange(nz) * c / nz
        z = np.where(z >= c / 2, z - c, z)
        order = np.argsort(z)
        return z[order], p[:, order]

    @staticmethod
    def z_weight(states, z1, z2):
        r"""
        Compute the weight of each state in the window :math:`z_1 \leq z \leq z_2` (2D systems), i.e. the integral of the
        planar average :math:`p_n(z)` (see :py:meth:`z_profile`) over the window, computed analytically. The window can
        cross the boundaries of the cell (the states are periodic along z), e.g. a window centered on a layer at z = 0.

        Args:
            states (:py:class:`dict`) : states returned by :py:meth:`get_states`
            z1, z2 (:py:class:`float`) : limits of the window (bohr), with z1 <= z2

        Returns:
            :py:class:`numpy.ndarray` : the weight of each state in the window (one for the whole cell)
        """
        lattice, c = WfcParser._layer_cell(states)
        if z2 < z1:
            raise ValueError('The window must have z1 <= z2')
        mill, coef = np.asarray(states['mill']), np.asarray(states['c'])
        m3 = mill[:, 2]
        M = np.abs(m3).max()
        pairs, ipair = np.unique(mill[:, :2], axis=0, return_inverse=True)
        ipair = np.ravel(ipair)
        g = 2 * np.pi / c
        d = np.arange(-2 * M, 2 * M + 1)                       # differences m' - m of the indexes along z
        dg = np.where(d == 0, 1., d * g)
        integral = np.where(d == 0, z2 - z1, (np.exp(1j * dg * z2) - np.exp(1j * dg * z1)) / (1j * dg))
        m = np.arange(2 * M + 1)
        K = integral[m[None, :] - m[:, None] + 2 * M]          # K[m, m'] = integral of exp(i(m'-m)gz)
        weights = np.zeros(len(coef))
        for ib, cb in enumerate(coef):
            C = np.zeros((len(pairs), 2 * M + 1), dtype=complex)
            C[ipair, m3 + M] = cb
            weights[ib] = np.einsum('pm,mn,pn->', C.conj(), K, C).real / c
        return weights

    @staticmethod
    def to_real_space(states, grid=None):
        r"""
        Compute the states in real space, :math:`\psi_n(\mathbf{r}) = \Omega^{-1/2}\sum_\mathbf{G}c_n(\mathbf{G})e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}`,
        on the uniform grid :math:`\mathbf{r}_{ijl} = (i/n_1)\mathbf{a}_1+(j/n_2)\mathbf{a}_2+(l/n_3)\mathbf{a}_3` of the cell,
        with :math:`\sum|\psi_n|^2\Delta V = 1` (:math:`\Delta V = \Omega/(n_1n_2n_3)`).

        Args:
            states (:py:class:`dict`) : states returned by :py:meth:`get_states`
            grid (:py:class:`tuple` or None) : the numbers of points (n1, n2, n3) along the lattice vectors. Default is
                None, that is the minimum grid that represents the states, 2 max|m_i| + 1 points along each vector. An
                element None of the tuple sets the minimum number of points along that vector, e.g. (36, 36, None)

        Returns:
            :py:class:`numpy.ndarray` : the states (bohr^-3/2) with shape (bands, n1, n2, n3)
        """
        b = np.asarray(states['b'], dtype=float)
        lattice = 2 * np.pi * np.linalg.inv(b).T
        mill, coef = np.asarray(states['mill']), np.asarray(states['c'])
        nmin = 2 * np.abs(mill).max(axis=0) + 1
        grid = nmin if grid is None else np.array([n0 if n is None else int(n) for n, n0 in zip(grid, nmin)])
        if np.any(grid < nmin):
            raise ValueError('The grid must have at least %s points' % (tuple(nmin),))
        kc = np.asarray(states['k']) @ lattice.T / (2 * np.pi)   # crystal coordinates of k
        idx = tuple((mill % grid).T)
        psi = np.zeros((len(coef),) + tuple(grid), dtype=complex)
        for ib, cb in enumerate(coef):
            A = np.zeros(tuple(grid), dtype=complex)
            A[idx] = cb
            psi[ib] = np.fft.ifftn(A) * np.prod(grid)
        s = [np.arange(n) / n for n in grid]
        phase = np.exp(2j * np.pi * (kc[0] * s[0][:, None, None] + kc[1] * s[1][None, :, None] + kc[2] * s[2][None, None, :]))
        return psi * phase[None] / np.sqrt(abs(np.linalg.det(lattice)))

    @staticmethod
    def inplane_map(states, z_window=None, grid=None):
        r"""
        Compute the density of each state integrated along z over a window, :math:`\int|\psi_n(x,y,z)|^2dz`, on the in-plane
        points of the real space grid (2D systems, see :py:meth:`to_real_space`). The integral is the sum over the planes
        of the grid with z (taken in the interval [-c/2,c/2)) inside the window, times their distance c/n3.

        Args:
            states (:py:class:`dict`) : states returned by :py:meth:`get_states`
            z_window (:py:class:`tuple` or None) : the limits (z1, z2) of the window (bohr). Default is None (the whole cell)
            grid (:py:class:`tuple` or None) : the grid (n1, n2, n3), see :py:meth:`to_real_space`

        Returns:
            :py:class:`numpy.ndarray` : the maps (bohr^-2) with shape (bands, n1, n2). The sum of a map times the area
            of a pixel, A/(n1 n2), gives the weight of the state in the window
        """
        lattice, c = WfcParser._layer_cell(states)
        psi = WfcParser.to_real_space(states, grid)
        n3 = psi.shape[-1]
        z = np.arange(n3) * c / n3
        z = np.where(z >= c / 2, z - c, z)
        sel = np.ones(n3, dtype=bool) if z_window is None else (z >= z_window[0]) & (z <= z_window[1])
        return (np.abs(psi[..., sel])**2).sum(axis=-1) * c / n3

    def plot_z_profile(self, states, bands=None, nz=None, ax=None):
        """
        Plot the planar average of the density of the states (see :py:meth:`z_profile`) as a function of z, in the
        interval [-c/2,c/2), with the planes of the atoms marked by vertical lines.

        Args:
            states (:py:class:`dict`) : states returned by :py:meth:`get_states`
            bands (:py:class:`list` or None) : bands to plot (1-based). Default is None (all the bands of the states)
            nz (:py:class:`int` or None) : number of points of the grid along z, see :py:meth:`z_profile`
            ax (:py:class:`matplotlib.axes.Axes` or None) : axes of the plot. Default is None (a new figure)

        Returns:
            :py:class:`matplotlib.axes.Axes` : the axes of the plot
        """
        import matplotlib.pyplot as plt
        states = self._select_bands(states, bands)
        z, p = self.z_profile(states, nz)
        c = self._layer_cell(states)[1]
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 3))
        for b, pb in zip(states['bands'], p):
            ax.plot(z, pb, label='band %d' % b)
        planes = {}
        for name, pos in self.atoms:
            zat = round((pos[2] + c / 2) % c - c / 2, 3)
            planes.setdefault(zat, []).append(name)
        for zat, names in planes.items():
            ax.axvline(zat, color='gray', ls=':', lw=0.8)
            ax.text(zat, 1.01, ' '.join(sorted(set(names))), transform=ax.get_xaxis_transform(), ha='center',
                    va='bottom', fontsize=8, color='gray')
        ax.set_xlabel('z (bohr)')
        ax.set_ylabel(r'$p_n(z)$ (bohr$^{-1}$)')
        ax.legend(fontsize=8)
        return ax

    def plot_inplane_map(self, states, band, z_window=None, grid=None, ncells=(3, 3), ax=None, cmap='viridis'):
        """
        Plot the density of a state integrated along z over a window (see :py:meth:`inplane_map`) on a few cells in
        the xy plane, with the borders of one cell and the atoms (filled markers for the atoms in the window).

        Args:
            states (:py:class:`dict`) : states returned by :py:meth:`get_states`
            band (:py:class:`int`) : band to plot (1-based)
            z_window (:py:class:`tuple` or None) : the limits (z1, z2) of the window (bohr). Default is None (the whole cell)
            grid (:py:class:`tuple` or None) : the grid (n1, n2, n3), see :py:meth:`to_real_space`
            ncells (:py:class:`tuple`) : number of cells shown along the first two lattice vectors
            ax (:py:class:`matplotlib.axes.Axes` or None) : axes of the plot. Default is None (a new figure)
            cmap (:py:class:`string`) : colormap

        Returns:
            :py:class:`matplotlib.axes.Axes` : the axes of the plot
        """
        import matplotlib.pyplot as plt
        states = self._select_bands(states, band)
        lattice, c = self._layer_cell(states)
        m = self.inplane_map(states, z_window, grid)[0]
        n1, n2 = m.shape
        tiled = np.tile(m, ncells)
        s1 = (np.arange(n1 * ncells[0] + 1) - 0.5) / n1          # corners of the pixels, centered on the grid points
        s2 = (np.arange(n2 * ncells[1] + 1) - 0.5) / n2
        X = s1[:, None] * lattice[0, 0] + s2[None, :] * lattice[1, 0]
        Y = s1[:, None] * lattice[0, 1] + s2[None, :] * lattice[1, 1]
        if ax is None:
            fig, ax = plt.subplots(figsize=(5, 4))
        mesh = ax.pcolormesh(X, Y, tiled, cmap=cmap, shading='flat')
        plt.colorbar(mesh, ax=ax, label=r'$\int|\psi|^2dz$ (bohr$^{-2}$)')
        cell = np.array([[0, 0], lattice[0, :2], lattice[0, :2] + lattice[1, :2], lattice[1, :2], [0, 0]])
        ax.plot(cell[:, 0], cell[:, 1], 'w-', lw=1)
        to_crystal = np.linalg.inv(lattice[:2, :2].T)            # in-plane crystal coordinates of the atoms
        for name, pos in self.atoms:
            zat = (pos[2] + c / 2) % c - c / 2
            inside = z_window is None or z_window[0] <= zat <= z_window[1]
            s12 = to_crystal @ pos[:2]
            s12 = s12 - np.floor(s12 + 1e-6)                        # atom folded in the cell
            for i in range(ncells[0]):
                for j in range(ncells[1]):
                    r = (s12[0] + i) * lattice[0, :2] + (s12[1] + j) * lattice[1, :2]
                    ax.plot(r[0], r[1], 'o', ms=4, mec='r', mfc='r' if inside else 'none')
        ax.set_aspect('equal')
        ax.set_xlabel('x (bohr)')
        ax.set_ylabel('y (bohr)')
        ax.set_title('band %d' % states['bands'][0], fontsize=10)
        return ax
