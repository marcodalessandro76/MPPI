r"""
Class to read the screening databases of Yambo (``ndb.em1s``, ``ndb.em1d`` and ``ndb.pp`` with their fragments), to
compute the quantities derived from them, to transform them under a rotation and a translation of the system and to
write a new database with a modified screening.

The databases contain, for each q point of the IBZ and each frequency, the symmetrized matrix

.. math::
    X_{\mathbf{G}\mathbf{G}'}(\mathbf{q},\omega) = v^{1/2}_\mathbf{G}(\mathbf{q})\,\chi_{\mathbf{G}\mathbf{G}'}(\mathbf{q},\omega)\,
    v^{1/2}_{\mathbf{G}'}(\mathbf{q}) = \epsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q},\omega) - \delta_{\mathbf{G}\mathbf{G}'}

where :math:`\chi` is the reducible response function and :math:`v_\mathbf{G}(\mathbf{q}) = 4\pi/b_\mathbf{G}(\mathbf{q})^2`
is the (possibly cut) Coulomb interaction, with :math:`b_\mathbf{G}(\mathbf{q})` the ``bare_qpg`` of Yambo
(:math:`|\mathbf{q}+\mathbf{G}|` without cutoff). The cut Coulomb interaction is read from the ``ndb.cutoff`` database
written by the run that computed it.

The module can be loaded in the notebook as follows

>>> from mppi.Parsers import YamboScreeningParser

>>> scr = YamboScreeningParser('run_dir/gw_ppa', 'run_dir/SAVE')
>>> X = scr.get_X(1)

"""
import os
import shutil
import filecmp
import numpy as np
from mppi.Utilities.Constants import HaToeV

KINDS = ('em1s', 'em1d', 'pp')

def _read_text(database, name):
    """
    Read a character variable of a netCDF database as a string.
    """
    return b''.join(np.ma.filled(database[name][:], b' ').ravel()).decode(errors='ignore').strip()

def _complex(x):
    """
    Complex array from an array whose last dimension contains the real and the imaginary parts.
    """
    x = np.asarray(x, dtype=float)
    return x[..., 0] + 1j * x[..., 1]

class YamboScreeningParser():
    r"""
    Class to read a screening database of Yambo (``ndb.em1s``, ``ndb.em1d`` or ``ndb.pp``) and to work on its matrices
    :math:`X = v^{1/2}\chi v^{1/2} = \epsilon^{-1} - 1`. One object describes one database. The global data are read
    when the object is built, the matrices X only when they are requested, one q point at a time (:py:meth:`get_X`).
    The q points are numbered from 1, as in Yambo (``ndb.pp_fragment_<iq>``).

    The matrices are returned as complex arrays with shape (number of frequencies, ng, ng) and order [w, G, G']. The
    methods that compute derived quantities accept a matrix X as optional argument, so they can be used also on a
    matrix built by the user (e.g. the X of a system moved by :py:meth:`transform`), which can be written in a new
    database with :py:meth:`write`.

    The first q point is q = 0. Yambo computes it at the small vector q0 (the 'Field direction' of the database, with
    norm 1e-5): q0 is used for the Coulomb interaction of G = 0 and its direction for the symmetry transformations.

    Args:
        folder (:py:class:`string`) : folder with the database (``ndb.em1s``, ``ndb.em1d`` or ``ndb.pp`` and its
            fragments). It must contain only one of the three databases
        save (:py:class:`string`) : the ``SAVE`` folder of the computation (or its ``ns.db1`` file), needed for the
            lattice and the lattice parameters of Yambo, which convert the G vectors and the q points
        cutoff (:py:class:`string`) : the ``ndb.cutoff`` database with the cut Coulomb interaction used by the run (or
            the folder that contains it). Needed for the methods that use v if the database has a Coulomb cutoff.
            Without cutoff v is computed as :math:`4\pi/|\mathbf{q}+\mathbf{G}|^2`
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        folder, kind : the folder and the kind of the database ('em1s', 'em1d' or 'pp')
        ng (:py:class:`int`) : number of G vectors of the matrices
        bands (:py:class:`tuple`) : first and last band used in the response function
        kernel (:py:class:`string`) : the TDDFT kernel ('HARTREE' is the RPA)
        cutoff_geometry (:py:class:`string`) : the Coulomb cutoff of the run ('none', 'slab z', ...)
        lattice (:py:class:`numpy.ndarray`) : lattice vectors as rows (bohr)
        reciprocal_lattice (:py:class:`numpy.ndarray`) : reciprocal lattice vectors as rows (bohr^-1)
        volume (:py:class:`float`) : volume of the cell (bohr^3)
        alat (:py:class:`numpy.ndarray`) : the lattice parameters of Yambo (alat_x, alat_y, alat_z), bohr
        syms (:py:class:`numpy.ndarray`) : the symmetry operations of the ``SAVE`` (cartesian matrices applied to the
            k points, as in :py:class:`YamboDftParser`), which can be passed to :py:meth:`transform`
        time_reversal (:py:class:`bool`) : True if the ``SAVE`` uses the time reversal symmetry
        G (:py:class:`numpy.ndarray`) : G vectors in cartesian coordinates (bohr^-1), shape (ng, 3)
        G_crystal (:py:class:`numpy.ndarray`) : Miller indices of the G vectors (integers), shape (ng, 3)
        q (:py:class:`numpy.ndarray`) : q points of the IBZ in cartesian coordinates (bohr^-1), shape (nq, 3)
        q_crystal (:py:class:`numpy.ndarray`) : q points in crystal coordinates (units of the b_i)
        nq (:py:class:`int`) : number of q points
        q0 (:py:class:`numpy.ndarray`) : the small vector used for q -> 0 (bohr^-1)
        frequencies (:py:class:`numpy.ndarray`) : the frequencies of the database (complex, eV), read from the first
            fragment
        nfreqs (:py:class:`int`) : number of frequencies
    """

    def __init__(self, folder, save, cutoff=None, verbose=True):
        from netCDF4 import Dataset
        from mppi.Parsers.YamboDftParser import YamboDftParser
        self.folder = folder
        self.verbose = verbose
        found = [k for k in KINDS if os.path.isfile(os.path.join(folder, 'ndb.' + k))]
        if len(found) != 1:
            raise IOError('The folder %s must contain one of the databases ndb.em1s, ndb.em1d, ndb.pp (found: %s)'
                          % (folder, ', '.join('ndb.' + k for k in found) if found else 'none'))
        self.kind = found[0]
        self.header = os.path.join(folder, 'ndb.' + self.kind)

        db1 = os.path.join(save, 'ns.db1') if os.path.isdir(save) else save
        if not os.path.isfile(db1):
            raise IOError('File %s not found' % db1)
        dft = YamboDftParser(db1, verbose=False)
        self.lattice = np.array(dft.lattice, dtype=float)
        self.alat = np.array(dft.alat, dtype=float)
        self.syms = np.array(dft.syms, dtype=float)
        self.time_reversal = dft.time_reversal
        self.volume = abs(np.linalg.det(self.lattice))
        self.reciprocal_lattice = 2 * np.pi * np.linalg.inv(self.lattice).T

        with Dataset(self.header) as d:
            pars = np.array(d['X_PARS_1'][:], dtype=float)
            self.ng = int(pars[0])
            self.bands = (int(pars[1]), int(pars[2]))
            self.kernel = _read_text(d, 'X_TDDFT_KERNEL')
            self.cutoff_geometry = _read_text(d, 'CUTOFF') if 'CUTOFF' in d.variables else 'none'
            self.q0 = np.array(np.ma.filled(d['X_PARS_3'][3:6], 0.), dtype=float)
            G_iku = np.array(d['X_RL_vecs'][:], dtype=float).T
            q_iku = np.array(d['HEAD_QPT'][:], dtype=float).T
        self.G = G_iku * 2 * np.pi / self.alat
        self.q = q_iku * 2 * np.pi / self.alat
        self.nq = len(self.q)
        G_crys = self.G @ self.lattice.T / (2 * np.pi)
        self.G_crystal = np.rint(G_crys).astype(int)
        if np.abs(G_crys - self.G_crystal).max() > 1e-3:
            raise ValueError('The G vectors of the database are not vectors of the reciprocal lattice of %s: wrong SAVE?'
                             % db1)
        self.q_crystal = self.q @ self.lattice.T / (2 * np.pi)
        self._gkey = {tuple(m): i for i, m in enumerate(self.G_crystal)}
        self.frequencies = self.get_frequencies(1)
        self.nfreqs = len(self.frequencies)

        self.cutoff = None
        self._bare_qpg = None
        if cutoff is not None:
            self._read_cutoff(cutoff)
        elif self.cutoff_geometry.lower() != 'none':
            print('Warning: the database uses the Coulomb cutoff "%s" but no ndb.cutoff is given: the methods that '
                  'need v are not available' % self.cutoff_geometry)
        if verbose:
            self.get_info()

    def get_info(self):
        """
        Provide information on the database
        """
        print('Screening database %s' % self.header)
        print('number of q points', self.nq, '- number of G vectors', self.ng, '- bands %d-%d' % self.bands,
              '- kernel', self.kernel, '- cutoff', self.cutoff_geometry)
        print('frequencies (eV):', np.round(self.frequencies, 4))

    def _xfile(self, iq, folder=None):
        """
        File of the folder (default the folder of the database) that contains the matrix X of the q point iq.
        """
        folder = self.folder if folder is None else folder
        fragment = os.path.join(folder, 'ndb.%s_fragment_%d' % (self.kind, iq))
        return fragment if os.path.isfile(fragment) else os.path.join(folder, 'ndb.' + self.kind)

    def _check_iq(self, iq):
        if not (1 <= int(iq) <= self.nq):
            raise ValueError('The index of the q point must be between 1 and %d' % self.nq)
        return int(iq)

    def get_frequencies(self, iq):
        """
        Frequencies of the q point iq (complex, eV).

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)
        """
        from netCDF4 import Dataset
        iq = int(iq)
        with Dataset(self._xfile(iq)) as d:
            return _complex(d['FREQ_sec_iq%d' % iq][:]) * HaToeV

    def check_frequencies(self):
        """
        Check that all the q points have the same frequencies (of the first one).

        Returns:
            :py:class:`bool` : True if all the q points have the same frequencies
        """
        same = all(np.allclose(self.get_frequencies(iq), self.frequencies, rtol=1e-6, atol=1e-6)
                   for iq in range(1, self.nq + 1))
        if self.verbose:
            print('All the q points have the same frequencies' if same else 'The frequencies depend on the q point')
        return same

    def get_X(self, iq):
        """
        Read the matrix X = v^1/2 chi v^1/2 = eps^-1 - 1 of the q point iq (the netCDF file is closed after reading).

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)

        Returns:
            :py:class:`numpy.ndarray` : complex array with shape (nfreqs, ng, ng), order [w, G, G']
        """
        from netCDF4 import Dataset
        iq = self._check_iq(iq)
        with Dataset(self._xfile(iq)) as d:
            x = np.array(d['X_Q_%d' % iq][:])
        # the Fortran array X(G, G', w) is seen by netCDF with the order [w, G', G]
        return np.transpose(_complex(x), (0, 2, 1))

    def _get(self, iq, X):
        """
        The matrix X given by the user (checked) or the one of the database.
        """
        if X is None:
            return self.get_X(iq)
        X = np.asarray(X)
        if X.shape[-2:] != (self.ng, self.ng):
            raise ValueError('The matrix X must have shape (..., %d, %d)' % (self.ng, self.ng))
        return X

    def eps_inv(self, iq, X=None):
        """
        Inverse dielectric matrix (symmetrized), eps^-1 = 1 + X.

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)
            X (:py:class:`numpy.ndarray`) : the matrix X to use (shape (..., ng, ng)). Default None (the X of the
                database)

        Returns:
            :py:class:`numpy.ndarray` : array with the shape of X
        """
        return np.eye(self.ng) + self._get(iq, X)

    def chi0(self, iq, X=None):
        r"""
        Symmetrized independent particle response function :math:`\tilde\chi_0 = v^{1/2}\chi_0 v^{1/2} = 1 - (1+X)^{-1}`,
        obtained by inverting the Dyson equation of the RPA in the block of the database. Available only for the
        databases computed with the 'HARTREE' kernel (RPA).

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)
            X (:py:class:`numpy.ndarray`) : the matrix X to use (shape (..., ng, ng)). Default None (the X of the
                database)

        Returns:
            :py:class:`numpy.ndarray` : array with the shape of X
        """
        if self.kernel.upper() != 'HARTREE':
            raise ValueError('chi0 can be recovered only for the RPA (HARTREE kernel), the kernel is %s' % self.kernel)
        I = np.eye(self.ng)
        return I - np.linalg.inv(I + self._get(iq, X))

    def _read_cutoff(self, cutoff):
        """
        Read the cut Coulomb interaction (bare_qpg) from the ndb.cutoff database.
        """
        from netCDF4 import Dataset
        filename = os.path.join(cutoff, 'ndb.cutoff') if os.path.isdir(cutoff) else cutoff
        if not os.path.isfile(filename):
            raise IOError('File %s not found' % filename)
        with Dataset(filename) as d:
            description = _read_text(d, 'DESCRIPTION') if 'DESCRIPTION' in d.variables else ''
            q_iku = np.array(d['HEAD_QPT'][:], dtype=float).T
            b = _complex(d['CUT_BARE_QPG'][:])          # shape (ng_cut, nq)
        if len(q_iku) != self.nq or not np.allclose(q_iku * 2 * np.pi / self.alat, self.q, atol=1e-5):
            raise ValueError('The q points of %s are not those of the screening database' % filename)
        if b.shape[0] < self.ng:
            raise ValueError('%s contains only %d G vectors, the screening database %d' % (filename, b.shape[0], self.ng))
        if description.lower() != self.cutoff_geometry.lower():
            print('Warning: the cutoff of %s is "%s", the screening database has "%s"' % (filename, description,
                                                                                          self.cutoff_geometry))
        self.cutoff = filename
        bqpg = b[:self.ng].T
        self._bare_qpg = bqpg.real if np.abs(bqpg.imag).max() == 0 else bqpg

    def bare_qpg(self, iq):
        r"""
        The ``bare_qpg`` of Yambo for the G vectors of the database, :math:`b_\mathbf{G}(\mathbf{q})` with
        :math:`v_\mathbf{G}(\mathbf{q}) = 4\pi/b_\mathbf{G}(\mathbf{q})^2`: read from ndb.cutoff, or
        :math:`|\mathbf{q}+\mathbf{G}|` without cutoff (at q = 0, :math:`|\mathbf{q}_0|` for G = 0, as in Yambo).

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)

        Returns:
            :py:class:`numpy.ndarray` : array with shape (ng,) (bohr^-1)
        """
        iq = self._check_iq(iq)
        if self._bare_qpg is not None:
            return self._bare_qpg[iq - 1]
        if self.cutoff_geometry.lower() != 'none':
            raise ValueError('The database uses the Coulomb cutoff "%s": give its ndb.cutoff to the constructor'
                             % self.cutoff_geometry)
        b = np.linalg.norm(self.q[iq - 1] + self.G, axis=1)
        if iq == 1:
            b[np.all(self.G_crystal == 0, axis=1)] = np.linalg.norm(self.q0)
        return b

    def coulomb(self, iq):
        r"""
        Coulomb interaction :math:`v_\mathbf{G}(\mathbf{q})` (cut if the run used a cutoff) for the G vectors of the
        database (Ha bohr^3, i.e. :math:`4\pi/|\mathbf{q}+\mathbf{G}|^2` without cutoff).

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)

        Returns:
            :py:class:`numpy.ndarray` : array with shape (ng,)
        """
        return 4 * np.pi / self.bare_qpg(iq)**2

    def W(self, iq, X=None):
        r"""
        Screened interaction :math:`W = v^{1/2}(1 + X)v^{1/2}` of the system described by the database (atomic units,
        Ha bohr^3, same normalization as v).

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)
            X (:py:class:`numpy.ndarray`) : the matrix X to use (shape (..., ng, ng)). Default None (the X of the
                database)

        Returns:
            :py:class:`numpy.ndarray` : array with the shape of X
        """
        s = np.sqrt(4 * np.pi) / self.bare_qpg(iq)
        return s[:, None] * self.eps_inv(iq, X) * s[None, :]

    def chi(self, iq, X=None):
        r"""
        Reducible response function :math:`\chi = v^{-1/2} X v^{-1/2}` (atomic units, bohr^-3 Ha^-1).

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)
            X (:py:class:`numpy.ndarray`) : the matrix X to use (shape (..., ng, ng)). Default None (the X of the
                database)

        Returns:
            :py:class:`numpy.ndarray` : array with the shape of X
        """
        s = self.bare_qpg(iq) / np.sqrt(4 * np.pi)
        return s[:, None] * self._get(iq, X) * s[None, :]

    def G_indexes(self, vectors, name='vectors'):
        """
        Index in the G list of the database of each of the given cartesian vectors (bohr^-1). Raise a ValueError if a
        vector is not in the list.

        Args:
            vectors (:py:class:`numpy.ndarray`) : cartesian vectors, shape (n, 3)
            name (:py:class:`string`) : name of the vectors used in the error message

        Returns:
            :py:class:`numpy.ndarray` : the indexes (0-based)
        """
        m = np.atleast_2d(vectors) @ self.lattice.T / (2 * np.pi)
        mi = np.rint(m).astype(int)
        bad = np.abs(m - mi).max(axis=1) > 1e-4
        if bad.any():
            raise ValueError('The %s %s are not vectors of the reciprocal lattice' % (name, mi[bad][:3].tolist()))
        idx = np.array([self._gkey.get(tuple(x), -1) for x in mi])
        if (idx < 0).any():
            raise ValueError('The %s with Miller indexes %s are not in the G list of the database'
                             % (name, mi[idx < 0][:5].tolist()))
        return idx

    def transform(self, X, iq, rotation=None, translation=None):
        r"""
        Matrix X of the system moved by :math:`\mathbf{r}\to R\mathbf{r}+\boldsymbol{\tau}`, i.e.
        :math:`\chi'(\mathbf{r},\mathbf{r}') = \chi(R^{-1}(\mathbf{r}-\boldsymbol{\tau}), R^{-1}(\mathbf{r}'-\boldsymbol{\tau}))`,
        which gives

        .. math::
            X'_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) = e^{-i(\mathbf{G}-\mathbf{G}')\cdot\boldsymbol{\tau}}\,
            X_{R^{-1}\mathbf{G},R^{-1}\mathbf{G}'}(R^{-1}\mathbf{q})

        Only the operations that map q onto itself or onto -q are implemented (the X outside the IBZ is not
        reconstructed): if :math:`R^{-1}\mathbf{q} = -\mathbf{q}` the time reversal (Onsager) relation
        :math:`X_{\mathbf{G}\mathbf{G}'}(-\mathbf{q}) = X_{-\mathbf{G}',-\mathbf{G}}(\mathbf{q})` is used. At q = 0 the
        direction of q0 is used. The G set of the database must be closed under R, and R must leave the Coulomb
        interaction invariant (e.g. R z = +-z for a slab cutoff along z, checked on v when it is available): then the
        symmetrized X transforms as chi.

        Args:
            X (:py:class:`numpy.ndarray`) : the matrix X of the q point iq, shape (..., ng, ng)
            iq (:py:class:`int`) : index of the q point (1-based)
            rotation (:py:class:`numpy.ndarray`) : orthogonal 3x3 matrix R in cartesian coordinates. Default None (identity)
            translation (:py:class:`numpy.ndarray`) : the vector tau in cartesian coordinates (bohr). Default None (zero)

        Returns:
            :py:class:`numpy.ndarray` : the transformed matrix, with the shape of X
        """
        iq = self._check_iq(iq)
        X = self._get(iq, X)
        R = np.eye(3) if rotation is None else np.asarray(rotation, dtype=float)
        tau = np.zeros(3) if translation is None else np.asarray(translation, dtype=float)
        if not np.allclose(R @ R.T, np.eye(3), atol=1e-6):
            raise ValueError('The rotation must be an orthogonal matrix')
        Rinv = R.T
        if iq == 1:
            if not np.linalg.norm(self.q0) > 0:
                raise ValueError('The direction q0 of the database is not defined')
            qv = self.q0 / np.linalg.norm(self.q0)
            same, opposite = np.allclose(Rinv @ qv, qv, atol=1e-6), np.allclose(Rinv @ qv, -qv, atol=1e-6)
        else:
            qc = self.q_crystal[iq - 1]
            to_crys = lambda v: v @ self.lattice.T / (2 * np.pi)
            same = np.allclose(to_crys(Rinv @ self.q[iq - 1]), qc, atol=1e-4)
            opposite = np.allclose(to_crys(Rinv @ self.q[iq - 1]), -qc, atol=1e-4)
        if same:
            m = self.G_indexes(self.G @ Rinv.T, name='rotated G vectors R^-1 G')
            Xm = X[..., m, :][..., :, m]
        elif opposite:
            m = self.G_indexes(-self.G @ Rinv.T, name='rotated G vectors -R^-1 G')
            Xm = np.swapaxes(X[..., m, :][..., :, m], -1, -2)
        else:
            raise NotImplementedError('R^-1 q is not q or -q at the q point %d: the X outside the IBZ is needed' % iq)
        if self._bare_qpg is not None or self.cutoff_geometry.lower() == 'none':
            b = self.bare_qpg(iq)
            if np.abs(b[m] - b).max() > 1e-4 * np.abs(b).max():
                raise ValueError('The rotation does not leave the Coulomb interaction invariant (cutoff "%s")'
                                 % self.cutoff_geometry)
        elif self.verbose:
            print('Warning: the invariance of the cut Coulomb interaction under the rotation is not checked (no '
                  'ndb.cutoff)')
        ph = np.exp(-1j * self.G @ tau)
        return ph[:, None] * Xm * ph.conj()[None, :]

    def write(self, dst, X, overwrite=False, tol=1e-7):
        """
        Write a new database in the folder dst, with this database as template: copy the header and all the fragments
        (only the ``ndb.<kind>*`` files, never ``ndb.RIM_W``, which Yambo builds from X and reuses without checking it)
        and replace the matrices X of the given q points. A fragment of dst is rewritten only if its X differs by more
        than tol from the new one, and it is read back and checked. The header (bands, number of G vectors, cutoff,
        frequencies) stays the one of the template: Yambo then assumes the parameters of the template.

        Args:
            dst (:py:class:`string`) : the destination folder (created if needed). It cannot be the folder of the
                template
            X (:py:class:`dict` or callable) : dictionary iq -> matrix X (shape (nfreqs, ng, ng), order [w, G, G']) or
                a function of iq that returns the matrix (or None to keep the X of the template). The q points
                missing in the dictionary keep the X of the template
            overwrite (:py:class:`boolean`) : allow writing in a folder that contains a database with a different
                header or a ``ndb.RIM_W`` file (which is removed)
            tol (:py:class:`float`) : tolerance on the maximum difference of X to rewrite a fragment that already
                exists in dst

        Returns:
            :py:class:`dict` : with the keys 'written' (q points with the new X written), 'copied' (fragments copied
            from the template) and 'unchanged' (fragments already up to date in dst), the lists of the q points
        """
        from netCDF4 import Dataset
        if self._xfile(1) == self.header:
            raise NotImplementedError('Only fragmented databases can be written')
        if os.path.abspath(dst) == os.path.abspath(self.folder):
            raise ValueError('The destination folder cannot be the folder of the template database')
        if os.path.isfile(dst):
            raise ValueError('%s is a file' % dst)
        os.makedirs(dst, exist_ok=True)
        header = os.path.join(dst, os.path.basename(self.header))
        others = [k for k in KINDS if k != self.kind and os.path.isfile(os.path.join(dst, 'ndb.' + k))]
        if others:
            raise ValueError('The folder %s contains another screening database (%s)' % (dst, ', '.join(others)))
        rim = os.path.join(dst, 'ndb.RIM_W')
        if os.path.isfile(rim):
            if not overwrite:
                raise ValueError('The folder %s contains a ndb.RIM_W, built from the old X: use overwrite=True to '
                                 'remove it' % dst)
            os.remove(rim)
        if os.path.isfile(header) and not filecmp.cmp(header, self.header, shallow=False):
            if not overwrite:
                raise ValueError('The folder %s contains a database with a different header: use overwrite=True'
                                 % dst)
            prefix = 'ndb.%s_fragment_' % self.kind
            for f in os.listdir(dst):
                if f.startswith(prefix):
                    os.remove(os.path.join(dst, f))
        if not os.path.isfile(header) or not filecmp.cmp(header, self.header, shallow=False):
            shutil.copy2(self.header, header)

        get = X if callable(X) else (lambda iq: X.get(iq))
        written, copied, unchanged = [], [], []
        for iq in range(1, self.nq + 1):
            src = self._xfile(iq)
            out = os.path.join(dst, os.path.basename(src))
            Xnew = get(iq)
            if Xnew is None:
                if not os.path.isfile(out) or not filecmp.cmp(src, out, shallow=False):
                    shutil.copy2(src, out)
                    copied.append(iq)
                else:
                    unchanged.append(iq)
                continue
            Xnew = np.asarray(Xnew)
            if Xnew.shape != (self.nfreqs, self.ng, self.ng):
                raise ValueError('The matrix X of the q point %d must have shape (%d, %d, %d)'
                                 % (iq, self.nfreqs, self.ng, self.ng))
            if os.path.isfile(out):
                try:
                    with Dataset(out) as d:
                        old = np.transpose(_complex(d['X_Q_%d' % iq][:]), (0, 2, 1))
                    if old.shape == Xnew.shape and np.abs(old - Xnew).max() <= tol:
                        unchanged.append(iq)
                        continue
                except (OSError, KeyError, IndexError):
                    pass
            shutil.copy2(src, out)
            x = np.transpose(Xnew, (0, 2, 1))
            with Dataset(out, 'r+') as d:
                d['X_Q_%d' % iq][:] = np.stack([x.real, x.imag], axis=-1).astype(np.float32)
            with Dataset(out) as d:
                back = np.transpose(_complex(d['X_Q_%d' % iq][:]), (0, 2, 1))
            if np.abs(back - Xnew).max() > 1e-6 * max(1., np.abs(Xnew).max()):
                raise IOError('The matrix X written in %s differs from the given one' % out)
            written.append(iq)
        if self.verbose:
            print('Database written in %s: X written at %d q points, %d fragments copied from the template, %d '
                  'unchanged' % (dst, len(written), len(copied), len(unchanged)))
        return {'written': written, 'copied': copied, 'unchanged': unchanged}

    def _layer_c(self):
        """
        Length c of the cell along z. Check that the first two lattice vectors are in the xy plane and that the third
        one is along z (with a positive z component).
        """
        L = self.lattice
        tol = 1e-6 * np.abs(L).max()
        if np.abs(L[:2, 2]).max() > tol or np.abs(L[2, :2]).max() > tol or L[2, 2] <= 0:
            raise ValueError('The methods for 2D systems need the first two lattice vectors in the xy plane and the '
                             'third one along z')
        return L[2, 2]

    def inplane_blocks(self):
        """
        Group the G vectors of the database by their in-plane part (2D systems: first two lattice vectors in the xy
        plane, third one along z). The matrices of a 2D system are block diagonal only in the limit of an isolated
        layer, but the blocks G_par = G_par' are the ones used by :py:meth:`planar`.

        Returns:
            :py:class:`dict` : the keys are the in-plane Miller indexes (m1, m2) of G_par, the values the indexes
            (0-based) of the G vectors with that in-plane part, sorted by G_z
        """
        self._layer_c()
        blocks = {}
        for i, m in enumerate(self.G_crystal):
            blocks.setdefault((int(m[0]), int(m[1])), []).append(i)
        return {k: np.array(sorted(v, key=lambda i: self.G[i, 2])) for k, v in blocks.items()}

    def inplane_Q(self, iq, G_par=(0, 0)):
        """
        In-plane vector Q = q + G_par (cartesian, bohr^-1) of a block of :py:meth:`inplane_blocks` at the q point iq.

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)
            G_par (:py:class:`tuple`) : in-plane Miller indexes (m1, m2) of G_par

        Returns:
            :py:class:`numpy.ndarray` : the vector Q (its z component is the one of q, zero for a 2D sampling)
        """
        iq = self._check_iq(iq)
        return self.q[iq - 1] + G_par[0] * self.reciprocal_lattice[0] + G_par[1] * self.reciprocal_lattice[1]

    def planar(self, M, iq, G_par, z, z2):
        r"""
        Mixed representation of a matrix M (e.g. X, W or chi) of a 2D system on the block G_par = G_par':

        .. math::
            M(\mathbf{Q}; z, z') = \frac{1}{c}\sum_{G_z, G_z'} e^{iG_z z} M_{\mathbf{G}\mathbf{G}'} e^{-iG_z' z'}

        with :math:`\mathbf{Q} = \mathbf{q} + \mathbf{G}_\parallel` (see :py:meth:`inplane_Q`), so that
        :math:`M(\mathbf{r},\mathbf{r}') = \frac{1}{A}\sum_\mathbf{Q} e^{i\mathbf{Q}\cdot(\boldsymbol{\rho}-\boldsymbol{\rho}')}
        M(\mathbf{Q}; z, z')` for the convention
        :math:`M(\mathbf{r},\mathbf{r}') = \frac{1}{V}\sum_\mathbf{q}\sum_{\mathbf{G}\mathbf{G}'}e^{i(\mathbf{q}+\mathbf{G})\cdot\mathbf{r}}
        M_{\mathbf{G}\mathbf{G}'}e^{-i(\mathbf{q}+\mathbf{G}')\cdot\mathbf{r}'}` (A area of the cell, V = A c,
        the off-diagonal blocks of the in-plane part are neglected). E.g. for the slab-cut Coulomb interaction
        :math:`v(\mathbf{Q}; z, z') = (2\pi/Q)e^{-Q|z-z'|}` for :math:`|z-z'| < c/2`.

        Args:
            M (:py:class:`numpy.ndarray`) : matrix with shape (..., ng, ng) (the leading dimensions, e.g. the
                frequencies, are kept)
            iq (:py:class:`int`) : index of the q point (1-based), used only to check the arguments
            G_par (:py:class:`tuple`) : in-plane Miller indexes (m1, m2) of the block
            z (:py:class:`float` or :py:class:`numpy.ndarray`) : values of z (bohr)
            z2 (:py:class:`float` or :py:class:`numpy.ndarray`) : values of z' (bohr)

        Returns:
            :py:class:`numpy.ndarray` : array with shape (..., len(z), len(z2))
        """
        self._check_iq(iq)
        c = self._layer_c()
        M = np.asarray(M)
        if M.shape[-2:] != (self.ng, self.ng):
            raise ValueError('The matrix must have shape (..., %d, %d)' % (self.ng, self.ng))
        blocks = self.inplane_blocks()
        key = (int(G_par[0]), int(G_par[1]))
        if key not in blocks:
            raise ValueError('No G vector of the database has the in-plane Miller indexes %s' % (key,))
        idx = blocks[key]
        Gz = self.G[idx, 2]
        left = np.exp(1j * np.outer(np.atleast_1d(z), Gz))
        right = np.exp(-1j * np.outer(Gz, np.atleast_1d(z2)))
        return left @ M[..., idx, :][..., :, idx] @ right / c
