r"""
Class to read the (possibly cut) Coulomb interaction used by a Yambo run, from the ``ndb.cutoff`` and ``ndb.RIM``
databases. Both depend only on the geometry of the run (cell, q grid, Coulomb cutoff, random points of the RIM), so
they can be shared by several runs, e.g. by a screening database and by a modified copy of it.

- ``ndb.cutoff`` contains the ``bare_qpg`` of Yambo, :math:`b_\mathbf{G}(\mathbf{q})`, with
  :math:`v_\mathbf{G}(\mathbf{q}) = 4\pi/b_\mathbf{G}(\mathbf{q})^2` (:math:`b = |\mathbf{q}+\mathbf{G}|` without
  cutoff), for the first G vectors of the global list of the ``SAVE`` (the database does not store the G vectors).
  These are the values used in the screening: the matrices X of the screening databases contain
  :math:`v^{1/2}_\mathbf{G}v^{1/2}_{\mathbf{G}'}`.
- ``ndb.RIM`` contains the Monte Carlo average of :math:`1/b_\mathbf{G}(\mathbf{q}+\mathbf{q}')^2` over the small cell
  of the BZ around each q point (random integration method), used by Yambo for the exchange self-energy, the head of
  the BSE kernel and to build ``ndb.RIM_W``, and needed at q = 0 where v diverges.

The two files may be in different folders and a run may have one of them, both or none. The module can be loaded in
the notebook as follows

>>> from mppi.Parsers import YamboCoulombParser

>>> coul = YamboCoulombParser('run_dir/SAVE', cutoff='run_dir/ndb_rim_cutoff', rim='run_dir/ndb_rim_cutoff')
>>> v = coul.coulomb(2, ng=100)

"""
import os
import numpy as np
from mppi.Parsers.YamboScreeningParser import _read_text, _complex

Q0_DEF_NORM = 1e-5      # norm of the small vector that replaces q = 0 in Yambo (q0_def_norm)

def _database_file(path, name):
    """
    The database ``name`` given as a file or as the folder that contains it.
    """
    filename = os.path.join(path, name) if os.path.isdir(path) else path
    if not os.path.isfile(filename):
        raise IOError('File %s not found' % filename)
    return filename

class YamboCoulombParser():
    r"""
    Class to read the Coulomb interaction of a Yambo run from the ``ndb.cutoff`` and/or ``ndb.RIM`` databases. The q
    points are numbered from 1, as in Yambo. All the quantities are in atomic units (bohr, Ha).

    Args:
        save (:py:class:`string`) : the ``SAVE`` folder of the computation (or its ``ns.db1`` file), for the lattice,
            the lattice parameters of Yambo and the global list of the G vectors
        cutoff (:py:class:`string`) : the ``ndb.cutoff`` database or the folder that contains it. Default None
        rim (:py:class:`string`) : the ``ndb.RIM`` database or the folder that contains it. Default None. The two paths
            are independent; at least one of them is required
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        cutoff, rim : the files read (None if not given)
        geometry (:py:class:`string`) : the Coulomb cutoff of ndb.cutoff ('slab z', ...), 'none' if only ndb.RIM is
            given (the CUTOFF written in ndb.RIM is not reliable and is not used)
        lattice (:py:class:`numpy.ndarray`) : lattice vectors as rows (bohr)
        alat (:py:class:`numpy.ndarray`) : the lattice parameters of Yambo (alat_x, alat_y, alat_z), bohr
        volume (:py:class:`float`) : volume of the cell (bohr^3)
        G (:py:class:`numpy.ndarray`) : the global list of the G vectors of the SAVE (cartesian, bohr^-1)
        G_crystal (:py:class:`numpy.ndarray`) : their Miller indexes
        q (:py:class:`numpy.ndarray`) : the q points of the IBZ (cartesian, bohr^-1), shape (nq, 3)
        q_crystal (:py:class:`numpy.ndarray`) : the q points in crystal coordinates
        nq, nqbz (:py:class:`int`) : number of q points in the IBZ and in the BZ
        ng_cut (:py:class:`int`) : number of G vectors of ndb.cutoff (None without it)
        rim_ng (:py:class:`int`) : number of G vectors of ndb.RIM (None without it)
        rim_is_diagonal (:py:class:`bool`) : True if ndb.RIM contains only the diagonal elements G = G'
        rim_points (:py:class:`int`) : number of random points of the RIM
        rim_RL_volume, RL_volume (:py:class:`float`) : volume of the BZ estimated by the Monte Carlo integration and
            exact one (bohr^-3)
    """

    def __init__(self, save, cutoff=None, rim=None, verbose=True):
        from netCDF4 import Dataset
        from mppi.Parsers.YamboDftParser import YamboDftParser
        if cutoff is None and rim is None:
            raise ValueError('Give the path of ndb.cutoff, of ndb.RIM or of both')
        self.verbose = verbose
        db1 = _database_file(save, 'ns.db1')
        dft = YamboDftParser(db1, verbose=False)
        self.lattice = np.array(dft.lattice, dtype=float)
        self.alat = np.array(dft.alat, dtype=float)
        self.volume = abs(np.linalg.det(self.lattice))
        self.reciprocal_lattice = 2 * np.pi * np.linalg.inv(self.lattice).T
        with Dataset(db1) as d:
            self.G = np.array(d['G-VECTORS'][:], dtype=float).T * 2 * np.pi / self.alat
        G_crys = self.G @ self.lattice.T / (2 * np.pi)
        self.G_crystal = np.rint(G_crys).astype(int)

        heads = {}
        self.cutoff = self.rim = None
        self.geometry = 'none'
        self.ng_cut = self.rim_ng = None
        self._bare_qpg = self._rim_qpg = None
        if cutoff is not None:
            self.cutoff = _database_file(cutoff, 'ndb.cutoff')
            with Dataset(self.cutoff) as d:
                self.geometry = _read_text(d, 'DESCRIPTION') if 'DESCRIPTION' in d.variables else \
                    _read_text(d, 'GEOMETRY')
                heads['ndb.cutoff'] = (np.array(d['HEAD_QPT'][:], dtype=float).T,
                                       np.array(d['HEAD_R_LATT'][:], dtype=float))
                b = _complex(d['CUT_BARE_QPG'][:]).T       # netCDF (ng_cut, nq, 2) -> (nq, ng_cut)
            self._bare_qpg = b.real if np.abs(b.imag).max() == 0 else b
            self.ng_cut = b.shape[1]
            if self.ng_cut > len(self.G):
                raise ValueError('ndb.cutoff has %d G vectors, the SAVE %d' % (self.ng_cut, len(self.G)))
        if rim is not None:
            self.rim = _database_file(rim, 'ndb.RIM')
            with Dataset(self.rim) as d:
                pars = np.array(d['PARS'][:], dtype=float)
                heads['ndb.RIM'] = (np.array(d['HEAD_QPT'][:], dtype=float).T,
                                    np.array(d['HEAD_R_LATT'][:], dtype=float))
                r = np.array(d['RIM_qpg'][:], dtype=float)
            self._rim_qpg = np.transpose(r, (2, 1, 0))     # Fortran RIM_qpg(iq, G, G') -> [iq, G, G']
            self.rim_ng = int(pars[0])
            self.rim_is_diagonal = bool(int(pars[1]))
            self.rim_points = int(pars[2])
            self.rim_RL_volume = pars[3]
            self.RL_volume = pars[4]
        if len(heads) == 2:
            (q1, r1), (q2, r2) = heads.values()
            if q1.shape != q2.shape or not np.allclose(q1, q2, atol=1e-5) or not np.array_equal(r1[:2], r2[:2]):
                raise ValueError('The q points of %s and %s do not agree: they are not from the same geometry'
                                 % (self.cutoff, self.rim))
        q_iku, rlatt = list(heads.values())[0]
        self.q = q_iku * 2 * np.pi / self.alat
        self.nq = len(self.q)
        self.nqbz = int(rlatt[1])
        self.q_crystal = self.q @ self.lattice.T / (2 * np.pi)
        if self.rim is None:
            self.RL_volume = (2 * np.pi)**3 / self.volume
        elif self.cutoff is None and verbose:
            print('Warning: without ndb.cutoff the point values of v (coulomb, bare_qpg) are computed without cutoff')
        if verbose:
            self.get_info()

    def get_info(self):
        """
        Provide information on the data read
        """
        print('Coulomb interaction: cutoff "%s"' % self.geometry)
        print('number of q points', self.nq, '(IBZ),', self.nqbz, '(BZ)')
        if self.cutoff is not None:
            print('ndb.cutoff: %s, %d G vectors' % (self.cutoff, self.ng_cut))
        if self.rim is not None:
            print('ndb.RIM: %s, %d G vectors%s, %d random points' % (self.rim, self.rim_ng,
                  ' (diagonal)' if self.rim_is_diagonal else '', self.rim_points))

    def _check_iq(self, iq):
        if not (1 <= int(iq) <= self.nq):
            raise ValueError('The index of the q point must be between 1 and %d' % self.nq)
        return int(iq)

    def bare_qpg(self, iq, ng=None, q0_norm=Q0_DEF_NORM):
        r"""
        The ``bare_qpg`` of Yambo, :math:`b_\mathbf{G}(\mathbf{q})` with
        :math:`v_\mathbf{G}(\mathbf{q}) = 4\pi/b_\mathbf{G}(\mathbf{q})^2`, on the first ng G vectors of the SAVE.
        Read from ndb.cutoff if it is given, otherwise (no cutoff) :math:`|\mathbf{q}+\mathbf{G}|`, with q = 0 at the
        first q point and q0_norm for G = 0, as in Yambo. Without ndb.cutoff the values are correct only if the run did
        not use a Coulomb cutoff (ndb.RIM does not tell it reliably).

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)
            ng (:py:class:`int`) : number of G vectors. Default None (all the G vectors of ndb.cutoff, or of the SAVE
                without it)
            q0_norm (:py:class:`float`) : the norm of q used for G = 0 at q = 0 without ndb.cutoff (bohr^-1). Default
                1e-5, the value of Yambo

        Returns:
            :py:class:`numpy.ndarray` : array with shape (ng,) (bohr^-1)
        """
        iq = self._check_iq(iq)
        if self._bare_qpg is not None:
            ng = self.ng_cut if ng is None else int(ng)
            if ng > self.ng_cut:
                raise ValueError('ndb.cutoff contains only %d G vectors' % self.ng_cut)
            return self._bare_qpg[iq - 1, :ng]
        ng = len(self.G) if ng is None else int(ng)
        b = np.linalg.norm(self.q[iq - 1] + self.G[:ng], axis=1)
        if iq == 1:
            b[np.all(self.G_crystal[:ng] == 0, axis=1)] = q0_norm
        return b

    def coulomb(self, iq, ng=None, q0_norm=Q0_DEF_NORM):
        r"""
        Coulomb interaction :math:`v_\mathbf{G}(\mathbf{q}) = 4\pi/b_\mathbf{G}(\mathbf{q})^2` (Ha bohr^3) on the first
        ng G vectors of the SAVE (see :py:meth:`bare_qpg`). These are the values contained in the matrices X of the
        screening databases.

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)
            ng (:py:class:`int`) : number of G vectors. Default None (see :py:meth:`bare_qpg`)
            q0_norm (:py:class:`float`) : the norm of q used for G = 0 at q = 0 without ndb.cutoff (bohr^-1)

        Returns:
            :py:class:`numpy.ndarray` : array with shape (ng,)
        """
        v = 4 * np.pi / self.bare_qpg(iq, ng, q0_norm)**2
        return v.real if np.iscomplexobj(v) and np.abs(v.imag).max() == 0 else v

    def coulomb_rim(self, iq):
        r"""
        Coulomb interaction averaged over the small cell of the BZ around q, from ndb.RIM (Ha bohr^3):
        :math:`\langle v\rangle_\mathbf{G}(\mathbf{q}) = 4\pi\,R_{\mathbf{G}\mathbf{G}}(\mathbf{q})/(2 w_q)`, where R is
        the ``RIM_qpg`` of Yambo and :math:`w_q = (\Omega_{BZ}/N_q)/(2\pi)^3` (as in ``scatter_Gamp.F``). At q = 0 and
        G = 0 the average includes the integrable divergence of v.

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)

        Returns:
            :py:class:`numpy.ndarray` : array with shape (rim_ng,) if ndb.RIM is diagonal, (rim_ng, rim_ng) otherwise
        """
        if self._rim_qpg is None:
            raise ValueError('ndb.RIM was not given to the constructor')
        iq = self._check_iq(iq)
        q_weight = self.RL_volume / self.nqbz / (2 * np.pi)**3
        r = self._rim_qpg[iq - 1]
        r = np.diag(r) if self.rim_is_diagonal else r
        return 4 * np.pi * r / (2 * q_weight)

    def slab_bare_qpg(self, iq, ng=None, rim=None):
        r"""
        Analytic ``bare_qpg`` of the slab cutoff of Yambo (``cutoff_slab.F``), on the first ng G vectors of the SAVE:

        .. math::
            b_\mathbf{G}(\mathbf{q}) = |\mathbf{q}+\mathbf{G}| / \sqrt{1 - e^{-Q_\parallel l}\cos((q+G)_\perp l)}

        with :math:`l = a_\perp/2` (alat of the normal axis), :math:`Q_\parallel \geq 10^{-5}`; at the first q point
        q = 0 and :math:`|\mathbf{q}+\mathbf{G}| = 10^{-5}` for G = 0, where without RIM Yambo uses an average over a
        circle of the area of the small cell of the BZ.

        Args:
            iq (:py:class:`int`) : index of the q point (1-based)
            ng (:py:class:`int`) : number of G vectors. Default None (ng_cut, or all the G vectors of the SAVE)
            rim (:py:class:`bool`) : True if the run used the RIM (then the value at q = 0, G = 0 is not replaced).
                Default None (True if ndb.RIM was given)

        Returns:
            :py:class:`numpy.ndarray` : array with shape (ng,) (bohr^-1)
        """
        words = self.geometry.lower().split()
        if len(words) < 2 or words[0] != 'slab' or words[1] not in ('x', 'y', 'z'):
            raise ValueError('The cutoff of the run is "%s", not a slab' % self.geometry)
        iq = self._check_iq(iq)
        ng = (self.ng_cut or len(self.G)) if ng is None else int(ng)
        rim = (self.rim is not None) if rim is None else rim
        a = 'xyz'.index(words[1])
        others = [i for i in range(3) if i != a]
        lcut = self.alat[a] / 2
        qG = (self.q[iq - 1] if iq > 1 else np.zeros(3)) + self.G[:ng]
        norm = np.linalg.norm(qG, axis=1)
        g0 = np.all(self.G_crystal[:ng] == 0, axis=1)
        if iq == 1:
            norm[g0] = Q0_DEF_NORM
        Qpar = np.maximum(np.linalg.norm(qG[:, others], axis=1), Q0_DEF_NORM)
        b = norm / np.sqrt(1 - np.exp(-Qpar * lcut) * np.cos(qG[:, a] * lcut))
        if iq == 1 and not rim:
            sup = self.RL_volume / self.nqbz / (2 * np.pi / self.alat[a])
            radius = np.sqrt(sup / np.pi)
            series, factorial = 0., 1.
            for n in range(1, 16):
                factorial *= n
                series += (-1)**(n + 1) / n / factorial * (radius * lcut)**n
            b[g0] = 1 / np.sqrt(2 * np.pi / sup * series)
        return b

    def check_slab(self, rim=None):
        """
        Compare the ``bare_qpg`` of ndb.cutoff with the analytic slab formula (:py:meth:`slab_bare_qpg`). This also
        checks the correspondence between the values of ndb.cutoff and the G list of the SAVE, which the database does
        not store.

        Args:
            rim (:py:class:`bool`) : True if the run used the RIM (see :py:meth:`slab_bare_qpg`). Default None (True if
                ndb.RIM was given)

        Returns:
            :py:class:`dict` : maximum relative deviations, with the keys 'q' (q points 2, ..., nq), 'q=0' (first q
            point, G != 0) and 'q=0,G=0'
        """
        if self._bare_qpg is None:
            raise ValueError('ndb.cutoff was not given to the constructor')
        dev = {'q': 0., 'q=0': 0., 'q=0,G=0': 0.}
        g0 = np.all(self.G_crystal[:self.ng_cut] == 0, axis=1)
        for iq in range(1, self.nq + 1):
            b, ref = self.bare_qpg(iq), self.slab_bare_qpg(iq, rim=rim)
            rel = np.abs(b - ref) / np.abs(ref)
            if iq == 1:
                dev['q=0'] = rel[~g0].max()
                dev['q=0,G=0'] = rel[g0].max()
            else:
                dev['q'] = max(dev['q'], rel.max())
        if self.verbose:
            print('ndb.cutoff vs analytic slab formula, max relative deviation: %.1e (q != 0), %.1e (q = 0, G != 0), '
                  '%.1e (q = 0, G = 0)' % (dev['q'], dev['q=0'], dev['q=0,G=0']))
        return dev
