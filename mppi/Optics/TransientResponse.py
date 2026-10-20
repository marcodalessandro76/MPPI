r"""
This module defines the tools to compute the transient response of a pump-probe experiment, i.e. the change of the
response to a delta-shaped probe induced by the pump, from the equilibrium non-linear susceptibilities extracted from
the real-time polarization with :py:class:`Xn_single_frequency` and :py:class:`Xn_frequency_mixing`.
The module can be loaded in the notebook as follows

>>> from mppi.Optics import Transient_Response

>>> tr = Transient_Response.from_Xn(mix, single, broadening=0.3)
>>> w, dchi = tr.dchi(taus, E0P, pump='accumulated', sigma=3.2)

"""
import numpy as np
from scipy.special import erf
from mppi.Utilities import Constants as C
from mppi.Optics.Utils import time_weighted

PUMPS = ['monochromatic', 'probe time', 'accumulated']

def envelope_weight(tau, t, sigma, power=2, model='accumulated'):
    r"""
    Weight of the envelope :math:`f(t) = e^{-t^2/2\sigma^2}` of a Gaussian pump in the transient response at the time t
    after a delta-shaped probe at the delay :math:`\tau` from the maximum of the pump. With model = 'probe time' the weight
    is :math:`f^p(\tau)`, with model = 'accumulated' it is the mean of :math:`f^p` between the probe and t,

    .. math::
        W_p(\tau,t) = \frac{1}{t}\int_0^t f^p(\tau+t')\,dt'

    and :math:`f^p(\tau)` for :math:`t \leq 0`. The power p is 2 for the third order terms (intensity envelope) and 1 for the
    second order ones (field envelope, used by :py:class:`Transient_Response` only with the 'probe time' model).

    Args:
        tau (:py:class:`float`): delay of the probe from the maximum of the pump
        t (:py:class:`numpy.ndarray`): times after the probe, in the units of tau
        sigma (:py:class:`float`): width of the field envelope, in the units of tau
        power (:py:class:`int`): power p of the envelope. Default is 2
        model (:py:class:`str`): 'probe time' or 'accumulated'. Default is 'accumulated'

    Returns:
        :py:class:`numpy.ndarray`: the weight at the times t
    """
    t = np.asarray(t, dtype=float)
    fp = np.exp(-power * tau**2 / (2 * sigma**2))
    if model == 'probe time':
        return fp * np.ones_like(t)
    if model != 'accumulated':
        raise ValueError("The model of the envelope must be 'probe time' or 'accumulated'")
    c = np.sqrt(power / 2) / sigma
    tt = np.where(t > 0., t, 1.)
    W = np.sqrt(np.pi) / (2 * c * tt) * (erf(c * (tau + tt)) - erf(c * tau))
    return np.where(t > 0., W, fp)

class Transient_Response():
    r"""
    Class to compute the transient response of a pump-probe experiment from the equilibrium non-linear susceptibilities.

    The probe is a delta-shaped field at the time :math:`t_p` and the pump is :math:`E_{0P}f(t-t_P)\sin(\omega_P(t-t_P))`,
    with :math:`f = 1` for a monochromatic pump and :math:`f(t) = e^{-t^2/2\sigma^2}` for a Gaussian pulse. The transient
    response is :math:`\chi^{neq}(\omega,\tau) = (P_{Pp}(\omega)-P_P(\omega))/E_p(\omega) = \chi^{(1)}(\omega)+\delta\chi(\omega,\tau)`,
    where :math:`P_{Pp}` and :math:`P_P` are the polarizations with both the fields and with the pump only and
    :math:`\tau = t_p - t_P` is the delay of the probe. For a monochromatic pump, up to the third order in the fields,

    .. math::
        \delta\chi(\omega,\tau) = \frac{iE_{0P}}{2}\left[\chi^{(2)}_+(\omega)e^{-i\omega_P\tau} - \chi^{(2)}_-(\omega)e^{i\omega_P\tau}\right]
        + \frac{E_{0P}^2}{4}\left[\chi_0(\omega) - \chi_+(\omega)e^{-2i\omega_P\tau} - \chi_-(\omega)e^{2i\omega_P\tau}\right]

    where the minus signs come from the sine shape of the pump and the susceptibilities follow from the keys of
    :py:class:`Xn_frequency_mixing` (which include the degeneracy factors) and from the linear susceptibility
    :math:`\chi^{(1)}` of a run with the probe only:

    .. math::
        \chi_0(\omega) = \frac{\chi_{(1,0)}(\omega)-\chi^{(1)}(\omega)}{E_{0P,runs}^2/4} , \quad
        \chi_\pm(\omega) = \chi_{(1,\pm2)}(\omega\mp2\omega_P) , \quad
        \chi^{(2)}_\pm(\omega) = \chi_{(1,\pm1)}(\omega\mp\omega_P)

    The key (1,m) computed at the probe frequency :math:`\omega-m\omega_P` gives the output at :math:`\omega`, so the step of
    the grid of the probe frequencies must divide the pump frequency of the runs (or twice it, if only the third order terms
    are used) and :math:`\delta\chi` is defined on the frequencies of the grid at a distance of at least :math:`2\omega_P`
    (:math:`\omega_P` for the second order terms only) from its edges.

    For a Gaussian pump two models are available (see the method :py:meth:`dchi`):

    * 'probe time': the pump acts with its envelope at the probe time, the second order terms are multiplied by
      :math:`f(\tau)` and the third order ones by :math:`f^2(\tau)`. It requires that the response started by the probe decays
      before the envelope changes, i.e. a width of the resonances :math:`\eta \gg \hbar/\sigma`;
    * 'accumulated': the response functions of the third order susceptibilities are multiplied by the mean intensity
      envelope between the probe and the emission time (see :py:func:`envelope_weight`): the response at the time t after
      the probe feels the pump at all the times between the probe and t. For the term :math:`\chi_0` this follows from a
      pump that shifts the transition energies following its intensity envelope (pump far below the resonances), for the
      terms :math:`\chi_\pm` it is an extension of this model, not derived. On the anharmonic oscillator it reproduces the
      exact third order response within a few percent. The model is not defined for the second order terms, which are
      available only with the other two pumps.

    Corrections that depend on the code used to compute the susceptibilities (e.g. the effects of the time integrator)
    are not included.

    Args:
        freqs (:py:class:`numpy.ndarray`): uniform grid of the probe frequencies (eV)
        chi1 (:py:class:`numpy.ndarray`): linear susceptibility on the grid, from runs with the probe only
        keys (:py:class:`dict`): susceptibilities of :py:class:`Xn_frequency_mixing` (keys (1,m)) on the grid, from runs with
            probe and pump. The keys (1,0) and (1,+-2) are needed for the third order terms, (1,+-1) for the second order ones
        E0P (:py:class:`float`): amplitude of the pump of the frequency mixing runs (au)
        wP (:py:class:`float`): frequency of the pump of the frequency mixing runs (eV)
        centrosymmetric (:py:class:`bool`): if True the second order terms are not included, even if the keys (1,+-1) are
            provided (in a centrosymmetric system they vanish and their computed values are numerical noise). Default is False
        verbose (:py:class:`bool`): define the amount of information provided on terminal

    Attributes:
        freqs (:py:class:`numpy.ndarray`): grid of the probe frequencies (eV)
        step (:py:class:`float`): step of the grid (eV)
        wP (:py:class:`float`): frequency of the pump of the frequency mixing runs (eV)
        E0P (:py:class:`float`): amplitude of the pump of the frequency mixing runs (au)
        orders (:py:class:`list`): orders in the fields of the terms included (2 and/or 3)
        w (:py:class:`numpy.ndarray`): frequencies (eV) on which the transient response is defined
    """

    def __init__(self, freqs, chi1, keys, E0P, wP, centrosymmetric=False, verbose=True):
        self.freqs = np.asarray(freqs, dtype=float)
        if len(self.freqs) < 2:
            raise ValueError('The grid of the probe frequencies needs at least two frequencies')
        self.step = (self.freqs[-1] - self.freqs[0]) / (len(self.freqs) - 1)
        if not np.allclose(np.diff(self.freqs), self.step, rtol=1e-3, atol=0.):
            raise ValueError('The grid of the probe frequencies must be uniform and increasing')
        self.wP = wP
        self.E0P = E0P
        self.chi1 = np.asarray(chi1)
        self.keys = {key: np.asarray(val) for key, val in keys.items()}
        self.centrosymmetric = centrosymmetric
        self.verbose = verbose
        self.orders = []
        self._shifts = {} # shift of the grid (number of steps) of the keys (1,m)
        if not centrosymmetric and (1, 1) in keys and (1, -1) in keys:
            self.orders.append(2)
            for m in [1, -1]:
                self._shifts[m] = self._grid_shift(m)
        if all(key in keys for key in [(1, 0), (1, 2), (1, -2)]):
            self.orders.append(3)
            for m in [0, 2, -2]:
                self._shifts[m] = self._grid_shift(m)
        if not self.orders:
            raise ValueError('The keys (1,0) and (1,+-2) (third order) or (1,+-1) (second order, not centrosymmetric systems) '
                             'are needed')
        self._nmax = max(abs(n) for n in self._shifts.values())
        if 2 * self._nmax >= len(self.freqs):
            raise ValueError('The grid of the probe frequencies is shorter than %s times the pump frequency'
                             % (4 if 3 in self.orders else 2))
        self.w = self.freqs[self._nmax:len(self.freqs) - self._nmax]
        if verbose:
            self.get_info()

    @classmethod
    def from_Xn(cls, mix, single, broadening=None, renormalize=True, direction=0, centrosymmetric=False, verbose=True):
        """
        Initialize the class from the analysis of the runs with probe and pump (:py:class:`Xn_frequency_mixing`) and of
        the runs with the probe only (:py:class:`Xn_single_frequency`), computed on the same grid of probe frequencies.
        The susceptibilities are computed with the same a posteriori broadening.

        Args:
            mix (:py:class:`Xn_frequency_mixing`): analysis of the runs with probe and pump
            single (:py:class:`Xn_single_frequency`): analysis of the runs with the probe only
            broadening (:py:class:`float`): Lorentzian broadening (eV) added to the susceptibilities linear in the probe,
                see the compute_Xn methods. Default is None (no broadening)
            renormalize (:py:class:`bool`): normalization of the Lorentzian kernel, see the compute_Xn methods. Default is True
            direction (:py:class:`int`): cartesian direction of the polarization (that of the fields). Default is 0 (x)
            centrosymmetric (:py:class:`bool`): if True the second order terms are not included. Default is False
            verbose (:py:class:`bool`): define the amount of information provided on terminal
        """
        f_single, f_mix = single.fields_freqs, mix.probe_freqs
        if len(f_single) != len(f_mix) or not np.allclose(f_single, f_mix, rtol=1e-5, atol=0.):
            raise ValueError('The runs with the probe only and the ones with probe and pump must have the same probe frequencies')
        keys = mix.compute_Xn(broadening=broadening, renormalize=renormalize)[direction]
        chi1 = single.compute_Xn(broadening=broadening, renormalize=renormalize)[direction][1]
        keys = {key: val for key, val in keys.items() if key[0] == 1}
        return cls(mix.probe_freqs * C.HaToeV, chi1, keys, mix.pump['amplitude'], mix.pump_freq * C.HaToeV,
                   centrosymmetric=centrosymmetric, verbose=verbose)

    def _grid_shift(self, m):
        """
        Number of steps of the grid that corresponds to m times the pump frequency.
        """
        n = int(np.round(m * self.wP / self.step))
        if abs(n * self.step - m * self.wP) > 1e-3 * self.step:
            raise ValueError('The step of the grid of the probe frequencies (%s eV) must divide %s times the pump frequency '
                             '(%s eV)' % (self.step, abs(m), self.wP))
        return n

    def _on_output(self, y, m):
        """
        Values on the output grid w of the quantity y computed at the probe frequencies of the key (1,m).
        """
        n = self._shifts[m]
        return y[..., self._nmax - n:len(self.freqs) - self._nmax - n]

    def get_info(self):
        """
        Provide information on the grid and on the terms included
        """
        print('Probe frequencies: %d, from %.3f to %.3f eV, step %.4f eV' % (len(self.freqs), self.freqs[0],
                                                                            self.freqs[-1], self.step))
        print('Pump of the frequency mixing runs: %.4f eV, amplitude %.4e au' % (self.wP, self.E0P))
        print('Orders of the terms included:', self.orders, '(centrosymmetric system)' if self.centrosymmetric else '')
        print('Transient response defined from %.3f to %.3f eV (%d frequencies)' % (self.w[0], self.w[-1], len(self.w)))

    def components(self):
        r"""
        Susceptibilities of the transient response on the output grid w.

        Returns:
            :py:class:`dict`: with the keys 'chi_0', 'chi_+', 'chi_-' (third order) and 'chi2_+', 'chi2_-' (second order)
            of the terms included, see the documentation of the class
        """
        comp = {}
        if 3 in self.orders:
            comp['chi_0'] = self._on_output((self.keys[(1, 0)] - self.chi1) / (self.E0P**2 / 4), 0)
            comp['chi_+'] = self._on_output(self.keys[(1, 2)], 2)
            comp['chi_-'] = self._on_output(self.keys[(1, -2)], -2)
        if 2 in self.orders:
            comp['chi2_+'] = self._on_output(self.keys[(1, 1)], 1)
            comp['chi2_-'] = self._on_output(self.keys[(1, -1)], -1)
        return comp

    def dchi(self, taus, E0P, pump='monochromatic', sigma=None, wP=None, orders=None):
        r"""
        Compute the change :math:`\delta\chi(\omega,\tau)` of the response to a delta-shaped probe induced by the pump
        :math:`E_{0P}f(t-t_P)\sin(\omega_P(t-t_P))` (the change of the dielectric function is :math:`4\pi\delta\chi`).

        Args:
            taus (:py:class:`float` or :py:class:`numpy.ndarray`): delays :math:`\tau = t_p - t_P` of the probe (fs). For a
                Gaussian pump :math:`t_P` is the time of its maximum
            E0P (:py:class:`float`): amplitude of the pump (au)
            pump (:py:class:`str`): 'monochromatic' (f = 1), or a Gaussian pump with the envelope at the probe time
                ('probe time') or with the accumulated envelope ('accumulated', third order terms only: with the second
                order terms in the attribute orders use orders=[3]). Default is 'monochromatic'
            sigma (:py:class:`float`): width of the envelope of the Gaussian pump field, :math:`f(t) = e^{-t^2/2\sigma^2}` (fs)
            wP (:py:class:`float`): frequency of the pump (eV) used in the phases :math:`e^{\mp i\omega_P\tau}`. Default is
                None, that is the pump frequency of the frequency mixing runs (which is always used to read the keys
                on the shifted grids: a small difference, irrelevant in the frequency arguments of the susceptibilities,
                can matter in the phases)
            orders (:py:class:`list` or None): orders of the terms to include (a subset of the attribute orders). Default
                is None (all the orders of the attribute orders)

        Returns:
            :py:class:`tuple`: the frequencies w (eV) and the array dchi with shape (len(taus), len(w)), in au
        """
        if pump not in PUMPS:
            raise ValueError('The pump must be one of %s' % PUMPS)
        if pump != 'monochromatic' and sigma is None:
            raise ValueError('The width sigma of the envelope is needed for a Gaussian pump')
        orders = self.orders if orders is None else list(orders)
        if not set(orders) <= set(self.orders):
            raise ValueError('The orders %s are not available, the terms included are of order %s' % (orders, self.orders))
        if pump == 'accumulated' and 2 in orders:
            raise ValueError('The accumulated envelope is defined only for the third order terms: use orders=[3]')
        taus = np.atleast_1d(np.asarray(taus, dtype=float)) * C.FsToAu
        wph = (self.wP if wP is None else wP) / C.HaToeV
        freqs = self.freqs / C.HaToeV
        sig = None if sigma is None else sigma * C.FsToAu
        terms = []  # (power of the envelope, m, susceptibility at the probe frequencies, coefficient times the phase)
        if 2 in orders:
            terms += [(1, 1, self.keys[(1, 1)], lambda tau: 1j * E0P / 2 * np.exp(-1j * wph * tau)),
                      (1, -1, self.keys[(1, -1)], lambda tau: -1j * E0P / 2 * np.exp(1j * wph * tau))]
        if 3 in orders:
            chi0 = (self.keys[(1, 0)] - self.chi1) / (self.E0P**2 / 4)
            terms += [(2, 0, chi0, lambda tau: E0P**2 / 4),
                      (2, 2, self.keys[(1, 2)], lambda tau: -E0P**2 / 4 * np.exp(-2j * wph * tau)),
                      (2, -2, self.keys[(1, -2)], lambda tau: -E0P**2 / 4 * np.exp(2j * wph * tau))]
        dchi = np.zeros((len(taus), len(self.w)), dtype=complex)
        for itau, tau in enumerate(taus):
            for power, m, chi, coeff in terms:
                if pump == 'monochromatic':
                    y = chi
                elif pump == 'probe time':
                    y = chi * envelope_weight(tau, 0., sig, power, 'probe time')
                else:
                    # the weight in time is applied on the grid of the key and the result is shifted: a shift in
                    # frequency is a phase in time, which does not change the product W R
                    y = time_weighted(freqs, chi, lambda t: envelope_weight(tau, t, sig, power, 'accumulated'))
                dchi[itau] += coeff(tau) * self._on_output(y, m)
        return self.w, dchi
