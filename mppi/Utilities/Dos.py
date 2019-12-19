"""
PyBigDFT
This module is used to display the Density of States (DoS).
"""

AU_eV = 27.21138386


class DiracSuperposition():
    """
    Defines as superposition of Dirac deltas which can be used to
    plot the density of states
    """

    def __init__(self, dos, wgts=[1.0]):
        """
        Parameters:
           dos: array containing the density of states per eack k-point.
                Should be of shape 2
           wgts: contains the weights of each of the k-points
        """
        import numpy as np
        self.dos = dos
        if isinstance(wgts, float):
            self.norm = [wgts]
        else:
            self.norm = wgts
        # set range for this distribution
        e_min = 1.e100
        e_max = -1.e100
        ddos = np.ravel(dos)
        if len(ddos) > 0:
            e_min = min(e_min, np.min(ddos) - 0.05 *
                        (np.max(ddos) - np.min(ddos)))
            e_max = max(e_max, np.max(ddos) + 0.05 *
                        (np.max(ddos) - np.min(ddos)))
        self.xlim = (e_min, e_max)

    def curve(self, xs, sigma, wgts=None):
        import numpy as np
        dos_g = 0.0
        idos = 0
        for norm, dos in zip(self.norm, self.dos):
            if wgts is not None:
                norms = wgts[idos]*norm
                idos += 1
            else:
                norms = np.ones(len(dos))*norm
            kptcurve = self.peaks(xs, dos, norms, sigma)
            dos_g += kptcurve
        return xs, dos_g

    def peak(self, omega, e, sigma):
        """
        Define if a peak is a Gaussian or a Lorenzian (temporarily only the
        gaussian is defined)
        """
        import numpy as np
        nfac = np.sqrt(2.0*np.pi)
        val = np.exp(- (omega - e)**2 / (2.0 * sigma**2))/(nfac*sigma)
        return val

    def peaks(self, xs, dos, norms, sigma):
        """
        Return the array of the whole set of peaks
        """
        curve = 0.0
        for e, nrm in zip(dos, norms):
            curve += self.peak(xs, e, sigma)*nrm
        return curve


# def _bandarray_to_data(jspin, bandarrays):
#     lbl = 'up' if jspin == 0 else 'dw'
#     kptlists = [[], []]
#     for orbs in bandarrays:
#         for ispin, norbs in enumerate(orbs.info):
#             if norbs == 0 or ispin != jspin:
#                 continue
#             # energy values
#             kptlists[0].append(orbs[ispin, :norbs])
#             # normalization
#             kptlists[1].append(orbs.kwgt*(1.0-2*ispin))
#             # print 'kpt',kptlists
#     return kptlists, lbl


class Dos():
    """
    Definition of the density of state class
    """

    def __init__(self, bandarrays=None, energies=None, evals=None, units='eV',
                 label='1', sigma=0.2/AU_eV, npts=2500, fermi_level=None,
                 norm=1.0, sdos=None):
        """
        Extract a quantity which is associated to the DoS, that can be plotted
        """
        import numpy as np
        self.ens = []
        self.labels = []
        # self.norms=[]
        self.npts = npts
        # if bandarrays is not None:
        #     self.append_from_bandarray(bandarrays, label)
        if evals is not None:
            self.append_from_dict(evals, label)
        if energies is not None:
            self.append(np.array([energies]), label=label, units=units,
                        norm=(np.array([norm])
                        if isinstance(norm, float) else norm))
        self.sigma = self.conversion_factor(units)*sigma
        self.fermi_level(fermi_level, units=units)
        if sdos is not None:
            self._embed_sdos(sdos)

    def _embed_sdos(self, sdos):
        self.sdos = []
        for i, xdos in enumerate(sdos):
            self.sdos.append({'coord': xdos['coord']})
            jdos = 0
            for subspin in xdos['dos']:
                if len(subspin[0]) == 0:
                    continue
                d = {'doslist': subspin}
                try:
                    self.ens[jdos]['sdos'].append(d)
                except KeyError:
                    self.ens[jdos]['sdos'] = [d]
                jdos += 1

    # def append_from_bandarray(self, bandarrays, label):
    #     """
    #     Add a new band array to the previous DoS. Important for kpoints DOS
    #     """
    #     import numpy as np
    #     for jspin in range(2):
    #         kptlists, lbl = _bandarray_to_data(jspin, bandarrays)
    #         self.append(np.array(kptlists[0]), label=label+lbl, units='AU',
    #                     norm=np.array(kptlists[1]))

    # def append_from_dict(self, evals, label):
    #     import numpy as np
    #     "Get the energies from the different flavours given by the dict"
    #     evs = [[], []]
    #     ef = None
    #     for ev in evals:
    #         occ = self.get_ev(ev, ['e_occ', 'e_occupied'])
    #         if occ:
    #             ef = max(occ)
    #         vrt = self.get_ev(ev, ['e_vrt', 'e_virt'])
    #         eigen = False
    #         if occ:
    #             eigen = occ
    #         if vrt:
    #             eigen = vrt
    #         if not eigen:
    #             eigen = self.get_ev(ev)
    #         if not occ and not vrt and eigen:
    #             ef = max(eigen)
    #         if not eigen:
    #             continue
    #         for i, e in enumerate(eigen):
    #             if e:
    #                 evs[i].append(e)
    #     for i, energs in enumerate(evs):
    #         if len(energs) == 0:
    #             continue
    #         self.append(np.array(energs), label=label,
    #                     units='AU', norm=1.0-2.0*i)
    #     if ef:
    #         self.fermi_level(ef, units='AU')

    def get_ev(self, ev, keys=None):
        "Get the correct list of the energies for this eigenvalue"
        res = False
        if keys is None:
            ener = ev.get('e')
            spin = ev.get('s')
            if ener and spin == 1:
                res = [ener]
            elif ener and spin == -1:
                res = [None, ener]
        else:
            for k in keys:
                if k in ev:
                    res = ev[k]
                    if not isinstance(res, list):
                        res = [res]
                    break
        return res

    def append(self, energies, label=None, units='eV', norm=1.0):
        if not isinstance(norm, float) and len(norm) == 0:
            return
        dos = self.conversion_factor(units)*energies
        self.ens.append({'dos': DiracSuperposition(dos, wgts=norm)})
        lbl = label if label is not None else str(len(self.labels)+1)
        self.labels.append(lbl)
        # self.norms.append(norm)
        self.range = self._set_range()

    def conversion_factor(self, units):
        if units == 'AU':
            fac = AU_eV
        elif units == 'eV':
            fac = 1.0
        else:
            raise ValueError('Unrecognized units ('+units+')')
        return fac

    def fermi_level(self, fermi_level, units='eV'):
        if fermi_level is not None:
            self.ef = fermi_level*self.conversion_factor(units)

    def _set_range(self, npts=None):
        import numpy as np
        if npts is None:
            npts = self.npts
        e_min = 1.e100
        e_max = -1.e100
        for dos in self.ens:
            mn, mx = dos['dos'].xlim
            e_min = min(e_min, mn)
            e_max = max(e_max, mx)
        return np.arange(e_min, e_max, (e_max-e_min)/npts)

    def curve(self, dos, norm, sigma=None):
        import numpy as np
        if sigma is None:
            sigma = self.sigma
        nrm = np.sqrt(2.0*np.pi)*sigma/norm
        dos_g = []
        for e_i in self.range:
            if len(dos.shape) == 2:
                nkpt = dos.shape[0]
                value = 0.0
                for ikpt in range(nkpt):
                    value += np.sum(np.exp(- (e_i - dos[ikpt, :])**2 /
                                           (2.0 * sigma**2))/nrm[ikpt])
            else:
                value = np.sum(
                    np.exp(- (e_i - dos[:])**2 / (2.0 * sigma**2))/nrm)
            # Append data corresponding to each energy grid
            dos_g.append(value)
        return np.array(dos_g)

    # def dump(self, sigma=None):
    #     "For Gnuplot"
    #     if sigma is None:
    #         sigma = self.sigma
    #     # data=[self.curve(dos,norm=self.norms[i],sigma=sigma)
    #     #       for i,dos in enumerate(self.ens)]
    #     data = [dos['dos'].curve(self.range, sigma=sigma)[1]
    #             for dos in self.ens]
    #
    #     for i, e in enumerate(self.range):
    #         print(e, ' '.join(map(str, [d[i] for d in data])))

    def plot(self, sigma=None, legend=True, xlmin=None, xlmax=None, ylmin=None,
             ylmax=None):
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Slider  # , Button, RadioButtons
        if sigma is None:
            sigma = self.sigma
        self.fig, self.ax1 = plt.subplots()
        self.plotl = []
        for i, dos in enumerate(self.ens):
            # self.plotl.append(self.ax1.plot(self.range,self.curve(dos,norm=self.norms[i],sigma=sigma),label=self.labels[i]))
            self.plotl.append(self.ax1.plot(
                *dos['dos'].curve(self.range, sigma=sigma),
                label=self.labels[i]))
        if xlmax is not None:
            plt.xlim(xmax=xlmax)
        if xlmin is not None:
            plt.xlim(xmin=xlmin)
        if ylmax is not None:
            plt.ylim(ymax=ylmax)
        if ylmin is not None:
            plt.ylim(ymin=ylmin)
        plt.xlabel('Energy [eV]', fontsize=18)
        plt.ylabel('DoS', fontsize=18)
        if self.ef is not None:
            plt.axvline(self.ef, color='k', linestyle='--')
            # self.ax1.annotate('Fermi level', xy=(self.ef,2),
            #        xytext=(self.ef, 10),
            #    arrowprops=dict(facecolor='white', shrink=0.05),
            # )
        if len(self.labels) > 1 and legend:
            plt.legend(loc='best')
        axcolor = 'lightgoldenrodyellow'
        # try:
        #     axsigma = plt.axes([0.2, 0.93, 0.65, 0.03], facecolor=axcolor)
        # except AttributeError:
        #     axsigma = plt.axes([0.2, 0.93, 0.65, 0.03], axisbg=axcolor)

        # self.ssig = Slider(axsigma, 'Smearing', 0.0, 0.4, valinit=sigma)
        # self.ssig.on_changed(self.update)
        # if hasattr(self, 'sdos') and self.sdos:
        #     self._set_sdos_selector()
        #     self._set_sdos()
        plt.show()

    # def _set_sdos_selector(self):
    #     import matplotlib.pyplot as plt
    #     from matplotlib.widgets import RadioButtons
    #     self.sdos_selector = RadioButtons(
    #         plt.axes([0.93, 0.05, 0.04, 0.11], axisbg='lightgoldenrodyellow'),
    #         ('x', 'y', 'z'), active=1)
    #     self.isdos = 1
    #     self.sdos_selector.on_clicked(self._update_sdos)

    # def _set_sdos(self):
    #     import numpy
    #     xs = self.sdos[self.isdos]['coord']
    #     self._set_sdos_sliders(numpy.min(xs), numpy.max(xs))
    #     self._update_sdos(0.0)  # fake value as it is unused

    # def _sdos_curve(self, sdos, vmin, vmax):
    #     import numpy
    #     xs = self.sdos[self.isdos]['coord']
    #     imin = numpy.argmin(numpy.abs(xs-vmin))
    #     imax = numpy.argmin(numpy.abs(xs-vmax))
    #     doslist = sdos[self.isdos]['doslist']
    #     # norms=self.sdos[self.isdos]['norms'][ispin]
    #     tocurve = [0.0 for i in doslist[imin]]
    #     for d in doslist[imin:imax+1]:
    #         tocurve = [t+dd for t, dd in zip(tocurve, d)]
    #     # tocurve=numpy.sum([ d[ispin] for d in doslist[imin:imax+1]],axis=0)
    #     return tocurve
    #     # float(len(xs))/float(imax+1-imin)*tocurve,norms

    # def _update_sdos(self, val):
    #     isdos = self.isdos
    #     if val == 'x':
    #         isdos = 0
    #     elif val == 'y':
    #         isdos = 1
    #     elif val == 'z':
    #         isdos = 2
    #     if isdos != self.isdos:
    #         self.isdos = isdos
    #         self._set_sdos()
    #
    #     vmin, vmax = (s.val for s in self.ssdos)
    #     if vmax < vmin:
    #         self.ssdos[1].set_val(vmin)
    #         vmax = vmin
    #     if vmin > vmax:
    #         self.ssdos[0].set_val(vmax)
    #         vmin = vmax
    #     # now plot the sdos curve associated to the given value
    #     sig = self.ssig.val
    #     curves = []
    #     for dos in self.ens:
    #         if 'sdos' not in dos:
    #             continue
    #         renorms = self._sdos_curve(dos['sdos'], vmin, vmax)
    #         curve = dos['dos'].curve(self.range, sigma=sig, wgts=renorms)
    #         curves.append(curve)
    #     if hasattr(self, '_sdos_plots'):
    #         for pl, curve in zip(self._sdos_plots, curves):
    #             pl[0].set_ydata(curve[1])
    #     else:
    #         self._sdos_plots = []
    #         for c in curves:
    #             self._sdos_plots.append(
    #                 self.ax1.plot(*c, label='sdos'))
    #     self.ax1.relim()
    #     self.ax1.autoscale_view()
    #     self.fig.canvas.draw_idle()

    # def _set_sdos_sliders(self, cmin, cmax):
    #     import matplotlib.pyplot as plt
    #     from futile.Figures import VertSlider
    #     if hasattr(self, 'ssdos'):
    #         self.ssdos[0].ax.clear()
    #         self.ssdos[0].__init__(
    #             self.ssdos[0].ax, 'SDos', cmin, cmax, valinit=cmin)
    #         self.ssdos[1].ax.clear()
    #         self.ssdos[1].__init__(self.ssdos[1].ax, '',
    #                                cmin, cmax, valinit=cmax)
    #     else:
    #         axcolor = 'red'
    #         axmin = plt.axes([0.93, 0.2, 0.02, 0.65], axisbg=axcolor)
    #         axmax = plt.axes([0.95, 0.2, 0.02, 0.65], axisbg=axcolor)
    #         self.ssdos = [
    #             VertSlider(axmin, 'SDos', cmin, cmax, valinit=cmin),
    #             VertSlider(axmax, '', cmin, cmax, valinit=cmax)]
    #     self.ssdos[0].valtext.set_ha('right')
    #     self.ssdos[1].valtext.set_ha('left')
    #     self.ssdos[0].on_changed(self._update_sdos)
    #     self.ssdos[1].on_changed(self._update_sdos)

    # def update(self, val):
    #     sig = self.ssig.val
    #     for i, dos in enumerate(self.ens):
    #         self.plotl[i][0].set_ydata(
    #             dos['dos'].curve(self.range, sigma=sig)[1])
    #         # self.plotl[i][0].set_ydata(self.curve(dos,norm=self.norms[i],sigma=sig))
    #
    #     self.ax1.relim()
    #     self.ax1.autoscale_view()
    #     self.fig.canvas.draw_idle()
