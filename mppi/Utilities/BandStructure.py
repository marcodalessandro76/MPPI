"""
This module contains some useful functions and a class for dealing with bands structures.
The module can be loaded in the notebook in one of the following way

>>> from mppi import Utilities as U

>>> U.BandStructure

or, for instance, to load only BandStructure

>>> from mppi.Utilities import BandStructure

>>> BandStructure

The path of a band structure computation with pw can be built with the build_pw_kpath function of
the mppi.Calculators.Tools module.

"""
import numpy as np

def parse_Ypp_output(data):
    """
    Extract the kpath, kpoints and bands from the dictionary with the columns of the o- file
    of a ypp band structure computation (as built by the YamboOutputParser). The columns of the file are
    the curvilinear abscissa along the path, the bands and the three components of the k points.

    Args:
        data (:py:class:`dict`) : dictionary with the columns col1, col2, ... of the o- file

    Returns:
        :py:class:`tuple` : the kpath, the kpoints (with shape (nk,3)) and the bands (with shape (nbands,nk))

    """
    ncols = len(data)
    kpath = np.array(data['col1'])
    kpoints = np.array([data['col%d'%ind] for ind in range(ncols-2,ncols+1)]).transpose()
    bands = np.array([data['col%d'%ind] for ind in range(2,ncols-2)])
    return kpath,kpoints,bands

class BandStructure():
    """
    Class to manage and plot the band structure.

    The class can be initialized in various way. Specific classmethods to perform the init using
    the output of both QuantumESPRESSO and Ypp are provided.

    Args:
        kpoints (:py:class:`array`) : array with the coordinates of the kpoints used to build the path.
            The coordinate used to express the kpoints are arbitrary, consistence with the ``high_sym_points``
            parameter is required
        bands (:py:class:`numpy.array`) : the element bands[i] contains the energies of the i-th
            band (in eV)
        kpath (:py:class:`array`) : array with the value of the curvilinear abscissa along the path.
            If this parameter is not provided the curvilinear abscissa is computed by the class assuming that
            kpoints are expressed in cartesian coordinates, so that their distance is built with the euclidean formula
        high_sym_points(:py:class:`dict`) : dictionary with the names and coordinates of the high_sym_points of the path
            If this parameter is not provided the high symmetry points are not marked when the plot method of the
            class is called

    """

    def __init__(self, kpoints, bands, kpath = None, high_sym_points = None):
        self.kpoints = np.array(kpoints)
        self.bands = np.array(bands)
        self.high_sym_points = high_sym_points
        if kpath is None : self.kpath = self.get_kpath()
        else : self.kpath = np.array(kpath)

    @classmethod
    def from_Pw(cls, results, high_sym_points = None, set_scissor = None, set_gap = None, set_direct_gap = None):
        """
        Initialize the BandStructure class from the result of a QuantumESPRESSO computation performed
        along a path. The class makes usage of the PwParser of this package. The kpoints are expressed in
        cartesian coordinates in units of 2pi/alat, so the high_sym_points have to be given in the same units.

        Args:
            results (:py:class:`string`) : the data-file-schema.xml that contains the result of the
                    QuantumESPRESSO computation
            high_sym_points(:py:class:`dict`) : dictionary with the names and the coordinates of the
                    high_sym_points of the path
            set_scissor (:py:class:`float`) : set the value of the scissor (in eV) that is added to the empty bands.
                    If a scissor is provided the set_gap and set_direct_gap parameters are ignored
            set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system. If set_gap is provided
                    the set_direct_gap parameter is ignored
            set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.

        """
        from mppi import Parsers as P
        data = P.PwParser(results,verbose=False)
        evals = data.get_evals(set_scissor=set_scissor,set_gap=set_gap,set_direct_gap=set_direct_gap)
        return cls(kpoints=data.kpoints,bands=evals.transpose(),high_sym_points=high_sym_points)

    @classmethod
    def from_Ypp(cls, results, high_sym_points = None, suffix = 'bands_interpolated'):
        """
        Initialize the BandStructure class from the results dictionary of a ypp band structure computation
        (as returned by the YamboCalculator). The curvilinear abscissa along the path is the one computed
        by ypp, and the kpoints are expressed in the coordinates of the cooOut variable of the ypp input, so the
        high_sym_points have to be given in the same coordinates.

        Args:
            results (:py:class:`dict`) : dictionary with the output of a Ypp computation
            high_sym_points(:py:class:`dict`) : dictionary with name and coordinates of the
                            high_sym_points of the path
            suffix (string) : specifies the suffix of the o- file use to build the bands

        """
        from mppi import Parsers as P
        data = P.YamboOutputParser(results['output'],verbose=False)
        kpath,kpoints,bands = parse_Ypp_output(data[suffix])
        return cls(kpath=kpath,kpoints=kpoints,bands=bands,high_sym_points=high_sym_points)

    @classmethod
    def from_Ypp_file(cls, file, high_sym_points = None):
        """
        Initialize the BandStructure class from the o- file written by a ypp band structure computation.

        Args:
            file (:py:class:`string`) : name of the o- file built by the Ypp postprocessing
            high_sym_points(:py:class:`dict`) : dictionary with name and coordinates of the
                            high_sym_points of the path

        """
        from mppi.Parsers import YamboOutputParser
        data = YamboOutputParser.from_file(file,verbose=False)
        kpath,kpoints,bands = parse_Ypp_output(list(data.values())[0])
        return cls(kpath=kpath,kpoints=kpoints,bands=bands,high_sym_points=high_sym_points)

    def get_kpath(self):
        """
        Compute the curvilinear abscissa along the path, assuming that the kpoints are expressed in
        cartesian coordinates.

        Returns:
            :py:class:`array` : values of the curvilinear abscissa along the path
        """
        steps = np.linalg.norm(np.diff(self.kpoints,axis=0),axis=1)
        return np.concatenate([[0.],np.cumsum(steps)])

    def get_high_sym_positions(self,atol=1e-4,rtol=1e-4):
        r"""
        Compute the position of the high_sym_points along the path. The method uses
        the numpy.allclose function to establish if the coordinates of a point on the
        path matches with an high_sym_points

        Args:
            atol (float) : absolute tolerance used by numpy.allclose
            rtol (float) : relative tolerance used by numpy.allclose

        Return:
            (tuple): tuple containing:
                (:py:class:`list`) : labels of the high symmetry points, in the order in which they are
                    found along the path. The label 'G' is converted to r'$\Gamma$' for a correct
                    rendering of the plot

                (:py:class:`list`) : coordinates on the path of the high symmetry points

        """
        if self.high_sym_points is None :
            return None

        found = []
        for point,coords in self.high_sym_points.items():
            for ind,k in enumerate(self.kpoints):
                if np.allclose(coords,k,rtol=rtol,atol=atol):
                    found.append((self.kpath[ind],r'$\Gamma$' if point == 'G' else point))
        found.sort()
        labels = [label for _,label in found]
        positions = [pos for pos,_ in found]
        return labels,positions

    def plot(self, plt, axes = None, selection = None, show_vertical_lines = True, **kwargs):
        """
        Plot the band structure.

        Args:
            plt (:py:class:`matplotlib.pyplot`) : the matplotlib object
            axes (:py:class:`matplotlib.pyplot.axes`) : the matplotlib axes object. If provided the plot
                is performed on the given axes
            selection (:py:class:`list`) : the list of bands that are plotted. If None all the
                bands are plotted. The band index starts from zero
            show_vertical_lines (:py:class:`bool`) : if True add the vertical lines with the positions
                of the high symmetry points on the path (if the high_sym_points variable is not None)
            kwargs : further parameter to edit the line style of the plot. The label (if given) is attributed
                only to the first plotted band, so that each band structure appears once in the legend

        """
        ax = axes if axes is not None else plt.gca()
        plotted_bands = range(len(self.bands)) if selection is None else selection

        label = kwargs.pop('label',None)
        for count,ind in enumerate(plotted_bands):
            ax.plot(self.kpath,self.bands[ind],label=label if count == 0 else None,**kwargs)

        high_sym_positions = self.get_high_sym_positions()
        if show_vertical_lines and high_sym_positions is not None :
            labels,positions = high_sym_positions
            for pos in positions:
                ax.axvline(pos,color='black',ls='--',lw=0.8)
            ax.set_xticks(positions)
            ax.set_xticklabels(labels,size=14)
