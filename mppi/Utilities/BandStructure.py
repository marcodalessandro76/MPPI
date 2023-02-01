"""
This module contains some useful functions and a class for dealing with bands structures.
The module can be loaded in the notebook in one of the following way

>>> from mppi import Utilities as U

>>> U.build_kpath

or, for instance, to load only BandStructure

>>> from mppi.Utilities import BandStructure

>>> BandStructure

"""

def _parse_Ypp_output(data):
    """
    Extract the kpath, kpoints and bands from the dictionary results['output'][type]
    given by YamboParser.
    """
    import numpy as np
    index_kx = len(data.keys())-3

    kpath = data['col0']
    kpoints = []
    for ind in [0,1,2]:
        kpoints.append(data['col'+str(index_kx+ind)])
    kpoints = np.array(kpoints).transpose()
    bands = []
    for ind in range(1,index_kx):
        bands.append(data['col'+str(ind)])
    bands = np.array(bands)

    return kpath,kpoints,bands

class BandStructure():
    """
    Class to manage and plot the band structure.

    The class can be initialized in various way. Specific classmethods to init using the output
    of both QuantumESPRESSO and Ypp are provided.

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
        self.kpoints = kpoints
        self.bands = bands
        self.high_sym_points = high_sym_points
        if kpath is None : self.kpath = self.get_kpath()
        else : self.kpath = kpath

    @classmethod
    def from_Pw(cls, results, high_sym_points = None, set_scissor = None, set_gap = None, set_direct_gap = None):
        """
        Initialize the BandStructure class from the result of a QuantumESPRESSO computation performed
        along a path. The class makes usage of the PwParser of this package.

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
        kpoints = data.kpoints
        evals = data.get_evals(set_scissor=set_scissor,set_gap=set_gap,set_direct_gap=set_direct_gap)
        bands = evals.transpose()
        return cls(kpoints=kpoints,bands=bands,high_sym_points=high_sym_points)

    @classmethod
    def from_Ypp(cls, results, high_sym_points = None, suffix = 'bands_interpolated'):
        """
        Initialize the BandStructure class from the data dictionay of the result of a Ypp postprocessing.
        The class make usage of the YamboParser of this package.

        Args:
            results (:py:class:`dict`) : dictionary with the output of a Ypp computation
            high_sym_points(:py:class:`dict`) : dictionary with name and coordinates of the
                            high_sym_points of the path
            suffix (string) : specifies the suffix of the o- file use to build the bands

        """
        from mppi import Parsers as P
        import numpy as np
        data = P.YamboParser(results).data
        kpath,kpoints,bands = _parse_Ypp_output(data[suffix])
        return cls(kpath=None,kpoints=kpoints,bands=bands,high_sym_points=high_sym_points)

    @classmethod
    def from_Ypp_file(cls, results, high_sym_points = None):
        """
        Initialize the BandStructure class from the a o-file built by the Ypp postprocessing.
        The class make usage of the YamboParser of this package.

        Args:
            results (:py:class:`string`) : name of the o-file built by the Ypp postprocessing
            high_sym_points(:py:class:`dict`) : dictionary with name and coordinates of the
                            high_sym_points of the path
            suffix (string) : specifies the suffix of the o- file use to build the bands

        """
        from mppi.Utilities.Utils import file_parser
        import numpy as np
        data = file_parser(results)
        kpoints = data[-3:,:].T
        bands = data[1:-3]
        return cls(bands=bands,kpoints=kpoints,high_sym_points=high_sym_points)

    def get_kpath(self):
        """
        Compute the curvilinear ascissa along the path.

        Returns:
            :py:class:`array` : values of the curvilinear ascissa along the path
        """
        import numpy as np
        kpoints = self.kpoints
        kpath = [0]
        distance = 0
        for nk in range(1,len(kpoints)):
            distance += np.linalg.norm(kpoints[nk]-kpoints[nk-1])
            kpath.append(distance)
        return np.array(kpath)

    def get_high_sym_positions(self,atol=1e-4,rtol=1e-4):
        """
        Compute the position of the high_sym_points along the path. The method uses
        the numpy.allclose function to establish if the coordinates of a point on the
        path matche with an high_sym_points

        Args:
            atol (float) : absolute tolerance used by numpy.allclose
            rtol (float) : relative tolerance used by numpy.allclose

        Return:
            (tuple): tuple containing:
                (:py:class:`list`) : labels of the high symmetry points as found along the path.
                    If the point 'G' is found its label is converted to r'$\Gamma' for a correct
                    rendering of the plot

                (:py:class:`tuple`) : coordinate on the path of the high symmetry points

        """
        import numpy as np

        if self.high_sym_points is None :
            return None

        high_sym = self.high_sym_points
        kpoints = self.kpoints
        kpath = self.kpath

        labels = []
        positions = []
        for point in high_sym:
            for ind,k in enumerate(kpoints):
                if np.allclose(high_sym[point],k,atol,rtol):
                    if point == 'G': labels.append(r'$\Gamma$')
                    else : labels.append(point)
                    positions.append(kpath[ind])
        return labels,positions

    def plot(self, plt, axes = None, selection = None, show_vertical_lines = True, **kwargs):
        """
        Plot the band structure.

        Args:
            plt (:py:class:`matplotlib.pyplot`) : the matplotlib object
            axes (:py:class:`matplotlib.pyplot.axes`) : the matplotlib axes object. If provided the plot
                is performed on the given axes
            selection (:py:class:`list`) : the list of bands that are plotted. If None all the
                bands computed by QuantumESPRESSO are plotted. The band index starts
                from zero
            show_vertical_lines (:py:class:`bool`) : if True add the vertical lines with the positions
                of the high symmetry points on the path (if the high_sym_points variable is not None)
            kwargs : further parameter to edit the line style of the plot

        """
        kpath = self.kpath
        nbands = len(self.bands)
        high_sym_positions = self.get_high_sym_positions()
        if high_sym_positions is not None:
            labels,positions = high_sym_positions

        plotted_bands = list(ind for ind in range(nbands)) if selection is None else selection

        if axes is not None:
            ax = axes
        else:
            ax = plt.gca()

        for ind in plotted_bands:
            ax.plot(kpath,self.bands[ind],**kwargs)
        if show_vertical_lines and high_sym_positions is not None :
            for pos in positions:
                ax.axvline(pos,color='black',ls='--')
                ax.set_xticks(positions)
                ax.set_xticklabels(labels,size=14)
