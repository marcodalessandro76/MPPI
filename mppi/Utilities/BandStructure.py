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
    Extract the kpoints and bands from the dictionary results['output'][type]
    given by YamboParser.
    """
    import numpy as np
    index_kx = len(data.keys())-3

    kpoints = []
    for ind in [0,1,2]:
        kpoints.append(data['col'+str(index_kx+ind)])
     # kpoints is converted to np.array to perfrom the transpose and convert back
     # to list since it is needed as a list by BandStructure
    kpoints = list(map(list,list(np.array(kpoints).transpose())))
    bands = []
    for ind in range(1,index_kx):
        bands.append(data['col'+str(ind)])
    bands = np.array(bands)

    return kpoints, bands

class BandStructure():
    """
    Class to manage and plot the band structure.

    The class can be initialized in various way. Specific classmethods to init using the output
    of both QuantumESPRESSO and Ypp are provided.

    Args:
        kpoints (:py:class:`list`) : list with the coordinates of the kpoints
            used to build the path
        bands (:py:class:`numpy.array`) : the element bands[i] contains the energies of the i-th
            band (in eV)
        high_sym_points(:py:class:`dict`) : dictionary with the names and coordinates of the
                        high_sym_points of the path

    """

    def __init__(self, bands, kpoints, high_sym_points):
        self.bands = bands
        self.kpoints = kpoints
        self.high_sym_points = high_sym_points

    @classmethod
    def from_Pw(cls, results, high_sym_points, set_scissor = None, set_gap = None, set_direct_gap = None):
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
        return cls(bands=bands,kpoints=kpoints,high_sym_points=high_sym_points)

    @classmethod
    def from_Ypp(cls,results,high_sym_points,suffix='bands_interpolated'):
        """
        Initialize the BandStructure class from the result of a Ypp postprocessing.
        The class make usage of the YamboParser of this package.

        Args:
            results (:py:class:`list`) : list that contiains the o- file provided as the output of a
                Ypp computation. results is the key ['output'][irun] of the run method of YamboCalculator
            high_sym_points(:py:class:`dict`) : dictionary with name and coordinates of the
                            high_sym_points of the path
            suffix (string) : specifies the suffix of the o- file use to build the bands

        """
        from mppi import Parsers as P
        import numpy as np
        data = P.YamboParser(results)
        kpoints, bands = _parse_Ypp_output(data[suffix])
        return cls(bands=bands,kpoints=kpoints,high_sym_points=high_sym_points)

    def get_bands(self):
        """
        Returns the band attribute of the class
        """
        return self.bands

    def get_path(self):
        """
        Compute the curvilinear ascissa along the path.

        Returns:
            kpath(array) : values of the curvilinear ascissa along the path
        """
        import numpy as np
        kpoints = np.array(self.kpoints)
        kpath = [0]
        distance = 0
        for nk in range(1,len(kpoints)):
            distance += np.linalg.norm(kpoints[nk-1]-kpoints[nk])
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
                (:py:class:`list`) : labels of the high symmetry points as found along the path

                (:py:class:`list`) : coordinate on the path of the high symmetry points

        """
        import numpy as np

        kpoints = self.kpoints
        high_sym = self.high_sym_points
        path = self.get_path()
        labels = []
        positions = []

        for point in high_sym:
            for ind,k in enumerate(kpoints):
                if np.allclose(high_sym[point],k,atol,rtol):
                    labels.append(point)
                    positions.append(path[ind])
        return labels,positions

    def plot(self, plt, selection = None, show_vertical_lines = True, **kwargs):
        """
        Plot the band structure.

        Args:
            plt(:py:class:`matplotlib.pyplot`) : the matplotlib object
            selection (list) : the list of bands that are plotted. If None all the
                bands computed by QuantumESPRESSO are plotted. The band index starts
                from zero
            kwargs : further parameter to edit the line style of the plot

        """
        path = self.get_path()
        labels,positions = self.get_high_sym_positions()
        nbands = len(self.bands)

        plotted_bands = list(ind for ind in range(nbands)) if selection is None else selection

        for ind in plotted_bands:
            plt.plot(path,self.bands[ind],**kwargs)
        if show_vertical_lines :
            for pos in positions:
                plt.axvline(pos,color='black',ls='--')

        ax = plt.gca()
        ax.set_xticks(positions)
        ax.set_xticklabels(labels,size=14)
