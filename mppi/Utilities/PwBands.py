"""
This file contains some useful functions fro managing and plotting bands structures
starting from a pw nscf calculation.
In this type of computation the complete list of kpoints along a given path is passed
to QuantumESPRESSO.

The module can be loaded in the notebook in one of the following way

>>> from mppi import Utilities as U

>>> U.build_kpath

or, for instance, to load only the PwBands

>>> from mppi.Utilities import PwBands

>>> PwBands

"""

def build_kpath(*kpoints,numstep=40):
    """
    Build a list of kpoints to be passed to the set_kpoints methods of the
    :class:`PwInput` for computing the band structure along a path.

    Args:
        klist(*list): specifies the high symmetry points along the k-path
        numstep(int): specifies the number of intermediate points used to build
                      the path

    Return:
        (list) : list of kpoints as nedded by pw in the nscf computation with
                tpiba_b option

    """
    klist = []
    for k in kpoints[:-1]:
        klist.append(k+[numstep])
    klist.append(kpoints[-1]+[0])
    return klist

class PwBands():
    """
    Class to manage and plot the band structure starting from the result of a nscf
    QuantumESPRESSO computation performed along a path. The class make usage of
    the PwParser of this package.

    Args:
        results (string) : the data-file-schema.xml that contains the result of the
                        nscf computation
        high_sym_points(:py:class:`dict`) : dictionary with name and coordinates of the
                        high_sym_points of the path
    """

    def __init__(self,results,high_sym_points):
        from mppi import Parsers as P
        self.data = P.PwParser(results,verbose=False)
        self.high_sym_points = high_sym_points

    def get_bands(self,convert_eV=True):
        """
        Convert the array data.evals into the the array bands, where bands[i] gives
        the energies of the i-th band along the path.
        The fermi level is used as the reference energy and the results are
        expressed in eV if convert_eV = True.
        """
        from mppi.Utilities import HaToeV
        import numpy as np
        evals = self.data.evals
        bands = []
        for b in range(len(evals[0])): #number of bands
            bands.append(evals[:,b])
        bands = np.array(bands) - self.data.fermi
        if convert_eV: bands *= HaToeV
        return bands

    def get_path(self):
        """
        Compute the curvilinear ascissa along the path.

        Returns:
        kpath(array) : values of the curvilinear ascissa along the path
        """
        import numpy as np
        kpoints = np.array(self.data.kpoints)
        kpath = [0]
        distance = 0
        for nk in range(1,len(kpoints)):
            distance += np.linalg.norm(kpoints[nk-1]-kpoints[nk])
            kpath.append(distance)
        return np.array(kpath)

    def get_high_sym_positions(self):
        """
        Compute the position of the high_sym_points along the path. Return two lists,
        the first one with the labels and the second one with the positions of the
        high symmetry points.

        """
        kpoints = self.data.kpoints
        high_sym = self.high_sym_points
        path = self.get_path()
        labels = []
        positions = []

        for point in high_sym:
            for ind,k in enumerate(kpoints):
                if high_sym[point] == k:
                    labels.append(point)
                    positions.append(path[ind])
        return labels,positions

    def plot(self,plt,convert_eV=True,selection = None):
        """
        Plot the band structure.

        Args:
            plt(:py:class:`matplotlib.pyplot`) : the matplotlib object
            selection (list) : the list of bands that are plotted. If none all the
                bands computed by QuantumESPRESSO are plotted. The band index starts
                from zero
        """

        bands = self.get_bands(convert_eV)
        path = self.get_path()
        labels,positions = self.get_high_sym_positions()

        plotted_bands = list(ind for ind in range(self.data.nbands)) if selection is None else selection

        for ind in plotted_bands:
            plt.plot(path,bands[ind])
        for pos in positions:
            plt.axvline(pos,color='black',ls='--')

        ax = plt.gca()
        ax.set_xticklabels(labels,size=14)
        ax.set_xticks(positions)
        plt.show()

        # Display the band gap
        self.data.get_gap()
