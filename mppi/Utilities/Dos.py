"""
This module defines the tools to build and manage the Density of States (Dos).
A Dos can be built starting from various inputs like the output of QuantumESPRESSO or a
Yambo computations or from a generic one-dimensional array. Different Dos can be managed together.
"""

def lorentzian(eta,x0,x):
    """
    Get the lorentzian function.
    """
    import numpy as np
    s2 = eta**2
    c = eta/np.pi
    x1 = (x-x0)**2
    return c/(x1+s2)

def gaussian(eta,x0,x):
    """
    Get the gaussian function.
    """
    import numpy as np
    c = 1.0/(eta*np.sqrt(2.0*np.pi))
    x1 = ((x-x0)/eta)**2
    return c*np.exp(-0.5*x1)

def build_histogram(values,weights=None,norm=1.0,minVal=None,maxVal=None,
                    step=0.01,eta=0.05,broad_kind=lorentzian):
    """
    This function builds the histogram associated to a generic one-dimensional array.
    If the weights are not specified a uniform array of weights, normalized to norm, is assumed.

    Args:
        values (:class:`numpy.array`) : one-dimensional array with the values used to build the Dos
        weights (:py:class:`numpy.array`) : one-dimensional array with the weight of each value.
            If None a uniform array normalized to one is used
        norm (:py:class:`float`) : the normalization of the (uniform) weights array
        minVal (:py:class:`float`) : values lower than this parameter are not included in the histogram
        maxVal (:py:class:`float`) : values higher than this parameter are not included in the histogram
        step (:py:class:`float`) : size of the bin (in the same units used for the values array)
        eta (:py:class:`float`) : magnitude of the broading parameter (in the same units used for the values array)
        broad_kind (:py:class:`string`) : type of broading function used (lorentzian, gaussian)
        label (:py:class:`string`) : label associated to the dos

    Return:
        (tuple) : tuple containing:
            (:py:class:`numpy.array`) : x axis of the histogram
            (:py:class:`numpy.array`) : histogram values

    """
    import numpy as np
    # set the weights and the range
    if weights is None:
        weights = np.ones(len(values))*norm/len(values)
    if minVal is None:
        minVal = values.min() - 100.*eta
    if maxVal is None:
        maxVal = values.max() + 100.*eta

    weights = weights[minVal < values]
    values = values[minVal < values]
    weights = weights[values < maxVal]
    values = values[values < maxVal]

    x = np.arange(minVal,maxVal,step)
    histo = np.zeros([len(x)])
    for v,w in zip(values,weights):
        histo += w*broad_kind(eta,v,x)

    return (x, histo)

def convert_PwData(evals,weights):
    """
    Convert the arrays with the structure of evals and weights of the PwParser class
    into the form suitable to be managed by the build_histogram function.

    Args:
        evals (:py:class:`numpy.array`) : array with the structure of the self.evals of
            PwParser
        weights (:py:class:`numpy.array`) : array with the structure of the self.weights of
            PwParser

    Return:
        (tuple) : tuple containing:
            (:py:class:`numpy.array`) : one-dimensioanal array with the energies
            (:py:class:`numpy.array`) : one-dimensional array with the associated weights

    """
    import numpy as np
    weights = np.ones(evals.shape)*weights
    return (evals.flatten(),weights.flatten())

class Dos():
    """
    Definition of the density of state class.
    The dos is normalized so that its integral is equal to the norm of the weights divided the step of
    x axis sampling.

    Attributes:
        dos (:py:class:`list`): list with the tuple (energies,histrogram) for each dos appended to the class
        labels (:py:class:`list`): list with the labels of the appended dos

    Args:
        energies (:class:`numpy.array`) : one-dimensional array with the values used to build the Dos
        weights (:py:class:`numpy.array`) : one-dimensional array with the weight of each value.
            If None a uniform array normalized to one is used
        minVal (float) : values lower than this parameter are not included in the histogram
        maxVal (float) : values higher than this parameter are not included in the histogram
        step (float) : size of the bin (in the same units used for the values array)
        eta (float) : magnitude of the broading parameter (in the same units used for the values array)
        broad_kind (string) : type of broading function used (lorentzian, gaussian)

    """

    def __init__(self, energies = None, weights = None, norm = 1.0, minVal = None, maxVal = None,
                 step = 0.01, eta = 0.05, broad_kind = lorentzian, label = None):
        self.dos = []
        self.labels = []
        if energies is not None:
            self.append(energies,weights=weights,norm=norm,minVal=minVal,maxVal=maxVal,
                        step=step,eta=eta,broad_kind=broad_kind,label=label)
    @classmethod
    def from_Pw(cls,results,set_gap=None,set_direct_gap=None, minVal = None, maxVal = None,
                step = 0.01, eta = 0.05, broad_kind = lorentzian,label=None):
        """
        Initialize the Dos class from the xml output file of a QuantumESPRESSO computation.
        The class makes usage of the PwParser of this package.

        Args:
            results (:py:class:`string`) : the data-file-schema.xml that contains the result of the
                            QuantumESPRESSO computation
            set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system
            set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.
                            If set_gap is provided this parameter is ignored
            minVal (float) : values lower than this parameter are not included in the histogram
            maxVal (float) : values higher than this parameter are not included in the histogram
            step (float) : size of the bin (in the same units used for the values array)
            eta (float) : magnitude of the broading parameter (in the same units used for the values array)
            broad_kind (string) : type of broading function used (lorentzian, gaussian)

        """
        from mppi import Parsers as P
        data = P.PwParser(results,verbose=False)
        evals = data.get_evals(set_gap,set_direct_gap)
        energies, weights = convert_PwData(evals,data.weights)
        return cls(energies,weights=weights,label=label,minVal=minVal,maxVal=maxVal,step =step,eta=eta,broad_kind=broad_kind)

    def append(self,energies, weights = None, norm = 1.0, minVal = None, maxVal = None,
               step = 0.01, eta = 0.05, broad_kind = lorentzian, label = None):
        """
        This method add the tuple (x,histo) generated by the function build_histogram
        to the dos members of the class. The label of the new dos is added to the labels
        member.

        Args:
            energies (:class:`numpy.array`) : one-dimensional array with the values used to build the Dos
            weights (:py:class:`numpy.array`) : one-dimensional array with the weight of each value.
                If None a uniform array normalized to one is used
            norm (:py:class:`float`) : the normalization of the (uniform) weights array
            minVal (:py:class:`float`) : values lower than this parameter are not included in the histogram
            maxVal (:py:class:`float`) : values higher than this parameter are not included in the histogram
            step (:py:class:`float`) : size of the bin (in the same units used for the values array)
            eta (:py:class:`float`) : magnitude of the broading parameter (in the same units used for the values array)
            broad_kind (:py:class:`string`) : type of broading function used (lorentzian, gaussian)
            label (:py:class:`string`) : label associated to the dos

        """
        self.dos.append(build_histogram(energies,weights=weights,norm=norm,minVal=minVal,maxVal=maxVal,
                   step=step,eta=eta,broad_kind=broad_kind))
        lbl = label if label is not None else str(len(self.labels)+1)
        self.labels.append(lbl)

    def append_fromPw(self,results,set_gap=None,set_direct_gap=None, minVal = None, maxVal = None,
                      step = 0.01, eta = 0.05, broad_kind = lorentzian,label=None):
        """
        Add one element to the Dos class starting from the xml output file of a QuantumESPRESSO
        computation.

        Args:
            results (:py:class:`string`) : the data-file-schema.xml that contains the result of the
                            QuantumESPRESSO computation
            set_gap (:py:class:`float`) : set the value of the gap (in eV) of the system
            set_direct_gap (:py:class:`float`) : set the value of the direct gap (in eV) of the system.
                            If set_gap is provided this parameter is ignored
            label (string) : the label of the appended dos
            minVal (float) : values lower than this parameter are not included in the histogram
            maxVal (float) : values higher than this parameter are not included in the histogram
            step (float) : size of the bin (in the same units used for the values array)
            eta (float) : magnitude of the broading parameter (in the same units used for the values array)
            broad_kind (string) : type of broading function used (lorentzian, gaussian)

        """
        from mppi import Parsers as P
        data = P.PwParser(results,verbose=False)
        evals = data.get_evals(set_gap,set_direct_gap)
        energies, weights = convert_PwData(evals,data.weights)
        self.append(energies,weights=weights,label=label,minVal=minVal,maxVal=maxVal,step =step,eta=eta,broad_kind=broad_kind)

    def append_fromPwData(self,evals,weights, minVal = None, maxVal = None,
                          step = 0.01, eta = 0.05, broad_kind = lorentzian, label = None):
        """
        Add one element to the Dos class starting from arrays with the structure of the
        evals and weights attributes of the PwParser class. This method can be used to
        build a JDos, using the transitions as evals.

        Args:
            evals (:py:class:`numpy.array`) : array with the structure of the self.evals of
                the PwParser
            weights (:py:class:`numpy.array`) : array with the structure of the self.weights of
                the PwParser
            label (string) : the label of the appended dos
            minVal (float) : values lower than this parameter are not included in the histogram
            maxVal (float) : values higher than this parameter are not included in the histogram
            step (float) : size of the bin (in the same units used for the values array)
            eta (float) : magnitude of the broading parameter (in the same units used for the values array)
            broad_kind (string) : type of broading function used (lorentzian, gaussian)

        """
        energies, weights = convert_PwData(evals,weights)
        self.append(energies,weights=weights,label=label,minVal=minVal,maxVal=maxVal,step =step,eta=eta,broad_kind=broad_kind)

    def plot(self, plt, rescale = False, include = None, legend = True):
        """
        Plot the elements of the Dos class

        Args:
            rescale (:py:class:`bool`) : if True all the dos are rescaled to the same maximum value equal to 1.0 (useful for comparison)
            include (:py:class:`list`) : list with the indexes (as appended to the dos member ) of the dos that are plotted
            legend (:py:class:`bool`) : if True show the legend of the plot

        """
        to_plot = include if include is not None else range(len(self.dos))
        for ind in to_plot:
            scale = 1.0
            if rescale: scale = max(self.dos[ind][1])
            plt.plot(self.dos[ind][0],self.dos[ind][1]/scale,label=self.labels[ind])
        if legend:
            plt.legend()
