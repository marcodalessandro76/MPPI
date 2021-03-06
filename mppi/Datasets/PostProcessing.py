"""
This module collects some useful post-processing functions that can be used in
the :class:`Dataset` class.
"""

def pw_parse_data(dataset):
    """
    Apply the PwParser to the elements of the results dictionary of the dataset.

    Args:
        dataset(:class:`Dataset`) : the instance of Dataset

    Returns:
        :py:class:`dict` : dictionary with the parsed data for all the (computed) runs
            of the dataset

    """
    from mppi import Parsers as P
    results = {}
    for run,data in dataset.results.items():
        results[run] = P.PwParser(data,verbose=False)
    return results

def pw_get_energy(dataset):
    """
    Extract the total energy from the results dictionary of the dataset.

    Args:
        dataset(:class:`Dataset`) : the instance of Dataset

    Returns:
        :py:class:`dict` : dictionary with the energy (in Hartree) for all the (computed) runs
            of the dataset
    """
    from mppi import Parsers as P
    energy = {}
    for run,data in dataset.results.items():
        results = P.PwParser(data,verbose=False)
        energy[run] = results.get_energy(convert_eV = False)
    return energy

def pw_get_gap(dataset):
    """
    Extract the value of the gap from the results dictionary of the dataset.

    Args:
        dataset(:class:`Dataset`) : the instance of Dataset

    Returns:
        :py:class:`dict` : dictionary with the gap (in eV) for all the (computed) runs
            of the dataset. Information on the nature of the gap (direct or indirect) are
            written on terminal
    """
    from mppi import Parsers as P
    gap = {}
    for run,data in dataset.results.items():
        results = P.PwParser(data,verbose=False)
        if results.get_gap() is not None:
            gap[run] = results.get_gap()['gap']
        else: gap[run] = 0
    return gap

def yambo_parse_data(dataset):
    """
    Apply the YamboParser to the elements of the results dictionary of the dataset.

    Args:
        dataset(:class:`Dataset`) : the instance of Dataset

    Returns:
        :py:class:`dict` : dictionary with the parsed data for all the (computed) runs
            of the dataset
    """
    from mppi import Parsers as P
    results = {}
    for run,data in dataset.results.items():
        results[run] = P.YamboParser(data,verbose=False)
    return results
