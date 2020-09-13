"""
This module collects some useful postprocessing functions that can be used in
the Dataset class.
"""
def QE_parse_data(dataset):
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
        results[run] = P.PwParser(data['output'],verbose=False)
    return results

def QE_get_energy(dataset):
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
        results = P.PwParser(data['output'],verbose=False)
        energy[run] = results.get_energy(convert_eV = False)
    return energy

def Yambo_parse_data(dataset):
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
        results[run] = P.YamboParser(data['output'],verbose=False)
    return results
