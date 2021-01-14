"""
Module that manages the parsing of all the elements of a Yambo computation.
This parser aim to deal with the output files, the ``ns.db1`` database written in the
SAVE folder, the ``dipoles`` and all the databases created by yambo in the $jobname
folder.
"""

class YamboParser():
    """
    Class that perform the parsing starting from the results :py:class:`dict` built
    by the :class:`YamboCalculator` class. In the actual implementation of the class the
    parser is able to deal with the o- files, the dipoles database, the ``ndb.RT_G_PAR``
    and the ``ns.db1`` database written in the SAVE folder.

    Attributes:
        data (:class:`YamboOutputParser`) : contains the instance of the :class:`YamboOutputParser`
            class that manage the parsing of the ``o-* files``
        dipoles (:class:`YamboDipolesParser`) : contains the instance of the :class:`YamboDipolesParser`
            tclass hat manages the parsing of the ``dipoles`` database
        dft (:class:`YamboDftParser`) : contains the instance of the :class:`YamboDftParser` that
            manages the parsing of the ``ns.db1`` database
        RTGreen (:class:`YamboRTGlesserParser`) : contains the instance of the
            :class:`YamboRTGlesserParser` that manages the parsing of the ``ndb.RT_G_PAR`
            database

    """

    def __init__(self,results, verbose = False, extendOut = True):
        """
        Initialize the data member of the class.

        Args:
            results (:py:class:`dict`): The dictionary of the results built by the
                :class:`YamboCalculator` class
            verbose (:py:class:`boolean`) : Determine the amount of information provided on terminal
            extendOut (:py:class:`boolean`) : Determine which dictionary is used as reference for the
                            names of the variables in the :class:`YamboOutputParser`

        """
        from mppi.Parsers import YamboOutputParser
        from mppi.Parsers import YamboDftParser
        from mppi.Parsers import YamboDipolesParser
        from mppi.Parsers import YamboRTGlesserParser

        if 'output' in results:
            self.data = YamboOutputParser(results['output'],verbose=verbose,extendOut=extendOut)
        else:
            print('There are no o- files in the %s dictionary. Please check...'%results)
        for key,value in results.items():
            if key == 'dipoles':
                self.dipoles = YamboDipolesParser(value,verbose=verbose)
            if key == 'dft':
                self.dft = YamboDftParser(value,verbose=verbose)
            if key == 'RT_G_PAR':
                self.RTGreen = YamboRTGlesserParser(value,verbose=verbose)

    def get_info(self):
        """
        Provide information on the structure of the attributes of the class
        """
        if hasattr(self,'data'):
            self.data.get_info()
            print(' ')
        if hasattr(self,'dipoles'):
            self.dipoles.get_info()
            print(' ')
        if hasattr(self,'dft'):
            self.dft.get_info()
            print(' ')
        if hasattr(self,'RTGreen'):
            self.RTGreen.get_info()
