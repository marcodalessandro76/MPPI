"""
Module that manages the parsing of all the elements of a Yambo computation.
This parser aim to deal with the output files, the ``ns.db1`` database written in the
SAVE folder, the ``dipoles`` and all the databases created by yambo in the $jobname
folder.
"""
from mppi import Parsers as P

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

        if 'output' in results:
            self.data = P.YamboOutputParser(results['output'],verbose=verbose,extendOut=extendOut)
        else:
            print('There are no o- files in the %s dictionary. Please check...'%results)
        for key,value in results.items():
            if key == 'dipoles':
                self.dipoles = P.YamboDipolesParser(value,verbose=verbose)
            if key == 'dft':
                self.dft = P.YamboDftParser(value,verbose=verbose)
            if key == 'RT_G_PAR':
                self.RTGreen = P.YamboRTGlesserParser(value,verbose=verbose)

    @classmethod
    def from_path(cls, run_dir, outputPath, dbsPath = None , verbose = True, extendOut = True):
        """
        Init the a :class:`YamboParser` instance using the 'o-' files found inside the
        outputPath, the ``ns.db1`` database in the SAVE folder and the databases found in the dbsPath.

        Args:
            run_dir (:py:class:`string`) : `run_dir` folder of the calculation
            outputPath (:py:class:`string`) : folder with the 'o-' files
            dbsPath (:py:class:`string`) : folder with the ndb databases. If it is
                None the databases are sought in the outputPath
            verbose (:py:class:`boolean`) : determine the amount of information provided on terminal
            extendOut (:py:class:`boolean`) : Determine which dictionary is used as reference for the
                            names of the variables in the :class:`YamboOutputParser`

        """
        from mppi.Calculators.YamboCalculator import build_results_dict
        results = build_results_dict(run_dir,outputPath,dbsPath=dbsPath,verbose=verbose)
        return cls(results,verbose=verbose)

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
