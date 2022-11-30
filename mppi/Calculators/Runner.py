class Runner():
    """
    This object is associated with the concept of execution of an action.
    It may be customized to be used inside workflows and datasets.
    The central functionality is the `run` method that can be customized on
    subclasses of `Runner`. In this object there are global and local options
    of a run method. All arguments passed at the instantiation are stored as
    global options. For each call to `run`, these global options may updated by
    the arguments of the run call.

    Args:
        **kwargs: global options of the runner. Deepcopied in the dictionary
          returned by :meth:`global_options`.

    Example:

        >>> torun=Runner(args1='one',args2='two')
        >>> print(torun.global_options())
        {'args1':'one','args2':'two'}
        >>> print(torun.get_global_option('args1'))
        'one'
    """

    def __init__(self, **kwargs):
        import copy
        self._global_options = copy.deepcopy(kwargs)

    def global_options(self):
        """
        Get all global options dict.

        Returns:
            :py:class:`dict`: The dictionary of the global options in its
            current status
        """
        return self._global_options

    def get_global_option(self, key):
        """

        Get one key in global options

        Args:
           key (string): the global option key

        Returns:
            The value of the global options labelled by ``key``

        """
        return self._global_options[key]

    def update_global_options(self, **kwargs):
        """
        Update the global options by providing keyword arguments.

        Args:
           **kwargs: arguments to be updated in the global options
        """
        self._global_options.update(kwargs)

    def pop_global_option(self, key):
        """
        Remove a given global option from the global option dictionary

        Args:
           key (string): the global option key

        Returns:
           The value of the global option
        """
        self._global_options.pop(key)

    def _run_options(self, **kwargs):
        """
        Create a local dictionary for a specific run.
        It combines the present status of global option with the local
        dictionary of the run and creates the run_options dictionary.
        This dictionary can be accessed during the definition of the
        process_run method. It contains all the relevant keys for the definition
        of the runner.
        """
        import copy
        self.run_options = copy.deepcopy(self._global_options)
        self.run_options.update(kwargs)

    def run(self, **kwargs):
        """
        Run method of the class. It performs the following actions:

         * Constructs the local dictionary to be passed as ``**kwargs`` to the
           `process_run` function;
         * Calls the :meth:`pre_processing` method (intended to prepare some
           actions associated to the :meth:`process_run` method);
         * Calls :meth:`process_run` with the dictionary returned by
           :meth:`pre_processing` as  `**kwargs`;
         * Update such dictionary with the results returned by
           :meth:`process_run` and call :meth:`post_processing`;
         * Returns the object passed by the call to :meth:`post_processing`
           class method

        Developers are therefore expected to override :meth:`pre_processing`
        :meth:`process_run` and :meth:`post_processing`,
        when subclassing :class:`Runner`.

        """
        from mppi.Utilities import Utils as f
        self._run_options(**kwargs)
        run_args = self.pre_processing()
        run_results = self.process_run(**run_args)
        f.dict_merge(dest=run_args, src=run_results)
        return self.post_processing(**run_args)

    def pre_processing(self):
        """
        Pre-treat the keyword arguments and the options, if needed.

        Returns:
           :py:class:`dict`: dictionary of the pre-treated keyword arguments
           that have to be actually considered by process_run.
        """
        return {}

    def process_run(self, **kwargs):
        """
        Main item of the runner, defines the information that have to be
        post_processed by post_processing.

        Args:
          **kwargs (:py:class:`dict`): keyword arguments as returned from the
            :meth:`pre_processing` method.

        Returns:
          :py:class:`dict`:
               dictionary objects to be passed to post_processing, once the
               dictionary returned by :meth:`pre_processing` has been updated
        """
        return kwargs

    def post_processing(self, **kwargs):
        """
        Post-processing, take the arguments as they are provided by the update
        of :meth:`process_run` and :meth:`pre_processing` methods.

        Returns:
           The final object that each call to the :meth:`run` method is
           supposed to provide.
        """
        return None
