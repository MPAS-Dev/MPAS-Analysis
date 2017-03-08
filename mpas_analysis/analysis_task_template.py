'''
This is an example analysis task to be used as a template for new tasks.
It should be copied into one of the component folders (`ocean`, `sea_ice`,
`land_ice`, etc.) and modified as needed.

Don't forget to remove this docstring. (It's not needed.)

Authors
-------
Xylar Asay-Davis
'''

# import python modules here

# import mpas_analysis module here (those with relative paths starting with
# dots)
from ..shared.analysis_task import AnalysisTask


class MyTask(AnalysisTask):  # {{{
    '''
    <Describe the analysis task here.>

    Authors
    -------
    <List of authors>
    '''

    def __init__(self, config, myArg='myDefaultValue'):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        myNewArg : str, optional
            <Describe the arg>

        Authors
        -------
        <List of authors>

        '''
        # first, call the constructor from the base class (AnalysisTask)
        super(MyTask, self).__init__(config)

        # next, name the task, the component (ocean, sea_ice, etc.) and the
        # tags (if any) of the component ('timeSeries', 'climatology',
        # horizontalMap', 'index', 'transect', etc.)
        self.taskName = 'myTask'
        self.componentName = 'component'
        self.tags = ['tag1', 'tag2']

        # then, store any additional arguments for use later on.  These would
        # likely include things like the name of a field, region, month,
        # season, etc. to be analyzed so that the same subclass of AnalysisTask
        # can perform several different tasks (potentially in parallel)
        self.myArg = myArg
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError: if myArg has an invalid value

        Authors
        -------
        <List of authors>
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(MyTask, self).setup_and_check()

        # then, perform additional checks specific to this analysis
        possibleArgs = ['blah', 'thing', 'stuff']
        if self.myArg not in possibleArgs:
            raise ValueError('MyTask must be constructed with argument myArg '
                             'having one of the values\n'
                             '{}.'.format(possibleArgs))

        # }}}

    def run(self):  # {{{
        '''
        Runs the analysis task.

        Individual tasks (children classes of this base class) should first
        call this method to perform any common steps in an analysis task,
        then, perform the steps required to run the analysis task.

        Authors
        -------
        <List of authors>
        '''

        # here is where the main 'meat' of the analysis task goes

        self._my_sub_task('someText', arg2='differentText')
        return
        # }}}

    # here is where you add helper methods that are meant to be non-public
    # (they start with an underscore), meaning you don't expect anyone to
    # access them outside of this file.  Typically you won't put as much in
    # the docstring as you would for a public function or method.
    #
    # you can either pass arguments (with or without defaults) or you can
    # 'save' arguments as member variables of `self` and then get them back
    # (like `self.myArg` here).
    def _my_sub_task(self, arg1, arg2=None):  # {{{
        '''
        <Performs my favority subtask>
        '''

        # perform the task
        print 'myArg:', self.myArg
        print 'arg1:', arg1
        if arg2 is not None:
            print 'arg2:', arg2
        # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
