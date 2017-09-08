'''
This is an example analysis task to be used as a template for new tasks.
It should be copied into one of the component folders (`ocean`, `sea_ice`,
`land_ice`, etc.) and modified as needed.

Specific parts or variable names that need to be modified by the developer are
explicitly indicated. All the others are to be left unchanged.  Anything in
comments or docstrings in triangle brackets <> indicates that it needs to be
replaced with the appropriate variable name, text, etc. for the task.

Instructions for creating a new analysis task:
1. create a new task by copying this file to the appropriate folder (ocean,
   sea_ice, etc.) and modifying it as described below.  Take a look at
   mpas_analysis/shared/analysis_task.py for additional guidance.
2. note, no changes need to be made to mpas_analysis/shared/analysis_task.py
3. modify config.default (and possibly any machine-specific config files in
   configs/<machine>)
4. import new analysis task in mpas_analysis/<component>/__init__.py
5. add new analysis task to run_analysis.py under build_analysis_list:
      analyses.append(<component>.MyTask(config, myArg='argValue'))
   This will add a new object of the MyTask class to a list of analysis tasks
   created in build_analysis_list.  Later on in run_analysis, it will first
   go through the list to make sure each task needs to be generated
   (by calling check_generate, which is defined in AnalysisTask), then, will
   call setup_and_check on each task (to make sure the appropriate AM is on
   and files are present), and will finally call run on each task that is
   to be generated and is set up properly.

Don't forget to remove this docstring. (It's not needed.)

Authors
-------
Xylar Asay-Davis
'''

# import python modules here

# import mpas_analysis module here (those with relative paths starting with
# dots)
from ..shared.analysis_task import AnalysisTask


# Everywhere in this template, change MyTask to the name of your task, starting
# with a capital letter and using mixed case but no underscores.  Leave
# "AnalysisTask" as it is.  This indicates that MyTask has AnalysisTask as its
# parent, "inhereting" several variables and methods (functions) from it that
# can be helpful in initializing, checking whether to perform analysis and
# performing the analysis.  See AnalysisTask in shared/analysis_task.py
class MyTask(AnalysisTask):  # {{{
    '''
    <Briefly describe the analysis task here.  Just a 1-2 sentence description
    of what the task does should be sufficient.>

    Authors
    -------
    <List of authors>
    '''

    # This function is the "constructor", which is called when you want to
    # create a new object from your class.  Typically, it will just store a few
    # useful pieces of information like the name of the task, the component and
    # maybe some tags describing the task.
    # Keep self and config arguments.
    # self is a way to access the object that the class MyTask is describing,
    #    letting you access member variables and methods (member functions of a
    #    class), as you will see in the examples below.  All methods in a
    #    python class start with the argument self, which is not included in
    #    the list of arguments when you call a method of an object (because it
    #    is always included automatically).
    # config is an MpasAnalysisConfigParser object that can be used to get
    #    configuration options stored in config.default or a custom config
    #    file specific to a given simulation.  See examples below or in
    #    existing analysis tasks.
    # myArg should either be modified or removed. An example might be if you
    #    changed "myArg" to be "fieldName" (probably with no default value).
    #    Then, you would store it in self.fieldName so you could access it
    #    later in the script (in setup_and_check, run, or one of your helper
    #    methods).  In this example, you would have:
    #        def __init__(self, config, fieldName):...
    #    and yu would then make a new task something like this:
    #        myTask = MyTask(config, fieldName='seaIceArea')
    def __init__(self, config, myArg='myDefaultValue'):  # {{{
        '''
        Construct the analysis task.
        <Add any additional description of what happens during construction>

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        myArg : str, optional
            <Describe the arg, or just remove remove it>

        Authors
        -------
        <List of authors>

        '''
        # first, call the constructor from the base class (AnalysisTask).
        # Modify MyTask, "myTask", "component" and ['tag1', 'tag2'] below:
        # taskName is the same as the class name but with a lowercase letter
        # componentName is one of 'ocean' or 'seaIce' (same as the name of
        #     the folder where the task resides)
        # tags are some useful names describing the task ('timeSeries',
        #     'climatology', 'horizontalMap', 'index', 'transect', etc.) that
        #     can be used in the 'generate' config option, e.g.
        #     'all_climatology' of 'no_transect'.  Tasks that include computing
        # time series, indices and/or climatologies should include
        # 'timeSeries', 'index' and/or 'climatology' tags, as these are helpful
        # used to compute and update the start and end dates of each of these
        # analysis tasks as part of setting up the task.
        #
        # super(MyTask, self).<method>() is a way of calling
        # AnalysisTask.<method>, since AnalysisTask is the "super" or parent
        # class of MyTask.  In this case, we first call Analysis.__init__(...)
        # before doing our own initialization.
        super(MyTask, self).__init__(config=config,
                                     taskName='myTask',
                                     componentName='component',
                                     tags=['tag1', 'tag2'])

        # then, store any additional arguments for use in setup_and_check, run
        # or helper methods.  Extra arguments would likely include things like
        # the name of a field, region, month, season, etc. to be analyzed so
        # that the same subclass of AnalysisTask can perform several different
        # tasks (potentially in parallel)
        # Modify or remove the folowing:
        self.myArg = myArg
        # Example:
        # self.fieldName = fieldName

        # }}}

    # this function will be called to figure out if the analysis task should
    # run.  It should check if the input arguments to the task are supported,
    # if appropriate analysis member(s) (AMs) were turned on, if the necessary
    # output files exist, etc.  If there is a problem, an exception should be
    # raised (see the example below) so the task will not be run.
    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError: <if myArg has an invalid value; modify as needed>

        Authors
        -------
        <List of authors>
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        # Here, change MyTask to the name of your task but don't make any other
        # changes
        # Again, super(MyTask, self).<method>() is a way of calling
        # AnalysisTask.<method>, since AnalysisTask is the "super" or parent
        # class of MyTask.  In this case, we first call
        # Analysis.setup_and_check() before doing our own setup and checking.
        super(MyTask, self).setup_and_check()

        # then, perform additional checks specific to this analysis.
        # Change or remove these as needed.

        # This checks if myArg was one of the allowed values, and raises an
        # exception if not.
        # From our example above, if you had fieldName as myArg, you could
        # check to see if fieldName was one of the allowed fields, say 'sst'
        # or 'sss' if our task only supported sea surface temperature or
        # sea surface height as input fields.
        possibleArgs = ['blah', 'thing', 'stuff']
        if self.myArg not in possibleArgs:
            raise ValueError('MyTask must be constructed with argument myArg '
                             'having one of the values\n'
                             '{}.'.format(possibleArgs))

        # the following checks to see if the timeSeriesStatsMonthly AM was on
        # and raises an exception if not
        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        # the following checks to see if the MOC AM was on
        # and stores the result for later use (e.g. computing something
        # in MPAS-Analysis direclty if it was not computed in an AM)
        self.mocAnalysisMemberEnabled = self.check_analysis_enabled(
            analysisOptionName='config_am_mocstreamfunction_enable',
            raiseException=False)

        # The following is an example of how you could get a list of input
        # files from a specific stream (timeSeriesStats).  AM output should
        # be taken from the historyStreams, not runStreams, to support ACME
        # short-term archiving capabilities

        # Get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        streamName = 'timeSeriesStatsMonthlyOutput'
        startDate = self.config.get('climatology', 'startDate')
        endDate = self.config.get('climatology', 'endDate')
        self.inputFiles = self.historyStreams.readpath(streamName,
                                                       startDate=startDate,
                                                       endDate=endDate,
                                                       calendar=self.calendar)

        # Each analysis task generates one or more plots and writes out an
        # associated xml file for each plot.  Once all tasks have finished,
        # the "master" task will run through all the tasks and look at
        # xmlFileNames to find out what XML files were written out.  Each task
        # should provide a list of files in the order that the corresponding
        # images should appear on the webpage.

        # Note: because of the way parallel tasks are handled in MPAS-Analysis,
        # we can't be sure that run() will be called (it might be launched
        # as a completely separate process) so it is not safe to store a list
        # of xml files from within run(). The recommended procedure is to
        # create a list of XML files here during setup_and_check() and possibly
        # use them during run()

        self.xmlFileNames = []

        # }}}

    def run(self):  # {{{
        '''
        The main method of the task that performs the analysis task.

        Authors
        -------
        <List of authors>
        '''

        # Add the main contents of the analysis task below

        # No need to call AnalysisTask.run() because it doesn't do anything,
        # so we don't call super(MyTask, self).run(), as we do for other
        # methods above.

        # Here is an example of a call to a local helper method (function):
        self._my_sub_task('someText', arg2='differentText')
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
        <Performs my favorite subtask>
        '''

        # perform the task
        # self.myArg is a copy of the argument we passed in to __init__ when we
        # built the task.  It is available in any method after that for us to
        # use as needed.
        print 'myArg:', self.myArg
        print 'arg1:', arg1
        if arg2 is not None:
            print 'arg2:', arg2
        # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
