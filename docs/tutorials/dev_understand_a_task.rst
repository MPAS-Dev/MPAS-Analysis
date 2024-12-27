.. _tutorial_understand_a_task:

Developers: Understanding an analysis task
==========================================

This tutorial walks a new developer through an existing analysis task to get
a more in-depth understanding of the code.  This tutorial is meant as the
starting point for the :ref:`tutorial_dev_add_task` tutorial.  It is a common
practice to find an existing analysis task that is as close as possible to the
new analysis, and to copy that existing task as a template for the new task.
This tutorial describes an existing analysis task, and
:ref:`tutorial_dev_add_task` uses it as a starting point for developing a new
task.

1. The big picture
------------------

MPAS-Analysis is meant to provide a first look at E3SM simulation output from
the MPAS-Ocean and MPAS-Seaice components.  The analysis is intended to be
robust and automated.  However, there is currently little effort to ensure that
the time period covered by the observations and model output are the same.
In other words, we often compare pre-industrial simulation results with
present-day observations.  The justification for this is twofold.  First, we
typically have few if any observations covering the pre-industrial period.
Second, we may be attempting to reduce biases that we assess to be much larger
than expected differences between pre-industrial and present-day climate
conditions.  Under these conditions, MPAS-Analysis provides us with a useful
first impression of how our simulation is doing.

1.1 MPAS output
~~~~~~~~~~~~~~~

The primary output from MPAS-Ocean and MPAS-Seaice are monthly and daily
averages of a large number of data fields.  Here are links to the list of:

* `MPAS-Ocean monthly fields <https://github.com/E3SM-Project/E3SM/blob/main/components/mpas-ocean/cime_config/buildnml#L1082-L1290>`_

* `MPAS-Ocean daily fields <https://github.com/E3SM-Project/E3SM/blob/main/components/mpas-ocean/cime_config/buildnml#L1028-L1062>`_

* `MPAS-Seaice monthly fields <https://github.com/E3SM-Project/E3SM/blob/main/components/mpas-seaice/cime_config/buildnml#L726-L823>`_

* `MPAS-Seaice daily fields <https://github.com/E3SM-Project/E3SM/blob/main/components/mpas-seaice/cime_config/buildnml#L695-L700>`_

The components also produce a smaller amount of more specialized output, such
as monthly maximum/minimum values.

MPAS data is provided on unstructured meshes, meaning that it isn't
particularly amenable to analysis with standard tools such as
`ESMValTool <https://www.esmvaltool.org/>`_.  Additionally, E3SM's science
campaigns require unique, sometimes regionally focused analysis not available
in existing tools.

1.2 Analysis tasks
~~~~~~~~~~~~~~~~~~

MPAS-Analysis is designed to a series of interdependent analysis tasks in
parallel with one another.  It builds up dependency graph between the tasks,
allowing independent tasks to run at the same time while putting dependent
tasks on hold until the tasks they depend on are completed.  Additionally,
MPAS-Analysis has some rudimentary tools for keeping track of the resources
that some computationally intensive tasks require to prevent the tool from
running out of memory.

Currently, nearly all operations in MPAS-Analysis must run on a single HPC
node.  (The exception is
`ncclimo <http://nco.sourceforge.net/nco.html#ncclimo-netCDF-Climatology-Generator>`_,
which is used to generate climatologies, and which can run in parallel across
up to 12 nodes if desired.) We hope to support broader task parallelism in
the not-too-distant future using the `parsl <https://parsl-project.org/>`_
python package.

Each analysis tasks is a class that descends from the
:py:class:`~mpas_analysis.shared.AnalysisTask` base class.  Tasks
can also have "subtasks" that do part of the work needed for the final
analysis.  A subtask might perform a computation on a specific region, period
of time, or season.  It might combine data from other subtasks into a single
dataset.  Or it might plot the data computed by a previous task.  The
advantages of dividing up the work are 1) that each subtask can potentially
run in parallel with other subtasks and 2) it can allow code reuse if the same
subtask can be used across multiple analysis tasks.


1.3 Shared framework
~~~~~~~~~~~~~~~~~~~~

MPAS-Analysis includes a shared framework used across analysis tasks.  The
framework is made up mostly of functions that can be called from within
analysis tasks but also includes some analysis tasks and subtasks that are
common to MPAS-Ocean, MPAS-Seaice and potentially other MPAS components
(notably MALI) that may be supported in the future.

This tutorial will not go though the shared framework in detail.  In addition
to the :py:class:`~mpas_analysis.shared.AnalysisTask` base class, the shared
framework includes the following
packages:

.. code-block:: none

    $ ls mpas_analysis/shared
    climatology
    constants
    generalized_reader
    html
    interpolation
    io
    mpas_xarray
    plot
    projection
    regions
    time_series
    timekeeping
    transects
    ...

A separate tutorial will explore the shared framework and how how to modify it.

2. Tour of an analysis task (``ClimatologyMapOHCAnomaly``)
----------------------------------------------------------

Aside from some code that takes care of managing analysis tasks and generating
web pages, MPAS-Analysis is made up almost entirely of analysis tasks and
shared functions they can call.  Since adding new analysis nearly always
mean creating a new class for the task, we start with a tour of an existing
analysis task as well as the :py:class:`~mpas_analysis.shared.AnalysisTask`
base class that it descends from.

We will use :py:class:`~mpas_analysis.ocean.ClimatologyMapOHCAnomaly` as an
example analysis task for this tour because it will turn out to be a useful
staring point for the analysis we want to add in :ref:`tutorial_dev_add_task`.
You can read more about :ref:`task_climatologyMapOHCAnomaly` in the User's
Guide.

It will be useful to open the following links in your browser to have a look
at the code directly:
`ClimatologyMapOHCAnomaly <https://github.com/MPAS-Dev/MPAS-Analysis/blob/main/mpas_analysis/ocean/climatology_map_ohc_anomaly.py>`_

..
   To do: switch the previous URL to https://github.com/MPAS-Dev/MPAS-Analysis/blob/develop

If you want to be a little more adventurous, you can also pull up the code
for the base class:
`AnalysisTask <https://github.com/MPAS-Dev/MPAS-Analysis/blob/develop/mpas_analysis/shared/analysis_task.py>`_

2.1 Attributes
~~~~~~~~~~~~~~

Classes can contain pieces of data called attributes.  In MPAS-Analysis, the
objects representing tasks share several attributes that they inherit from
the :py:class:`~mpas_analysis.shared.AnalysisTask` class.  A few of the most
important attributes of an analysis task are:

* ``config`` - an object for getting the values of config options

* ``namelist`` - an object for getting namelist options from the E3SM
  simulation

* ``runStreams`` - an object for finding MPAS output files in the ``run``
  directory.  In practice, this is always a restart file used to get the
  MPAS mesh and, for MPAS-Ocean, the vertical coordinate.

* ``historyStreams`` - an object for finding MPAS history streams (often
  ``timeSeriesStatsMonthlyOutput``).

* ``calendar`` - the name of the calendar that was used in the MPAS run
  (in practice always ``'noleap'`` or until recently ``'gregorian_noleap'``).

* ``xmlFileNames`` - a list of XML files associated with plots produced by this
  analysis task.  As we will discuss, these are used to help populate the
  web page showing the analysis.

* ``logger`` - and object that keeps track of sending output to log files
  (rather than the terminal) when the analysis is running.  During the
  ``run_task()`` phase of the analysis when tasks are running in parallel with
  each other, make sure to use ``logger.info()`` instead of ``print()`` to
  send output to the log file.

Within the methods of analysis task class, these attributes can be accessed
using the ``self`` object, e.g. ``self.config``.  It is often helpful to make
a local reference to the object to make the code more compact, e.g.:

.. code-block:: python

    config = self.config
    seasons = config.getexpression('climatologyMapOHCAnomaly', 'seasons')

The analysis task we're looking at, :py:class:`~mpas_analysis.ocean.ClimatologyMapOHCAnomaly`
has some attributes of its own:

* ``mpasClimatologyTask`` - the task that produced the climatology to be
  remapped and plotted

* ``refYearClimatologyTask`` - The task that produced the climatology from the
  first year to be remapped and then subtracted from the main climatology
  (since we want to plot an anomaly from the beginning of the simulation)

2.2 Constructor
~~~~~~~~~~~~~~~

Almost all classes have "constructors", which are methods for making a new
object of that class.  In python, the constructor is called ``__init__()``.
In general, the ``__`` (double underscore) is used in python to indicate a
function or method with special meaning.

The constructor of a subclass (such as
:py:class:`~mpas_analysis.ocean.ClimatologyMapOHCAnomaly`) always calls the
constructor of the superclass (:py:class:`~mpas_analysis.shared.AnalysisTask`
in this case).  So we'll talk about the constructor for
:py:class:`~mpas_analysis.ocean.ClimatologyMapOHCAnomaly` first and then get
to :py:class:`~mpas_analysis.shared.AnalysisTask`.

The constructor for :py:class:`~mpas_analysis.ocean.ClimatologyMapOHCAnomaly`
starts off like this:

.. code-block:: python

    def __init__(self, config, mpas_climatology_task,
                 ref_year_climatology_task, control_config=None):

As with all methods, it takes the ``self`` object as the first argument.
Then, it takes a ``config`` object, which is true of all analysis tasks.  Then,
it has some other arguments that are more specific to the analysis being
performed.  Here, we have 2 other analysis tasks as arguments:
``mpasClimatologyTask`` and ``refYearClimatologyTask``.  As described in
the previous section, these are tasks for computing climatologies that will
later be remapped to a comparison grid for plotting.  A little later in the
constructor, we store references to these tasks as attributes:

.. code-block:: python

        self.mpas_climatology_task = mpas_climatology_task
        self.ref_year_climatology_task = ref_year_climatology_task

Returning to the constructor above, the first thing we do it to call the
super class's ``__init__()`` method:

.. code-block:: python

    def __init__(self, config, mpas_climatology_task,
                 ref_year_climatology_task, control_config=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
            The task that produced the climatology to be remapped and plotted

        ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
            The task that produced the climatology from the first year to be
            remapped and then subtracted from the main climatology

        control_config : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """

        field_name = 'deltaOHC'
        # call the constructor from the base class (AnalysisTask)
        super().__init__(config=config, taskName='climatologyMapOHCAnomaly',
                         componentName='ocean',
                         tags=['climatology', 'horizontalMap', field_name,
                               'publicObs', 'anomaly'])

We're passing along the ``config`` options to the base class so it can store
them.  Then, we're giving the task a unique ``taskName`` (the same as the class
name except that it starts with a lowercase letter).  We're saying that the
MPAS ``componentName`` is the ocean.

Then, we giving the task a number of ``tags`` that can be helpful in
determining whether or not to generate this particular analysis based on the
:ref:`config_generate`.  The tags are used to describe various aspects of the
analysis.  Here, we will produce plots of a ``climatology`` (as opposed to a
time series). The plot will be a ``horizontalMap``.  It will involve the
variable ``deltaOHC``.  This analysis doesn't involve any observations, but we
include a ``publicObs`` tag to indicate that it doesn't require any proprietary
observational data sets that we do not have the rights to make public.
(Currently, we have a few such data sets for things like Antarctic melt rates.)
Finally, the analysis involves an ``anomaly`` computed relative to the
beginning of the simulation.

From there, we get the values of some config options, raising errors if we
find something unexpected:

.. code-block:: python

        section_name = self.taskName

        # read in what seasons we want to plot
        seasons = config.getexpression(section_name, 'seasons')

        if len(seasons) == 0:
            raise ValueError(f'config section {section_name} does not contain '
                             f'valid list of seasons')

        comparison_grid_names = config.getexpression(section_name,
                                                     'comparisonGrids')

        if len(comparison_grid_names) == 0:
            raise ValueError(f'config section {section_name} does not contain '
                             f'valid list of comparison grids')

        depth_ranges = config.getexpression('climatologyMapOHCAnomaly',
                                            'depthRanges',
                                            use_numpyfunc=True)

By default, these config options look like this:

.. code-block:: ini

    [climatologyMapOHCAnomaly]
    ## options related to plotting horizontally remapped climatologies of
    ## ocean heat content (OHC) against control model results (if available)

    ...

    # Months or seasons to plot (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct,
    # Nov, Dec, JFM, AMJ, JAS, OND, ANN)
    seasons =  ['ANN']

    # comparison grid(s) ('latlon', 'antarctic') on which to plot analysis
    comparisonGrids = ['latlon']

    # A list of pairs of minimum and maximum depths (positive up, in meters) to
    # include in the vertical sums.  The default values are the equivalents of the
    # default ranges of the timeSeriesOHCAnomaly task, with a value of -10,000 m
    # intended to be well below the bottom of the ocean for all existing MPAS-O
    # meshes.
    depthRanges = [(0.0, -10000.0), (0.0, -700.0), (-700.0, -2000.0), (-2000.0, -10000.0)]

We plot only the annual mean OHC anomaly and we plot it only on a global
latitude-longitude grid.  The range of depths is:

* the full ocean column

* sea surface to 700 m depth

* 700 m to 2000 m depth

* 2000 m to the seafloor

A user would be free to change any of these config options, and the analysis
should run correctly.  They could choose to plot on a different comparison
grid, add new seasons, or change the depth range.  As long as they ran the
analysis in a fresh directory (or purged output from a previous analysis run),
this should work correctly.

Next, we store some values that will be useful later:

.. code-block:: python

        mpas_field_name = 'deltaOHC'

        variable_list = ['timeMonthly_avg_activeTracers_temperature',
                         'timeMonthly_avg_layerThickness']

This particular analysis involves 4 different depth ranges over which we
compute the ocean heat content.  The remainder of the analysis is performed
separately for each of these depth ranges in subtask.  We loop over the
depth range and add a subtask that will first compute the ocean heat content
(OHC) and then remap it to the comparison grids (``RemapMpasOHCClimatology``):

.. code-block:: python

        for min_depth, max_depth in depth_ranges:
            depth_range_string = \
                f'{np.abs(min_depth):g}-{np.abs(max_depth):g}m'
            remap_climatology_subtask = RemapMpasOHCClimatology(
                mpas_climatology_task=mpas_climatology_task,
                ref_year_climatology_task=ref_year_climatology_task,
                parent_task=self,
                climatology_name=f'{field_name}_{depth_range_string}',
                variable_list=variable_list,
                comparison_grid_names=comparison_grid_names,
                seasons=seasons,
                min_depth=min_depth,
                max_depth=max_depth)

            self.add_subtask(remap_climatology_subtask)

            ...

We will explore the ``RemapMpasOHCClimatology`` subtask later in the tutorial
so we will not discuss it further here.

Still within the loop over depth range, we then add a subtask
(``PlotClimatologyMapSubtask``) for plot we want to create, one for each each
comparison grid and season. (By default, there is only one comparison grid
and one "season": the full year, ``ANN``.)

.. code-block:: python

        for min_depth, max_depth in depth_ranges:
            ...
            out_file_label = f'deltaOHC_{depth_range_string}'
            remap_observations_subtask = None
            if control_config is None:
                ref_title_label = None
                ref_field_name = None
                diff_title_label = 'Model - Observations'

            else:
                control_run_name = control_config.get('runs', 'mainRunName')
                ref_title_label = f'Control: {control_run_name}'
                ref_field_name = mpas_field_name
                diff_title_label = 'Main - Control'

            for comparison_grid_name in comparison_grid_names:
                for season in seasons:
                    # make a new subtask for this season and comparison grid
                    subtask_name = f'plot{season}_{comparison_grid_name}_{depth_range_string}'

                    subtask = PlotClimatologyMapSubtask(
                        self, season, comparison_grid_name,
                        remap_climatology_subtask, remap_observations_subtask,
                        controlConfig=control_config, subtaskName=subtask_name)

                    subtask.set_plot_info(
                        outFileLabel=out_file_label,
                        fieldNameInTitle=f'$\\Delta$OHC over {depth_range_string}',
                        mpasFieldName=mpas_field_name,
                        refFieldName=ref_field_name,
                        refTitleLabel=ref_title_label,
                        diffTitleLabel=diff_title_label,
                        unitsLabel=r'GJ m$^{-2}$',
                        imageCaption=f'Anomaly in Ocean Heat Content over {depth_range_string}',
                        galleryGroup='OHC Anomaly',
                        groupSubtitle=None,
                        groupLink='ohc_anom',
                        galleryName=None)

                    self.add_subtask(subtask)

First, we make sure the subtask has a unique name.  If two tasks or subtasks
have the same ``taskName`` and ``subtaskName``, MPAS-Analysis will only run
the last one and the task manager may become confused.

Then, we create a ``subtask`` object that is an instance of the
:py:class:`~mpas_analysis.ocean.plot_climatology_map_subtask.PlotClimatologyMapSubtask`
class.  This class is shared between several ocean analysis tasks for plotting
climatologies as horizontal maps.  It can plot just MPAS output, remapped to
one or more comparison grids and averaged over one or more seasons.  It can
also plot that data against an observational field that has been remapped to
the same comparison grid and averaged over the same seasons.  In this case,
there are no observations available for comparison
(``remap_observations_subtask = None``).  A user may have provided a
"control" run of MPAS-Analysis to compare with this analysis run (a so-called
"model vs. model" comparison).  If so, ``control_config`` will have config
options describing the other analysis run.  If not, ``control_config`` is
``None``.

Next, We call the
:py:meth:`~mpas_analysis.ocean.plot_climatology_map_subtask.PlotClimatologyMapSubtask.set_plot_info`
method of :py:class:`~mpas_analysis.ocean.plot_climatology_map_subtask.PlotClimatologyMapSubtask`
to provide things like the title and units for the plot and the field to plot.
We also provide information needed for the final analysis web page such as the
name of the gallery group.  (We do not provide a gallery name within the
gallery group because there will be no other galleries within this group.)
All the plots for a given comparison grid will end up in the same gallery,
with different depths and seasons one after the other.

Finally, we call :py:meth:`~mpas_analysis.shared.AnalysisTask.add_subtask()`
to add the ``subtask`` to this task.

2.3 ``setup_and_check()`` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``setup_and_check()`` method of an analysis task is called when it is clear
that this particular analysis has been requested (but before the analysis is
actually ready to run).  This is in contrast to the constructor, which is
run for *every* analysis task everytime MPAS-Analysis runs because we need
information from the analysis task (its name, component and tags) in order to
determine if it should run or not.

In this method, we would typically perform checks to make sure the simulation
has been configured properly to run the analysis.  For example, is the
necessary analysis member enabled.

.. code-block:: python

    def setup_and_check(self):
        """
        Checks whether analysis is being performed only on the reference year,
        in which case the analysis will not be meaningful.

        Raises
        ------
        ValueError: if attempting to analyze only the reference year
        """

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super().setup_and_check()

        start_year, end_year = self.mpas_climatology_task.get_start_and_end()
        ref_start_year, ref_end_year = \
            self.ref_year_climatology_task.get_start_and_end()

        if (start_year == ref_start_year) and (end_year == ref_end_year):
            raise ValueError('OHC Anomaly is not meaningful and will not work '
                             'when climatology and ref year are the same.')

In this particular case, we first call the super class' version of the
:py:meth:`~mpas_analysis.shared.AnalysisTask.setup_and_check()` method.  This
takes care of some important setup.

Then, we use this method to check if the user has specified meaningful values
for the climatology start and end year and the reference year.  If they happen
to be the same, it doesn't really make sense to run the analysis and it will
raise an error so the analysis gets skipped.

The ``ClimatologyMapOHCAnomaly`` has delegated all its work to its subtasks
so it doesn't define a ``run_task()`` method. Tasks or subtasks that actually
do the work typically need to define this method, as we will explore below.

3. Tour of a subtask (``RemapMpasOHCClimatology``)
--------------------------------------------------

The class ``RemapMpasOHCClimatology`` is, in some ways, more complicated than
its "parent" task :py:class:`~mpas_analysis.ocean.ClimatologyMapOHCAnomaly`.
It descends not from the :py:class:`~mpas_analysis.shared.AnalysisTask` base
class but from another subtask,
:py:class:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask`.
This tutorial won't attempt to cover
:py:class:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask` in
all its detail.  The basics are that that class starts with MPAS climatology
data over one or more ``seasons`` that has previously been computed by an
:py:class:`~mpas_analysis.shared.climatology.MpasClimatologyTask` task.  It
remaps that data from the MPAS mesh to one or more comparison grids (e.g.
global latitude-longitude or Antarctic stereographic) where it can be plotted
and compared with observations or another MPAS-Analysis run.

Here, we are not just using
:py:class:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask`
directly because we need to add to its functionality.  We need to compute the
OHC, which is not available straight from MPAS-Ocean output, from the
monthly-mean temperature and layer thickness.

3.1 Attributes
~~~~~~~~~~~~~~

The docstring indicates the attributes that ``RemapMpasOHCClimatology``
includes.  (It also has all the attributes of its super class,
:py:class:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask`,
and that class' super class, :py:class:`~mpas_analysis.shared.AnalysisTask`,
but we don't redundantly document these in the docstring in part because that
would be a maintenance nightmare.)

.. code-block:: python

    class RemapMpasOHCClimatology(RemapMpasClimatologySubtask):
        """
        A subtask for computing climatologies of ocean heat content from
        climatologies of temperature

        Attributes
        ----------
        ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
            The task that produced the climatology from the first year to be
            remapped and then subtracted from the main climatology

        min_depth, max_depth : float
            The minimum and maximum depths for integration
        """

The attributes are a task for computing the climatology over the reference
year (usually the start of the simulation), ``ref_year_climatology_task``,
and the minimum and maximum depth over which the ocean heat content will be
integrated.


3.2 Constructor
~~~~~~~~~~~~~~~

.. code-block:: python

        def __init__(self, mpas_climatology_task, ref_year_climatology_task,
                     parent_task, climatology_name, variable_list, seasons,
                     comparison_grid_names, min_depth, max_depth):

            """
            Construct the analysis task and adds it as a subtask of the
            ``parent_task``.

            Parameters
            ----------
            mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
                The task that produced the climatology to be remapped

            ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
                The task that produced the climatology from the first year to be
                remapped and then subtracted from the main climatology

            parent_task :  mpas_analysis.shared.AnalysisTask
                The parent task, used to get the ``taskName``, ``config`` and
                ``componentName``

            climatology_name : str
                A name that describes the climatology (e.g. a short version of
                the important field(s) in the climatology) used to name the
                subdirectories for each stage of the climatology

            variable_list : list of str
                A list of variable names in ``timeSeriesStatsMonthly`` to be
                included in the climatologies

            seasons : list of str, optional
                A list of seasons (keys in ``shared.constants.monthDictionary``)
                to be computed or ['none'] (not ``None``) if only monthly
                climatologies are needed.

            comparison_grid_names : list of {'latlon', 'antarctic'}
                The name(s) of the comparison grid to use for remapping.

            min_depth, max_depth : float
                The minimum and maximum depths for integration
            """

            depth_range_string = f'{np.abs(min_depth):g}-{np.abs(max_depth):g}m'
            subtask_name = f'remapMpasClimatology_{depth_range_string}'
            # call the constructor from the base class
            # (RemapMpasClimatologySubtask)
            super().__init__(
                mpas_climatology_task, parent_task, climatology_name,
                variable_list, seasons, comparison_grid_names,
                subtaskName=subtask_name)

            self.ref_year_climatology_task = ref_year_climatology_task
            self.run_after(ref_year_climatology_task)
            self.min_depth = min_depth
            self.max_depth = max_depth

Most of the arguments to the constructor are passed along to the constructor
of :py:class:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask`.
These include a reference to the class for computing MPAS climatologies
(used to find the input files and to make sure this task waits until that
task is finished), a reference to the "parent"
:py:class:`~mpas_analysis.ocean.ClimatologyMapOHCAnomaly` task for some of its
attributes, the name of the climatology supplied by the parent (something like
``deltaOHC_0-700m``, depending on the depth range), a list of the variables
that go into computing the OHC, the season(s) over which the climatology was
requested, the comparison grid(s) to plot on and a unique name for this
subtask.

The ``ref_year_climatology_task`` that computing the climatology over the
reference year is retained as an attribute to the class along with
the depth range.  These attributes will all be needed later when we compute the
OHC.  We indicate that this task must wait for the reference climatology to be
available by calling the :py:meth:`~mpas_analysis.shared.AnalysisTask.run_after()`.
The super class will do the same for the ``mpas_climatology_task`` task.  It
will also add this task as a subtask of the parent task.

3.3 ``setup_and_check()`` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As in the parent task, we need to define the ``setup_and_check()`` method.

.. code-block:: python

        def setup_and_check(self):
            """
            Perform steps to set up the analysis and check for errors in the setup.
            """

            # first, call setup_and_check from the base class
            # (RemapMpasClimatologySubtask), which will set up remappers and add
            # variables to mpas_climatology_task
            super().setup_and_check()

            # don't add the variables and seasons to mpas_climatology_task until
            # we're sure this subtask is supposed to run
            self.ref_year_climatology_task.add_variables(self.variableList,
                                                         self.seasons)

In this particular case, we first call the super class' version of the
:py:meth:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask.setup_and_check()`
method.  This takes care of some important setup including adding the variables
and season(s) we need to the ``mpas_climatology_task``.

Then, we use this method to add variables we need
and the requested season(s) to the task for computing the climatology over the
reference year (``ref_year_climatology_task``).  We don't do this in the
constructor because if we did, we would always be asking for the variables
needed to compute the OHC even if we don't actually end up computing it.  This
could be a big waste of time and disk space.  The super class
:py:class:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask` can't
take care of this for us because it isn't designed for computing anomalies,
just "normal" climatologies over a range of years.

.. _tutorial_understand_a_task_subtask_run_task:

3.4 ``run_task()`` method
~~~~~~~~~~~~~~~~~~~~~~~~~

Normally, the main work of a task happens in the ``run_task()`` method.
The ``RemapMpasOHCClimatology`` class doesn't define this method because it is
happy to inherit the
:py:meth:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask.run_task()`
method from its super class,
:py:class:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask`.

An abbreviated version of that method looks like this:

.. code-block:: python

    def run_task(self):
        """
        Compute the requested climatologies
        """
        ...
        for season in self.seasons:
            self._mask_climatologies(season, dsMask)
        ...

It calls a private helper method:

.. code-block:: python

    def _mask_climatologies(self, season, dsMask):
        """
        For each season, creates a masked version of the climatology
        """
        ...
        if not os.path.exists(maskedClimatologyFileName):
            ...

            # customize (if this function has been overridden)
            climatology = self.customize_masked_climatology(climatology,
                                                            season)

            write_netcdf(climatology, maskedClimatologyFileName)

This private method (the leading underscore indicates that it is private), in
turn, calls the ``customize_masked_climatology()`` method, which is our chance
to make changes to the climatology before it gets remapped.  That's where
we will actually compute the OHC from variables available from MPAS output.

3.5 ``customize_masked_climatology()`` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is how we compute the OHC itself:

.. code-block:: python

        def customize_masked_climatology(self, climatology, season):
            """
            Compute the ocean heat content (OHC) anomaly from the temperature
            and layer thickness fields.

            Parameters
            ----------
            climatology : xarray.Dataset
                the climatology data set

            season : str
                The name of the season to be masked

            Returns
            -------
            climatology : xarray.Dataset
                the modified climatology data set
            """

            ohc = self._compute_ohc(climatology)

            ...

We call a private helper method to do the actual work, so let's take a look
at that before we continue with ``customize_masked_climatology()``.

.. code-block:: python

        def _compute_ohc(self, climatology):
            """
            Compute the OHC from the temperature and layer thicknesses in a given
            climatology data sets.
            """
            ds_restart = xr.open_dataset(self.restartFileName)
            ds_restart = ds_restart.isel(Time=0)

            # specific heat [J/(kg*degC)]
            cp = self.namelist.getfloat('config_specific_heat_sea_water')
            # [kg/m3]
            rho = self.namelist.getfloat('config_density0')

            units_scale_factor = 1e-9

            n_vert_levels = ds_restart.sizes['nVertLevels']

            z_mid = compute_zmid(ds_restart.bottomDepth, ds_restart.maxLevelCell-1,
                                 ds_restart.layerThickness)

            vert_index = xr.DataArray.from_dict(
                {'dims': ('nVertLevels',), 'data': np.arange(n_vert_levels)})

            temperature = climatology['timeMonthly_avg_activeTracers_temperature']
            layer_thickness = climatology['timeMonthly_avg_layerThickness']

            masks = [vert_index < ds_restart.maxLevelCell,
                     z_mid <= self.min_depth,
                     z_mid >= self.max_depth]
            for mask in masks:
                temperature = temperature.where(mask)
                layer_thickness = layer_thickness.where(mask)

            ohc = units_scale_factor * rho * cp * layer_thickness * temperature
            ohc = ohc.sum(dim='nVertLevels')
            return ohc

This function uses a combination of mesh information taken from an MPAS
restart file (available from the ``self.restartFileName`` attribute inherited
from :py:class:`~mpas_analysis.shared.climatology.RemapMpasClimatologySubtask`),
namelist options available from the ``self.namelist`` reader (inherited from
:py:class:`~mpas_analysis.shared.AnalysisTask`), and ``temperature`` and
``layer_thickness`` from the ``climatology`` dataset itself.  As the
docstring for ``customize_masked_climatology()`` states, ``climatology`` is
and :py:class:`xarray.Dataset`.  We know it has variables
``timeMonthly_avg_activeTracers_temperature`` and
``timeMonthly_avg_layerThickness`` because we requested them back in the
constructor of :py:class:`~mpas_analysis.ocean.ClimatologyMapOHCAnomaly`.
We compute the ``ohc`` as an :py:class:`xarray.DataArray` that we return from
this helper method.

Back to ``customize_masked_climatology()``, we have:

.. code-block:: python

        def customize_masked_climatology(self, climatology, season):
            ...
            ohc = self._compute_ohc(climatology)

            ref_file_name = self.ref_year_climatology_task.get_file_name(season)
            ref_year_climo = xr.open_dataset(ref_file_name)
            if 'Time' in ref_year_climo.dims:
                ref_year_climo = ref_year_climo.isel(Time=0)
            ref_ohc = self._compute_ohc(ref_year_climo)

            climatology['deltaOHC'] = ohc - ref_ohc
            climatology.deltaOHC.attrs['units'] = 'GJ m^-2$'
            start_year = self.ref_year_climatology_task.startYear
            climatology.deltaOHC.attrs['description'] = \
                f'Anomaly from year {start_year} in ocean heat content'
            climatology = climatology.drop_vars(self.variableList)

            return climatology

We use the same helper function to compute the ``ref_ohc`` using the
climatology for the reference year.  Then, we compute the anomaly (the
difference between these two, ``deltaOHC``) and we add some attributes,
``units`` and ``description``, to make the NetCDF output that will go into the
analysis output directory a little more useful.

4. The full code for posterity
------------------------------

Since the ``ClimatologyMapOHCAnomaly`` analysis task may evolve in the future,
here is the full analysis task as described in this tutorial:

.. code-block:: python

    # This software is open source software available under the BSD-3 license.
    #
    # Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
    # Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
    # reserved.
    # Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
    #
    # Additional copyright and license information can be found in the LICENSE file
    # distributed with this code, or at
    # https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
    import xarray as xr
    import numpy as np

    from mpas_analysis.shared import AnalysisTask
    from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask
    from mpas_analysis.ocean.plot_climatology_map_subtask import \
        PlotClimatologyMapSubtask
    from mpas_analysis.ocean.utility import compute_zmid


    class ClimatologyMapOHCAnomaly(AnalysisTask):
        """
        An analysis task for comparison of the anomaly from a reference year
        (typically the start of the simulation) of ocean heat content (OHC)

        Attributes
        ----------
        mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
            The task that produced the climatology to be remapped and plotted

        ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
            The task that produced the climatology from the first year to be
            remapped and then subtracted from the main climatology
        """

        def __init__(self, config, mpas_climatology_task,
                     ref_year_climatology_task, control_config=None):
            """
            Construct the analysis task.

            Parameters
            ----------
            config : mpas_tools.config.MpasConfigParser
                Configuration options

            mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
                The task that produced the climatology to be remapped and plotted

            ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
                The task that produced the climatology from the first year to be
                remapped and then subtracted from the main climatology

            control_config : mpas_tools.config.MpasConfigParser, optional
                Configuration options for a control run (if any)
            """

            field_name = 'deltaOHC'
            # call the constructor from the base class (AnalysisTask)
            super().__init__(config=config, taskName='climatologyMapOHCAnomaly',
                             componentName='ocean',
                             tags=['climatology', 'horizontalMap', field_name,
                                   'publicObs', 'anomaly'])

            self.mpas_climatology_task = mpas_climatology_task
            self.ref_year_climatology_task = ref_year_climatology_task

            section_name = self.taskName

            # read in what seasons we want to plot
            seasons = config.getexpression(section_name, 'seasons')

            if len(seasons) == 0:
                raise ValueError(f'config section {section_name} does not contain '
                                 f'valid list of seasons')

            comparison_grid_names = config.getexpression(section_name,
                                                         'comparisonGrids')

            if len(comparison_grid_names) == 0:
                raise ValueError(f'config section {section_name} does not contain '
                                 f'valid list of comparison grids')

            depth_ranges = config.getexpression('climatologyMapOHCAnomaly',
                                                'depthRanges',
                                                use_numpyfunc=True)

            mpas_field_name = 'deltaOHC'

            variable_list = ['timeMonthly_avg_activeTracers_temperature',
                             'timeMonthly_avg_layerThickness']

            for min_depth, max_depth in depth_ranges:
                depth_range_string = \
                    f'{np.abs(min_depth):g}-{np.abs(max_depth):g}m'
                remap_climatology_subtask = RemapMpasOHCClimatology(
                    mpas_climatology_task=mpas_climatology_task,
                    ref_year_climatology_task=ref_year_climatology_task,
                    parent_task=self,
                    climatology_name=f'{field_name}_{depth_range_string}',
                    variable_list=variable_list,
                    comparison_grid_names=comparison_grid_names,
                    seasons=seasons,
                    min_depth=min_depth,
                    max_depth=max_depth)

                self.add_subtask(remap_climatology_subtask)

                out_file_label = f'deltaOHC_{depth_range_string}'
                remap_observations_subtask = None
                if control_config is None:
                    ref_title_label = None
                    ref_field_name = None
                    diff_title_label = 'Model - Observations'

                else:
                    control_run_name = control_config.get('runs', 'mainRunName')
                    ref_title_label = f'Control: {control_run_name}'
                    ref_field_name = mpas_field_name
                    diff_title_label = 'Main - Control'

                for comparison_grid_name in comparison_grid_names:
                    for season in seasons:
                        # make a new subtask for this season and comparison grid
                        subtask_name = f'plot{season}_{comparison_grid_name}_{depth_range_string}'

                        subtask = PlotClimatologyMapSubtask(
                            self, season, comparison_grid_name,
                            remap_climatology_subtask, remap_observations_subtask,
                            controlConfig=control_config, subtaskName=subtask_name)

                        subtask.set_plot_info(
                            outFileLabel=out_file_label,
                            fieldNameInTitle=f'$\\Delta$OHC over {depth_range_string}',
                            mpasFieldName=mpas_field_name,
                            refFieldName=ref_field_name,
                            refTitleLabel=ref_title_label,
                            diffTitleLabel=diff_title_label,
                            unitsLabel=r'GJ m$^{-2}$',
                            imageCaption=f'Anomaly in Ocean Heat Content over {depth_range_string}',
                            galleryGroup='OHC Anomaly',
                            groupSubtitle=None,
                            groupLink='ohc_anom',
                            galleryName=None)

                        self.add_subtask(subtask)

        def setup_and_check(self):
            """
            Checks whether analysis is being performed only on the reference year,
            in which case the analysis will not be meaningful.

            Raises
            ------
            ValueError: if attempting to analyze only the reference year
            """

            # first, call setup_and_check from the base class (AnalysisTask),
            # which will perform some common setup, including storing:
            #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
            #     self.namelist, self.runStreams, self.historyStreams,
            #     self.calendar
            super().setup_and_check()

            start_year, end_year = self.mpas_climatology_task.get_start_and_end()
            ref_start_year, ref_end_year = \
                self.ref_year_climatology_task.get_start_and_end()

            if (start_year == ref_start_year) and (end_year == ref_end_year):
                raise ValueError('OHC Anomaly is not meaningful and will not work '
                                 'when climatology and ref year are the same.')


    class RemapMpasOHCClimatology(RemapMpasClimatologySubtask):
        """
        A subtask for computing climatologies of ocean heat content from
        climatologies of temperature

        Attributes
        ----------
        ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
            The task that produced the climatology from the first year to be
            remapped and then subtracted from the main climatology

        min_depth, max_depth : float
            The minimum and maximum depths for integration
        """

        def __init__(self, mpas_climatology_task, ref_year_climatology_task,
                     parent_task, climatology_name, variable_list, seasons,
                     comparison_grid_names, min_depth, max_depth):

            """
            Construct the analysis task and adds it as a subtask of the
            ``parent_task``.

            Parameters
            ----------
            mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
                The task that produced the climatology to be remapped

            ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
                The task that produced the climatology from the first year to be
                remapped and then subtracted from the main climatology

            parent_task :  mpas_analysis.shared.AnalysisTask
                The parent task, used to get the ``taskName``, ``config`` and
                ``componentName``

            climatology_name : str
                A name that describes the climatology (e.g. a short version of
                the important field(s) in the climatology) used to name the
                subdirectories for each stage of the climatology

            variable_list : list of str
                A list of variable names in ``timeSeriesStatsMonthly`` to be
                included in the climatologies

            seasons : list of str, optional
                A list of seasons (keys in ``shared.constants.monthDictionary``)
                to be computed or ['none'] (not ``None``) if only monthly
                climatologies are needed.

            comparison_grid_names : list of {'latlon', 'antarctic'}
                The name(s) of the comparison grid to use for remapping.

            min_depth, max_depth : float
                The minimum and maximum depths for integration
            """

            depth_range_string = f'{np.abs(min_depth):g}-{np.abs(max_depth):g}m'
            subtask_name = f'remapMpasClimatology_{depth_range_string}'
            # call the constructor from the base class
            # (RemapMpasClimatologySubtask)
            super().__init__(
                mpas_climatology_task, parent_task, climatology_name,
                variable_list, seasons, comparison_grid_names,
                subtaskName=subtask_name)

            self.ref_year_climatology_task = ref_year_climatology_task
            self.run_after(ref_year_climatology_task)
            self.min_depth = min_depth
            self.max_depth = max_depth

        def setup_and_check(self):
            """
            Perform steps to set up the analysis and check for errors in the setup.
            """

            # first, call setup_and_check from the base class
            # (RemapMpasClimatologySubtask), which will set up remappers and add
            # variables to mpas_climatology_task
            super().setup_and_check()

            # don't add the variables and seasons to mpas_climatology_task until
            # we're sure this subtask is supposed to run
            self.ref_year_climatology_task.add_variables(self.variableList,
                                                         self.seasons)

        def customize_masked_climatology(self, climatology, season):
            """
            Compute the ocean heat content (OHC) anomaly from the temperature
            and layer thickness fields.

            Parameters
            ----------
            climatology : xarray.Dataset
                the climatology data set

            season : str
                The name of the season to be masked

            Returns
            -------
            climatology : xarray.Dataset
                the modified climatology data set
            """

            ohc = self._compute_ohc(climatology)
            ref_file_name = self.ref_year_climatology_task.get_file_name(season)
            ref_year_climo = xr.open_dataset(ref_file_name)
            if 'Time' in ref_year_climo.dims:
                ref_year_climo = ref_year_climo.isel(Time=0)
            ref_ohc = self._compute_ohc(ref_year_climo)

            climatology['deltaOHC'] = ohc - ref_ohc
            climatology.deltaOHC.attrs['units'] = 'GJ m^-2'
            start_year = self.ref_year_climatology_task.startYear
            climatology.deltaOHC.attrs['description'] = \
                f'Anomaly from year {start_year} in ocean heat content'
            climatology = climatology.drop_vars(self.variableList)

            return climatology

        def _compute_ohc(self, climatology):
            """
            Compute the OHC from the temperature and layer thicknesses in a given
            climatology data sets.
            """
            ds_restart = xr.open_dataset(self.restartFileName)
            ds_restart = ds_restart.isel(Time=0)

            # specific heat [J/(kg*degC)]
            cp = self.namelist.getfloat('config_specific_heat_sea_water')
            # [kg/m3]
            rho = self.namelist.getfloat('config_density0')

            units_scale_factor = 1e-9

            n_vert_levels = ds_restart.sizes['nVertLevels']

            z_mid = compute_zmid(ds_restart.bottomDepth, ds_restart.maxLevelCell-1,
                                 ds_restart.layerThickness)

            vert_index = xr.DataArray.from_dict(
                {'dims': ('nVertLevels',), 'data': np.arange(n_vert_levels)})

            temperature = climatology['timeMonthly_avg_activeTracers_temperature']
            layer_thickness = climatology['timeMonthly_avg_layerThickness']

            masks = [vert_index < ds_restart.maxLevelCell,
                     z_mid <= self.min_depth,
                     z_mid >= self.max_depth]
            for mask in masks:
                temperature = temperature.where(mask)
                layer_thickness = layer_thickness.where(mask)

            ohc = units_scale_factor * rho * cp * layer_thickness * temperature
            ohc = ohc.sum(dim='nVertLevels')
            return ohc
