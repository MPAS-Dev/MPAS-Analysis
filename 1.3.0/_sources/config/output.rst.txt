.. |n~| unicode:: U+00F1
   :trim:

.. _config_output:

Output
======

The ``[output]`` section of a configuration file contains options used to
specify the location of the "main" simulation::

  [output]
  ## options related to writing out plots, intermediate cached data sets, logs,
  ## etc.

  # directory where analysis should be written
  # NOTE: This directory path must be specific to each test case.
  baseDirectory = /dir/for/analysis/output

  # subdirectories within baseDirectory for analysis output
  plotsSubdirectory = plots
  logsSubdirectory = logs
  mpasClimatologySubdirectory = clim/mpas
  mappingSubdirectory = mapping
  timeSeriesSubdirectory = timeseries
  # provide an absolute path to put HTML in an alternative location (e.g. a web
  # portal)
  htmlSubdirectory = html

  # a list of analyses to generate.  Valid names can be seen by running:
  #   mpas_analysis --list
  # This command also lists tags for each analysis.
  # Shortcuts exist to generate (or not generate) several types of analysis.
  # These include:
  #   'all' -- all analyses will be run
  #   'all_publicObs' -- all analyses for which observations are availabe on the
  #                      public server (the default)
  #   'all_<tag>' -- all analysis with a particular tag will be run
  #   'all_<component>' -- all analyses from a given component (either 'ocean'
  #                        or 'seaIce') will be run
  #   'only_<component>', 'only_<tag>' -- all analysis from this component or
  #                                       with this tag will be run, and all
  #                                       analysis for other components or
  #                                       without the tag will be skipped
  #   'no_<task_name>' -- skip the given task
  #   'no_<component>', 'no_<tag>' -- in analogy to 'all_*', skip all analysis
  #                                   tasks from the given compoonent or with
  #                                   the given tag.  Do
  #                                      mpas_analysis --list
  #                                   to list all task names and their tags
  # an equivalent syntax can be used on the command line to override this
  # option:
  #    mpas_analysis config.analysis --generate \
  #         only_ocean,no_timeSeries,timeSeriesSST
  generate = ['all_publicObs']

Output Directories
------------------

The output directories are specified through a base directory and a set of
subdirectories within that base.  You will always need to set ``baseDirectory``
in a custom configuration file.

The subdirectories of output can be renamed if desired or an absolute path
can be specified if output to another location is desired.  This is
particularly relevant to ``htmlSubdirectory``, which may be pointed to a
a public space such as a web portal for display on the web.  Note:
MPAS-Analysis does not change the HTML output to be world readable so you
will need to do this manually after a run has completed (or inside of a job
script) to see the results on a public web page.

.. _config_generate:

Generate Option
---------------

The ``generate`` option is used to control which tasks run.  The simplest
choice is::

  generate = ['all']

in which case MPAS-Analysis will attempt to run all analysis tasks.  In this
mode, some tasks may fail with a warning printed to the screen during their
:py:meth:`~mpas_analysis.shared.analysis_task.AnalysisTask.setup_and_check`
phase if the simulation was not configured to support that task.  All tasks
that pass the
:py:meth:`~mpas_analysis.shared.analysis_task.AnalysisTask.setup_and_check`
phase will be run.

The next simplest option is to specify a single task name::

  generate = ['climatologyMapSST']

or a list of task names::

  generate = ['climatologyMapSST', 'climatologyMapSSS']

in which case only the listed tasks are run.

A third way to determine which tasks to generate is to make use of "tags" for
each task.  To see what tags each task has, run::

  mpas_analysis --list

This will show all available tasks together with the component they belong to
and the tags for each.  To run only those analysis tasks with a particular tag,
set, for example ``climatology``::

  generate = ['all_climatology']

This will generate only those tasks that make use of climatologies.

A useful tag for the public release of MPAS-Analysis is the ``publicObs`` tag,
which is found on all tasks that will run successfully if you have downloaded
the observations from the `E3SM public data repository`_.  Some MPAS-Analysis
tasks make use of data sets that are only available after registering with a
data portal or by contacting the authors of that data set directly, so that
these data have not been included in the data repository.  The default is to
generate only the tasks with observations in in the repository::

  generate = ['all_publicObs']

The names of components (``ocean`` or ``seaIce``) can also be used as tags.

There are also ways to specify that a given tag should not be present
(``no_<tag>``) or that only analysis with the given tag should be run
(``only_<tag>``).  These options are useful when combined in a series with
other generate options.  For example, to generate all tasks with publicly
available observation except those for the ``seaIce`` component, you could
specify::

  generate = ['all_publicObs', 'no_seaIce']

To specify that you wanted to plot climatologies from the ocean component, you
could use::

  generate = ['all_publicObs', 'only_climatologies', 'only_ocean']

If you wanted to plot all tasks with publicly available data sets that used
either climatologies or time series, you could use::

  generate = ['all_climatologies', 'all_timeSeries', 'only_publicObs']

Finally, we note that the ``generate`` option in the configuration file can
be overridden by specifying the ``--generate`` option on the command line::

  mpas_analysis --generate=all_publicObs,no_index,no_climatologyMapSST \
      config.my_run

This example would override whatever ``generate`` option was specified in
``config.my_run`` with a directive to generate only tasks that support the
publicly available observations, skipping those using climate indices (e.g.
El Ni |n~| o 3.4) and also skipping ``climatologyMapSST``.


.. _`E3SM public data repository`: https://web.lcrc.anl.gov/public/e3sm/diagnostics/