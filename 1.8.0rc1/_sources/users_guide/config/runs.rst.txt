.. _config_runs:

Runs
====

The ``[runs]`` section of a configuration file contains options used to name
the "main" run, a preprocessed E3SM v0 run (if any) and to point to analysis
of a control E3SM v1 or standalone MPAS run (if any)::

  [runs]
  ## options related to the run to be analyzed and control runs to be
  ## compared against

  # mainRunName is a name that identifies the simulation being analyzed.
  mainRunName = runName

  # preprocessedReferenceRunName is the name of a control run that has been
  # preprocessed to compare against (or None to turn off comparison).  Reference
  # runs of this type would have preprocessed results because they were not
  # performed with MPAS components (so they cannot be easily ingested by
  # MPAS-Analysis)
  preprocessedReferenceRunName = None

  # config file for a control run to which this run will be compared.  The
  # analysis should have already been run to completion once with this config
  # file, so that the relevant MPAS climatologies already exist and have been
  # remapped to the comparison grid.  Leave this option commented out if no
  # control run is desired.
  # controlRunConfigFile = /path/to/config/file

  # config file for a main run on which the analysis was already run to
  # completion.  The relevant MPAS climatologies already exist and have been
  # remapped to the comparison grid and time series have been extracted.
  # Leave this option commented out if the analysis for the main run should be
  # performed.
  # mainRunConfigFile = /path/to/config/file

The name of the "main" run (as opposed to a control run, if any) can be any
identifier that will be used in figure titles, legends, web pages and file
names to identify this run.  It does not need to be the name of the simulation
as specified in E3SM::

  mainRunName = runName

A few of the time series plots in MPAS-Analysis can be compared against a
preprocessed control run from E3SM v0 (which was similar to the CESM, the
Community Earth System Model).  If these data are available and the comparison
to these runs is desired, the name of the control run should be specified
here and the paths to the data set should be specified (see
:ref:`config_preprocessed`).  If not this name should be left as ``None``::

  preprocessedReferenceRunName = None

MPAS-Analysis supports comparison between the "main" run and a control run
from either E3SM or a standalone MPAS component.  By default, this feature is
disabled by commenting out the configuration option::

  # controlRunConfigFile = /path/to/config/file

To specify a control run, first run MPAS analysis on the control run.  Be
sure that:

  * the start and end year for climatologies, time series and climate indices
    is covered by the simulation output.
  * most configuration options for the control run are the same as for the
    main run.  The exceptions are contents of the ``[run]``, ``[input]`` and
    ``[output]`` sections.  The range of years for climatologies can be
    different, but this is discouraged.

Once the analysis has been run on the control run, a comparison is made by
uncommenting ``controlRunConfigFile`` and specifying the path to the
configuration file use in this analysis, e.g.::

  controlRunConfigFile = config.control_run

If analysis has already been run on the "main" run in a "main vs ref"
comparison, some time can be saved in performing the comparison
(particularly for higher resolution output, for which a lot of the
computation time goes into computing climatologies and extracting time
series).  By default, this feature is disabled by commenting out the
configuration option::

  # mainRunConfigFile = /path/to/config/file

To specify a main run, first run MPAS analysis on the main run.  The
"comparison" config file should be nearly identical to the "main" config
file except that:

  * The output ``baseDirectory`` should be different.
  * the start and end year for climatologies, time series and climate indices
    must be the actual range used if output data was not available to span
    the requested range in the "main" run.

All configuration information for the "main" run in the "main vs ref"
comparison is taken from the "comparsion" config file, not the "main" config.
Only the output directories and subdirectories for climatologies, time series,
mapping files and mask files (if these latter 2 were generated on the fly)
will be taken from the "main" config file.  Symbolic links will be made to
these directories so the comparison analysis run can reuse this data.
Specify the path to the configuration file use in "main" analysis by
uncommenting the option and providing a relative or absolute path to the
config file::

  mainRunConfigFile = config.main_run


