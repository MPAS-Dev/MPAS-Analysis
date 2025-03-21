.. _tutorial_getting_started:

User: Getting Started
=====================

This tutorial walks a new user through the basics of using MPAS-Analysis.

.. warning::

   Beware of hyphens and underscores. The conda package name and some
   environment names use a hyphen (``mpas-analysis``). The python package and
   the command name use an underscore (``mpas_analysis``).

1 Setting up a Conda Environment
---------------------------------

MPAS-Analysis relies on several packages that are only available as conda
packages from the ``conda-forge`` channel.  The first step for running
MPAS-Analysis is to create a conda environment with all the needed packages.

1.1 Installing Miniforge3
~~~~~~~~~~~~~~~~~~~~~~~~~

If you have not yet installed Anaconda, Miniconda or Miniforge, you will need
to begin there.  The concept behind Anaconda is that just about everything you
would need for a typical python workflow is included.  The concept behind
Miniconda and Miniforg3 is that you create different environments for different
purposes.  This allows for greater flexibility and tends to lead to fewer
conflicts between incompatible packages, particularly when using a channel
other than the ``defaults`` supplied by Anaconda.  Since we will use the
``conda-forge`` channel, the Miniforge approach is strongly recommended.

First download the `Miniforge3 installer`_ for your operating system, then run
it:

.. code-block:: bash

   $ /bin/bash Miniforge3-Linux-x86_64.sh

.. note::

   MPAS-Analysis and many of the packages it depends on support OSX and Linux
   but not Windows.  Typically, HPC machines need the Linux-x86_64 version.

In this tutorial, we assume that Miniforge3 is installed in the default
location, ``~/miniforge3``.  If you choose to install it somewhere else, just
make sure to make the appropriate substitution whenever you see a reference to
this path below.

You will see prompt like this::

   Do you wish the installer to initialize Miniforge3
   by running conda init? [yes|no]

You may wish to skip the step (answer "no") if you are working on a system
where you will also be using other conda environments.  This means you have to
manually activate ``conda`` whenever you need it.  For ``bash`` and similar
shells, this is:

.. code-block:: bash

   $ source ~/miniforge3/etc/profile.d/conda.sh

If you use ``csh``, ``tcsh`` or related shells, this becomes:

.. code-block:: csh

   > source ~/miniforge3/etc/profile.d/conda.csh

1.2 Creating a conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we need to add the `conda-forge channel`_ and make sure it always takes
precedence for packages available on that channel:

.. code-block:: bash

   $ conda config --add channels conda-forge
   $ conda config --set channel_priority strict

Then, you can create a new conda environment called ``mpas-analysis`` and
install the latest version of the ``mpas-analysis`` package into it:

.. code-block:: bash

   $ conda create -n mpas-analysis python=3.13 mpas-analysis

1.3 Activating the environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each time you open a new terminal window, to activate the ``mpas-analysis``
environment, you will need to run either for ``bash``:

.. code-block:: bash

   $ source ~/miniforge3/etc/profile.d/conda.sh
   $ conda activate mpas-analysis

or for ``csh``:

.. code-block:: csh

   > source ~/miniforge3/etc/profile.d/conda.csh
   > conda activate mpas-analysis

You can skip the ``source`` command if you chose to initialize Miniforge3 so it
loads automatically.

2 Downloading observations
---------------------------

MPAS-Analysis can only run a very limited set of analysis without observational
datasets for comparison.  Many analysis tasks also require "mapping files" that
describe interpolation between different grids and meshes.  Some tasks also
need masks that define regions of interest on different grids and meshes.

Most of the observational datasets as well as some common mapping files and
region masks are available for public download by calling:

.. code-block:: bash

   $ download_analysis_data -o </path/to/mpas_analysis/diagnostics>

Substitute ``</path/to/mpas_analysis/diagnostics>`` with a convenient location
for storing the diagnostics data.  You will need this path later on when you
set up a configuration file for running the analysis.

.. note::
   The observational data is approximately **45GB** in size and may take some
   time to download.  If the download gets interrupted, only incomplete or
   missing files will be downloaded if you run ``download_analysis_data`` again.

.. note::

   If you are on an E3SM supported system, the data has already
   been downloaded for you to use.  Please look at the ``baseDirectory`` config
   option in the ``[diagnostics]`` of an example `config file`_
   for your supported machine.

.. note::

   We do not currently support downloading a subset of the data but we plan to
   add this capability in the future.  In a pinch, you could manually download
   files.  For reference, here is the `full list of files`_ that are
   automatically downloaded from the `diagnostics web server`_.


3 Downloading an example run
-----------------------------

We provide data for an "ultra-low-resolution" E3SM simulation (480 km grid cells
for the ocean and sea-ice components) as a relatively small example dataset for
testing MPAS-Analysis.  The simulation is too coarse and too short to be of
scientific interest, but we wanted to keep the size manageable.  The ocean and
sea-ice output are about 950 MB.

Download the `Ultra-low-res ocean and sea-ice dataset`_ and untar it into a
convenient directory:

.. code-block:: bash

   $ mkdir A_WCYCL1850.ne4_oQU480.anvil
   $ cd A_WCYCL1850.ne4_oQU480.anvil
   $ tar -xvf ../20200305.A_WCYCL1850.ne4_oQU480.anvil.ocean_seaice.tar.gz
   $ cd ..

4 Configuring MPAS-Analysis
----------------------------

Configuration is with Python ``cfg`` (also called ``ini``) files:

.. code-block:: ini

   [runs]
   # mainRunName is a name that identifies the simulation being analyzed.
   mainRunName = runName

   [execute]
   # the number of parallel tasks (1 means tasks run in serial, the default)
   parallelTaskCount = 1
   # the parallelism mode in ncclimo ("serial" or "bck")
   ncclimoParallelMode = serial
   ...

The `default config file`_ contains over 1,000 config options, which gives a lot
of flexibility to MPAS-Analysis but can be more than bit overwhelming to new
users.

You can modify the default options with one or more custom config files.  We
provide an `example config file`_ as a starting point. This file contains the
most common config options that a user might want to customize.  The values are
mostly the same as in the `default config file`_.

If you are on an E3SM supported machine, you can use the `E3SM example config file`_ instead.  If you specify the name of the supported machine with the
``--machine`` flag when you call ``mpas_analysis``, there are several config
options that will be set for you automatically.  E3SM supported machines aren't
the primary focus of this tutorial, the following will assume you're working
on an unknown machine (or simply not taking advantage of known paths on a
supported machine).

First, you should copy this file to a new name for a specific run (say
``myrun.cfg``).  Then, you should modify any config options you want to change
in your new config file. At a minimum, you need to specify:

* ``[runs]/mainRunName``:  A name for the run to be included plot titles
  and legends
* ``[diagnostics]/base_path``: The base directory for observations,
  mapping files and region masks
* ``[input]/baseDirectory``: The directory for the simulation results
  to analyze
* ``[input]/mpasMeshName``: The name of the MPAS ocean/sea ice mesh
* ``[output]/baseDirectory``: The directory for the analysis results

We will cover these and a few other common options in this tutorial.  With the
exception of a few paths that you will need to provide, the config options
displayed below are the ones appropriate for the example E3SM simulation from
Section 3.

4.1 [runs]
~~~~~~~~~~

The ``[runs]`` section contains options related to which E3SM simulation(s) are
being analyzed:

.. code-block:: ini

    [runs]
    ## options related to the run to be analyzed and control runs to be
    ## compared against

    # mainRunName is a name that identifies the simulation being analyzed.
    mainRunName = A_WCYCL1850.ne4_oQU480.anvil

The ``mainRunName`` can be any useful name that will appear at the top of each
web page of the analysis output and in the legends or titles of the figures.
Often, this is the full name of the E3SM simulation but sometimes it is
convenient to have a shorter name.  In this case, we use part of the run name
but leave off the date of the simulation to keep it a little shorter.

4.2 [execute]
~~~~~~~~~~~~~

The ``[execute]`` section contains options related to serial or parallel
execution of the individual "tasks" that make up an MPAS-Analysis run.

.. code-block:: ini

    [execute]
    ## options related to executing parallel tasks

    # the number of parallel tasks (1 means tasks run in serial, the default)
    parallelTaskCount = 4

    # the parallelism mode in ncclimo ("serial", "bck" or "mpi")
    # Set this to "bck" (background parallelism) if running on a machine that can
    # handle 12 simultaneous processes, one for each monthly climatology.
    # Set to "mpi" to run one MPI task on each node and however many threads per
    # node to reach 12 total threads.
    ncclimoParallelMode = bck

``parallelTaskCount`` should typically be equal to the number of cores on a
laptop or a desktop.  In a high-performance computing (HPC) node, it is
typically not possible to run 30 to 60 tasks in parallel even if the node has
that many cores.  This is because many tasks consume enough memory that the
job will crash with too many tasks running at the same time.  We have found that
``parallelTaskCount`` should typically be somewhere between 6 and 12 for the
HPC machines we use for E3SM.

``ncclimoParallelMode`` indicates how `ncclimo`_ should be run to make
climatologies used in many MPAS-Analysis plots.  Typically, we recommend
``bck``, meaning ``ncclimo`` runs with 12 threads at once on the same node.
In circumstances where ``ncclimo`` is crashing and it appears to be running out
of memory, it is worth exploring ``serial`` or ``mpi`` modes, or using the
``xarray`` and ``dask`` instead to compute climatologies by setting
``[climatology]/useNcclimo = False``

For this tutorial, we suggest starting with 4 parallel tasks and ``ncclimo`` in
``bck`` mode.

4.3 [diagnostics]
~~~~~~~~~~~~~~~~~

The ``diagnostics`` section is used to supply the directory where you downloaded
observations in Section 2.

.. code-block:: ini

    [diagnostics]
    ## config options related to observations, mapping files and region files used
    ## by MPAS-Analysis in diagnostics computations.

    # The base path to the diagnostics directory.  Typically, this will be a shared
    # directory on each E3SM supported machine (see the example config files for
    # its location).  For other machines, this would be the directory pointed to
    # when running "download_analysis_data.py" to get the public observations,
    # mapping files and region files.
    base_path = /path/to/diagnostics

For ``base_path``, supply the path where you downloaded the data
``</path/to/mpas_analysis/diagnostics>``.

4.4 [input]
~~~~~~~~~~~

The ``[input]`` section provides paths to the E3SM simulation data and the name
of the MPAS-Ocean and MPAS-Seaice mesh.

.. code-block:: ini

    [input]
    ## options related to reading in the results to be analyzed

    # directory containing model results
    baseDirectory = /dir/for/model/output

    # Note: an absolute path can be supplied for any of these subdirectories.
    # A relative path is assumed to be relative to baseDirectory.
    # In this example, results are assumed to be in <baseDirecory>/run

    # subdirectory containing restart files
    runSubdirectory = run
    # subdirectory for ocean history files
    oceanHistorySubdirectory = archive/ocn/hist
    # subdirectory for sea ice history files
    seaIceHistorySubdirectory = archive/ice/hist

    # names of namelist and streams files, either a path relative to baseDirectory
    # or an absolute path.
    oceanNamelistFileName = run/mpaso_in
    oceanStreamsFileName = run/streams.ocean
    seaIceNamelistFileName = run/mpassi_in
    seaIceStreamsFileName = run/streams.seaice

    # name of the ocean and sea-ice mesh (e.g. EC30to60E2r2, WC14to60E2r3,
    # ECwISC30to60E2r1, SOwISC12to60E2r4, oQU240, etc.)
    mpasMeshName = oQU480

The ``baseDirectory`` is the path where you untarred the example run.

The ``mpasMeshName`` is the standard E3SM name for the MPAS-Ocean and
MPAS-Seaice mesh.  In this example, this is ``oQU480``, meaning the
quasi-uniform 480-km mesh for the ocean and sea ice.

The ``runSubdirectory`` must contain valid MPAS-Ocean and MPAS-Seaice restart
files, used to get information about the MPAS mesh and the ocean vertical grid.

The ``oceanHistorySubdirectory`` must contain MPAS-Ocean monthly mean output
files, typically named::

   mpaso.hist.am.timeSeriesStatsMonthly.YYYY-MM-DD.nc

Similarly, ``seaIceHistorySubdirectory`` contains the MPAS-Seaice monthly mean
output::

   mpassi.hist.am.timeSeriesStatsMonthly.YYYY-MM-DD.nc

Finally, MPAS-Analysis needs a set of "namelists" and "streams" files that
provide information on the E3SM configuration for MPAS-Ocean and MPAS-Seaice,
and about the output files, respectively.  These are typically also found in
the ``run`` directory.

For the example data, only ``baseDirectory`` and ``mpasMeshName`` need to be
set, the other options can be left as the defaults from the
`example config file`_.

For this tutorial, you just need to set ``baseDirectory`` to the place where
you untarred the simulation results.  The other config options should be as
they are in the code block above.

4.5 [output]
~~~~~~~~~~~~

The ``[output]`` section provides a path where the output from the analysis run
will be written, the option to output the results web pages to another
location, and a list of analysis to be generated (or explicitly skipped).

.. code-block:: ini

    [output]
    ## options related to writing out plots, intermediate cached data sets, logs,
    ## etc.

    # directory where analysis should be written
    # NOTE: This directory path must be specific to each test case.
    baseDirectory = /dir/for/analysis/output

    # provide an absolute path to put HTML in an alternative location (e.g. a web
    # portal)
    htmlSubdirectory = html

    # a list of analyses to generate.  Valid names can be seen by running:
    #   mpas_analysis --list
    # This command also lists tags for each analysis.
    # Shortcuts exist to generate (or not generate) several types of analysis.
    # These include:
    #   'all' -- all analyses will be run
    #   'all_publicObs' -- all analyses for which observations are available on the
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
    #                                   tasks from the given component or with
    #                                   the given tag.  Do
    #                                      mpas_analysis --list
    #                                   to list all task names and their tags
    # an equivalent syntax can be used on the command line to override this
    # option:
    #    mpas_analysis analysis.cfg --generate \
    #         only_ocean,no_timeSeries,timeSeriesSST
    generate = ['all_publicObs']

``baseDirectory`` is any convenient location for the output.

``htmlSubdirectory`` can simply be the ``<baseDirectory>/html``, the default or
an absolute path to another location.  The later is useful for HPC machines that
have a web portal.

Finally, the ``generate`` option provides a python list of flags that can be
used to determine which analysis will be generated.  For this tutorial, we will
stick with the default, ``'all_publicObs'``, indicating that we will only run
analysis where the observations are included on the public server and which
were downloaded in Section 2 (or analysis that does not require observations).

4.6 [climatology], [timeSeries] and [index]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These options determine the start and end years of climatologies (time averages
over a particular month, season or the full year), time series or the El Niño
climate index.

.. code-block:: ini

    [climatology]
    ## options related to producing climatologies, typically to compare against
    ## observations and previous runs

    # the first year over which to average climatalogies
    startYear = 3
    # the last year over which to average climatalogies
    endYear = 5

    [timeSeries]
    ## options related to producing time series plots, often to compare against
    ## observations and previous runs

    # start and end years for timeseries analysis.  Out-of-bounds values will lead
    # to an error.
    startYear = 1
    endYear = 5

    [index]
    ## options related to producing nino index.

    # start and end years for timeseries analysis.  Out-of-bounds values will lead
    # to an error.
    startYear = 1
    endYear = 5

For each of these, options a full year of data must exist for that year to
be included in the analysis.

For the example E3SM simulation that we downloaded in Section 3, only 5 years of
simulation data are available, so we are doing a climatology over the last 3
years (3 to 5) and displaying time series and the El Niño index over the full
5 years.

5 Running MPAS-Analysis
-----------------------

The hard work is done.  Now that we have a config file, we are ready to run.

On many file systems, MPAS-Analysis and other python-based software that used
NetCDF files based on the HDF5 file structure can experience file access errors
unless the following environment variable is set as follows in bash:

.. code-block:: bash

    $ export HDF5_USE_FILE_LOCKING=FALSE

or under csh:

.. code-block:: csh

    > setenv HDF5_USE_FILE_LOCKING FALSE

Then, running MPAS-Analysis is as simple as:

.. code-block:: bash

    $ mpas_analysis myrun.cfg

.. note::
   Adding ``--purge`` and/or ``--verbose`` flags to the ``mpas_analysis``
   command can be helpful, see the troubleshooting section below.

Typical output is the analysis is running correctly looks something like:

.. code-block:: none

    running: /home/xylar/Desktop/miniconda3/envs/mpas-analysis/bin/ESMF_RegridWeight
    Gen --source /tmp/tmph58_hgz4/src_mesh.nc --destination /tmp/tmph58_hgz4/dst_mes
    h.nc --weight /home/xylar/Desktop/analysis_test/analysis/A_WCYCL1850.ne4_oQU480.
    anvil/mapping/map_oQU480_to_0.5x0.5degree_bilinear.nc --method bilinear --netcdf
    4 --no_log --src_regional --ignore_unmapped
    running: /home/xylar/Desktop/miniconda3/envs/mpas-analysis/bin/ESMF_RegridWeight
    Gen --source /tmp/tmpxt8x1h_6/src_mesh.nc --destination /tmp/tmpxt8x1h_6/dst_mes
    h.nc --weight /home/xylar/Desktop/analysis_test/analysis/A_WCYCL1850.ne4_oQU480.
    anvil/mapping/map_obs_eke_0.25x0.25degree_to_0.5x0.5degree_bilinear.nc --method
    bilinear --netcdf4 --no_log --src_regional --ignore_unmapped
    running: /home/xylar/Desktop/miniconda3/envs/mpas-analysis/bin/ESMF_RegridWeight
    Gen --source /tmp/tmp3_7gpndz/src_mesh.nc --destination /tmp/tmp3_7gpndz/dst_mes
    h.nc --weight /home/xylar/Desktop/analysis_test/analysis/A_WCYCL1850.ne4_oQU480.
    anvil/mapping/map_oQU480_to_6000.0x6000.0km_10.0km_Antarctic_stereo_bilinear.nc
    --method bilinear --netcdf4 --no_log --src_regional --dst_regional --ignore_unma
    pped
    Preprocessing SOSE transect data...
      temperature
      salinity
      potentialDensity
      zonalVelocity
      meridionalVelocity
      velMag
      Done.
    running: /home/xylar/Desktop/miniconda3/envs/mpas-analysis/bin/ESMF_RegridWeight
    Gen --source /tmp/tmpt9n4vb5n/src_mesh.nc --destination /tmp/tmpt9n4vb5n/dst_mes
    h.nc --weight /home/xylar/Desktop/analysis_test/analysis/A_WCYCL1850.ne4_oQU480.
    anvil/mapping/map_oQU480_to_SOSE_transects_5km_bilinear.nc --method bilinear --n
    etcdf4 --no_log --src_regional --dst_regional --ignore_unmapped

    Running tasks:   2% |                                          | ETA:   0:09:04
    Running tasks:  52% |######################                    | ETA:   0:06:13
    Running tasks: 100% |##########################################| Time:  0:18:50

    Log files for executed tasks can be found in /home/xylar/Desktop/analysis_test/a
    nalysis/A_WCYCL1850.ne4_oQU480.anvil/logs
    Total setup time: 0:03:11.74
    Total run time: 0:22:02.33
    Generating webpage for viewing results...
    Done.

The first part of the output, before the progress bar, is the "setup" phase
where MPAS-Analysis is checking if the requested analysis can be run on the
simulation results.  The specific output shown here is related to creating
so-called mapping files that are used to interpolate between the ``oQU480`` mesh
and the various grids MPAS-Analysis uses to compare with observations.  Since
MPAS-Analysis didn't know about that ``oQU480`` mesh ahead of time, it is
creating mapping files and regions masks for this mesh on the fly.

The command-line tool has several more options you can explore with

.. code-block:: bash

    $ mpas_analysis --help

These include listing the available analysis tasks and their tags, purging a
previous analysis run before running the analysis again, plotting all available
color maps, and outputting verbose python error messages when the analysis fails
during the setup phase (before a progress bar appears).

6 Viewing the Output
--------------------

The primary output from MPAS-Analysis is a set of web pages, each containing
galleries of figures.  The output can be found in the directory you provided in
Section 4.5, which is the ``html`` subdirectory of the base output directory by
default.  If the web page is incomplete, it presumably means there was an error
during the analysis run, since the web page is generated as the final step.

The main web page has links to the ocean and sea-ice web pages as well as some
"provenance" information about which version of MPAS-Analysis you were using
and how it was configured.

The web page generated by this tutorial should look something like
`this example output`_.

.. _tutorial_getting_started_trouble:

7 Troubleshooting
-----------------

This section briefly describes strategies for diagnosing errors in
MPAS-Analysis.  This tutorial cannot hope to provide a comprehensive guide to
troubleshooting these errors.  Please search the documentation, Google the error
online, or get in touch with the MPAS-Analysis developer team (by
`posting an issue`_ on GitHub) if you are experiencing an error.


.. note::
   When all else fails, manually remove the entire MPAS-Analysis data output
   directories and the public html pages.

7.1 Purging old Analysis
~~~~~~~~~~~~~~~~~~~~~~~~

One thing you might want to try first if you are experiencing problems is to
delete any analysis you may already have in your output directory:

.. code-block:: bash

    $ mpas_analysis --purge myrun.cfg

This will first delete existing analysis and then run the analysis again.

7.2 Errors During Setup
~~~~~~~~~~~~~~~~~~~~~~~

If an error occurs during setup, by default the full python traceback is
suppressed.  This is because some tasks fail because the run being analyzed was
not configured for that analysis.  In such cases, many users want the analysis
to continue, simply skipping the tasks that can't be run.

However, this means that sometime the analysis is not configured properly and
as a results most or all tasks are not running.  To find out why, you will
probably need to run:

.. code-block:: bash

    $ mpas_analysis --verbose myrun.cfg

This will give you a detailed python stack trace.  Even if this is not helpful
to you, it might help developers to troubleshoot your issue.

7.2 Errors in Tasks
~~~~~~~~~~~~~~~~~~~

If you see the progress bar start but error occur during running of analysis
tasks, the error messages will not be displayed to the screen.  Instead, they
will be in log files (as stated in the short error message letting you know
that a task has failed).  The contents of these log files may help you to
determine the cause of the error.  If not, please include them if you are
`posting an issue`_ on GitHub.


.. _`Miniforge3 installer`: https://github.com/conda-forge/miniforge?tab=readme-ov-file#miniforge3
.. _`conda-forge channel`: https://conda-forge.org/
.. _`config file`: https://github.com/MPAS-Dev/MPAS-Analysis/tree/main/configs
.. _`Ultra-low-res ocean and sea-ice dataset`: https://web.lcrc.anl.gov/public/e3sm/diagnostics/test_output/20200305.A_WCYCL1850.ne4_oQU480.anvil/20200305.A_WCYCL1850.ne4_oQU480.anvil.ocean_seaice.tar.gz
.. _`full list of files`: https://github.com/MPAS-Dev/MPAS-Analysis/blob/main/mpas_analysis/obs/analysis_input_files
.. _`diagnostics web server`: https://web.lcrc.anl.gov/public/e3sm/diagnostics/
.. _`default config file`: https://github.com/MPAS-Dev/MPAS-Analysis/blob/main/mpas_analysis/default.cfg
.. _`example config file`: https://github.com/MPAS-Dev/MPAS-Analysis/blob/main/example.cfg
.. _`E3SM example config file`: https://github.com/MPAS-Dev/MPAS-Analysis/blob/main/example_e3sm.cfg
.. _`ncclimo`: http://nco.sourceforge.net/nco.html#ncclimo-netCDF-Climatology-Generator
.. _`this example output`: https://mpas-dev.github.io/MPAS-Analysis/examples/QU480
.. _`posting an issue`: https://github.com/MPAS-Dev/MPAS-Analysis/issues
