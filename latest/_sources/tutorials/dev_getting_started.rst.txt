.. _tutorial_dev_getting_started:

Developer: Getting Started
==========================

This mini-tutorial is meant as the starting point for other tutorials for
developers.  It describes the process for creating a fork of the MPAS-Analysis
repo, cloning the repository (and your fork) locally, making a git worktree for
development, and creating a conda environment that includes the
``mpas_analysis`` package and all of its dependencies, installed in a mode
appropriate for development.

1. Getting started on GitHub
----------------------------

1.1 Forking MPAS-Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

If you would like to contribute to MPAS-Analysis, you will need to create your
own fork of the `repository <https://github.com/MPAS-Dev/MPAS-Analysis>`_.  Go
to the link and click on ``Fork`` near the top right corner of the page.  The
Owner should be your GitHub username and the Repository name should be
``MPAS-Analysis``.  Check the box for "Copy the develop branch only".  Click
"Create fork".

1.2 Adding SSH keys
~~~~~~~~~~~~~~~~~~~

If you have not already done so, you should add SSH keys to GitHub that allow
you to push to your fork from the machine(s) where you will do your
development.  Instructions can be found
`here <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_.

1.3 Local git configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~

It will be convenient to have some basic configuration for ``git`` taken care
of before we clone the repository.  Here are some recommended config options
to set.  Edit your ``~/.gitconfig`` (create it if it doesn't exist).

.. code-block:: ini

    [user]
            name = Xylar Asay-Davis
            email = xylarstorm@gmail.com
    [core]
            editor = vim
    [color]
            ui = true
    [alias]
            logg = log --graph --oneline --decorate
    [rebase]
            autosquash = true

Obviously, change ``[user]`` config options to appropriate values for you.
You *must* use the email address associated with your GitHub account.
Otherwise, your commits will not be associated with your GitHub user name.

2. Cloning the repository
-------------------------

You will want to clone both the main MPAS-Analysis repository and your own
fork.  The MPAS-Analysis development tutorials assume that you will be
developing branches in different worktrees and recommend a directory structure
appropriate for this approach.

.. note::

    If you are on a machine with an old version of ``git``, you may need to
    add:

    .. code-block:: bash

        module load git

    to your ``.bashrc``.  You want a pretty recent version of ``git`` so you
    have the ``git worktree`` command.

Begin by creating a "base" directory for development in a convenient location
for keeping code.  This should not be on a "scratch" or other temporary drive
on an HPC machine.  The base directory should be named ``MPAS-Analysis``,
``mpas-analysis`` or something similar.

.. code-block:: bash

    $ mkdir mpas-analysis
    $ cd mpas-analysis

Within the base directory, clone the main repository into a directory called
``develop`` (the default branch is the ``develop`` branch):

.. code-block:: bash

    $ git clone git@github.com:MPAS-Dev/MPAS-Analysis.git develop
    $ cd develop

Add your fork as a "remote":

.. code-block:: bash

    $ git remote add <username>/MPAS-Analysis git@github.com:<username>/MPAS-Analysis.git

Make sure to replace ``<username>`` with your GitHub username.


3. Making a worktree
--------------------

To do your development, first make sure you are in the ``develop`` directory
within your base directory (e.g. ``mpas-analysis/develop``).  Then, "fetch" and
changes that might have happened on the ``develop`` branch so you are using
the latest version as a starting point:

.. code-block:: bash

    $ git fetch --all -p

This will fetch all branches from both the main repository and your fork.  It
will also prune (``-p``) any branches you might have deleted.

Then, make a worktree for developing your new feature:

.. code-block:: bash

    $ git worktree add ../add_my_fancy_task

The last argument (``add_my_fancy_task`` in this example) is both the name of
a directory within the base directory (``mpas-analysis``) and the name of the
branch you will be developing.

Go into that directory to do your development:

.. code-block:: bash

    $ cd ../add_my_fancy_task

4. Making a conda environment
-----------------------------

MPAS-Analysis relies on several packages that are only available as conda
packages from the ``conda-forge`` channel.  The first step for running
MPAS-Analysis is to create a conda environment with all the needed packages.

4.1 Installing Mambaforge
~~~~~~~~~~~~~~~~~~~~~~~~~

If you have not yet installed Anaconda, Miniconda or Mambaforge, you will need
to begin there.  The concept behind Anaconda is that just about everything you
would need for a typical python workflow is included.  The concept behind
Miniconda and Mambaforge is that you create different environments for
different purposes.  This allows for greater flexibility and tends to lead to
fewer conflicts between incompatible packages, particularly when using a
channel other than the ``defaults`` supplied by Anaconda.  Since we will use
the ``conda-forge`` channel and the ``mamba`` tools to speed up installation,
the Mambaforge approach is strongly recommended.  The main advantage of
Mambaforge over Miniconda is that it automatically takes care of a few steps
that we otherwise need to do manually.

First download the
`Mambaforge installer <https://github.com/conda-forge/miniforge#mambaforge>`_
for your operating system, then run it:

.. code-block:: bash

   $ /bin/bash Mambaforge-Linux-x86_64.sh

.. note::

   MPAS-Analysis and many of the packages it depends on support OSX and Linux
   but not Windows.

If you are on an HPC system, you can still install Miniconda into your home
directory.  Typically, you will need the Linux version.

.. note::

    At this time, we don't have experience with installing or running
    MPAS-Analysis on ARM or Power8/9 architectures.

You will be asked to agree to the terms and conditions. Type ``yes`` to
continue.

You will be prompted with a location to install. In this tutorial, we assume
that Mambaforge is installed in the default location, ``~/mambaforge``.  If
you are using Miniconda or chose to install Mambaforge somewhere else, just
make sure to make the appropriate substitution whenever you see a reference to
this path below.

.. note::

    On some HPC machines (particularly at LANL Institutional Computing and
    NERSC) the space in your home directory is quite limited.  You may want to
    install Mambaforge in an alternative location to avoid running out of
    space.

You will see prompt like this:

.. code-block::

    Do you wish the installer to initialize Mambaforge
    by running conda init? [yes|no]
    [no] >>>

You may wish to skip the step (answer ``no``) if you are working on a system
where you will also be using other conda environments, most notably
E3SM-Unified (which has its own Miniconda installation).  If you do not run
conda init, you have to manually activate ``conda`` whenever you need it.
For ``bash`` and similar shells, this is:

.. code-block:: bash

   $ source ~/mambaforge/etc/profile.d/conda.sh
   $ conda activate

If you use ``csh``, ``tcsh`` or related shells, this becomes:

.. code-block:: csh

   > source ~/mambaforge/etc/profile.d/conda.csh
   > conda activate

You may wish to create an alias in your ``.bashrc`` or ``.cshrc`` to make
this easier.  For example:

.. code-block:: bash

   alias init_conda="source ~/mambaforge/etc/profile.d/conda.sh; conda activate"


4.2 One-time Miniconda setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you installed Miniconda, rather than Mambaforge, you will need to add the
`conda-forge channel <https://conda-forge.org/>`_ and make sure it always takes
precedence for packages available on that channel:

.. code-block:: bash

   $ conda config --add channels conda-forge
   $ conda config --set channel_priority strict

Then, you will need to install the ``mamba`` package:

.. code-block:: bash

   $ conda install -y mamba

If you installed Mambaforge, these steps will happen automatically.

4.3 Create a development environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can create a new conda environment called ``mpas_dev`` and install the
dependencies that MPAS-Analysis needs by running the following in the worktree
where you are doing your development:

.. code-block:: bash

   $ mamba create -y -n mpas_dev --file dev-spec.txt "esmf=*=nompi_*"

The last argument is only needed on HPC machines because the conda version of
MPI doesn't work properly on these machines.  You can omit it if you're
setting up the conda environment on your laptop.

Then, you can activate the environment and install MPAS-Analysis in "edit"
mode by running:

.. code-block:: bash

   $ conda activate mpas_dev
   $ python -m pip install -e .

In this mode, any edits you make to the code in the worktree will be available
in the conda environment.  If you run ``mpas_analysis`` on the command line,
it will know about the changes.

.. note::

    If you add or remove files in the code, you will need to re-install
    MPAS-Analysis in the conda environment by rerunning

    .. code-block:: bash

       python -m pip install -e .

.. _tutorial_dev_get_started_activ_env:

4.4 Activating the environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each time you open a new terminal window, to activate the ``mpas_dev``
environment, you will need to run either for ``bash``:

.. code-block:: bash

   $ source ~/mambaforge/etc/profile.d/conda.sh
   $ conda activate mpas_dev

or for ``csh``:

.. code-block:: csh

   > source ~/mambaforge/etc/profile.d/conda.csh
   > conda activate mpas_dev

You can skip the ``source`` command if you chose to initialize Mambaforge or
Miniconda3 so it loads automatically.  You can also use the ``init_conda``
alias for this step if you defined one.

4.5 Switching worktrees
~~~~~~~~~~~~~~~~~~~~~~~

If you switch to a different worktree, it is safest to rerun the whole
process for creating the ``mpas_dev`` conda environment.  If you know that
the dependencies are the same as the worktree used to create ``mpas_dev``,
You can just reinstall ``mpas_analysis`` itself by rerunning

.. code-block:: bash

    python -m pip install -e .

in the new worktree.  If you forget this step, you will find that changes you
make in the worktree don't affect the ``mpas_dev`` conda environment you are
using.

5. Editing code
---------------

You may, of course, edit the MPAS-Analysis code using whatever tool you like.
I strongly recommend editing on your laptop and using
`PyCharm community edition <https://www.jetbrains.com/pycharm/download/>`_
to do the editing.  PyCharm provides many features including flagging
deviations from preferred coding style guidelines known as
`PEP8 <https://peps.python.org/pep-0008/>`_ and syntax error detection using
the ``mpas_dev`` conda environment you created.

6. Running MPAS-Analysis on a laptop
------------------------------------

If you wish to run MPAS-Analysis on your laptop (or desktop machine), you will
need to follow steps 2-6 of the :ref:`tutorial_getting_started` tutorial.

7. Running MPAS-Analysis on an E3SM supported machine
-----------------------------------------------------

7.1 Configuring MPAS-Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We configure MPAS-Analysis is with Python ``cfg`` (also called ``ini``) files:

.. code-block:: ini

   [runs]
   # mainRunName is a name that identifies the simulation being analyzed.
   mainRunName = runName

   [execute]
   ...

The `default config file <https://github.com/MPAS-Dev/MPAS-Analysis/blob/main/mpas_analysis/default.cfg>`_
contains thousands of config options, which gives a lot of flexibility to
MPAS-Analysis but can be more than bit overwhelming to new users and
developers.

The file `example_e3sm.cfg <https://github.com/MPAS-Dev/MPAS-Analysis/blob/develop/example_e3sm.cfg>`_
provides you with an example with some of the most common config options you
might need to change on an E3SM supported machine.  If you specify the name of
the supported machine with the ``--machine`` (or ``-m``) flag when you call
``mpas_analysis``, there are several config options that will be set for you
automatically.

First, you should copy this file to a new name for a specific run (say
``myrun.cfg``).  Then, you should modify any config options you want to change
in your new config file. At a minimum, you need to specify:

* ``mainRunName`` in ``[runs]``:  A name for the run to be included plot titles
  and legends (best if it's not super long)
* ``baseDirectory`` in ``[input]``: The directory for the simulation results
  to analyze
* ``mpasMeshName`` in ``[input]``: The official name of the MPAS-Ocean and
  -Seaice mesh
* ``baseDirectory`` in ``[output]``: The directory for the analysis results

We will cover these and a few other common options in this tutorial.  With the
exception of a few paths that you will need to provide, the config options
displayed below are the ones appropriate for the example E3SM simulation using
the QU480 MPAS mesh.

7.1.1 [runs]
++++++++++++

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

7.1.2 [execute]
+++++++++++++++

The ``[execute]`` section contains options related to serial or parallel
execution of the individual "tasks" that make up an MPAS-Analysis run.  For
the most part, you can let MPAS-Analysis take care of this on supported
machines.  The exception is that, in a development conda environment, you will
be using a version of ESMF that cannot run in parallel so you will need the
following:

.. code-block:: ini

    [execute]
    ## options related to executing parallel tasks

    # the number of MPI tasks to use in creating mapping files (1 means tasks run in
    # serial, the default)
    mapMpiTasks = 1

    # "None" if ESMF should perform mapping file generation in serial without a
    # command, or one of "srun" or "mpirun" if it should be run in parallel (or in
    # serial but with a command)
    mapParallelExec = None

If you are running into trouble with MPAS-Analysis, such as running out of
memory, you may want to explore other config options from this section.

7.1.3 [input]
+++++++++++++

The ``[input]`` section provides paths to the E3SM simulation data and the name
of the MPAS-Ocean and MPAS-Seaice mesh.

.. code-block:: ini

    [input]
    ## options related to reading in the results to be analyzed

    # directory containing model results
    baseDirectory = /lcrc/group/e3sm/ac.xylar/acme_scratch/anvil/20200305.A_WCYCL1850.ne4_oQU480.anvil

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

The ``baseDirectory`` is the path for the E3SM simulation. Here are paths to
some very low resolution simulations you can use on various supported machines:

Anvil or Chrysalis:

.. code-block::

    /lcrc/group/e3sm/ac.xylar/acme_scratch/anvil/20200305.A_WCYCL1850.ne4_oQU480.anvil
    /lcrc/group/e3sm/ac.xylar/acme_scratch/anvil/20201025.GMPAS-IAF.T62_oQU240wLI.anvil

Cori and Perlmutter:

.. code-block::

    /global/cfs/cdirs/e3sm/xylar/20200305.A_WCYCL1850.ne4_oQU480.anvil

Compy:

.. code-block::

    /compyfs/asay932/analysis_testing/test_output/20200305.A_WCYCL1850.ne4_oQU480.anvil


The ``mpasMeshName`` is the official name of the MPAS-Ocean and -Seaice mesh
used in the simulation, which should be in the simulation name and must be a
directory on the
`inputdata <https://web.lcrc.anl.gov/public/e3sm/inputdata/ocn/mpas-o/>`_
server  In this example, this is ``oQU480``, meaning the quasi-uniform 480-km
mesh for the ocean and sea ice.

The ``runSubdirectory`` must contain valid MPAS-Ocean and MPAS-Seaice restart
files, used to get information about the MPAS mesh and the ocean vertical grid.

The ``oceanHistorySubdirectory`` must contain MPAS-Ocean monthly mean output
files, typically named::

   mpaso.hist.am.timeSeriesStatsMonthly.YYYY-MM-DD.nc

Similarly, ``seaIceHistorySubdirectory`` contains the MPAS-Seaice monthly mean
output::

   mpassi.hist.am.timeSeriesStatsMonthly.YYYY-MM-DD.nc

In this example, we are using a run where short-term archiving has been used
so the output is not in the ``run`` directory.

Finally, MPAS-Analysis needs a set of "namelists" and "streams" files that
provide information on the E3SM configuration for MPAS-Ocean and MPAS-Seaice,
and about the output files, respectively.  These are typically also found in
the ``run`` directory.

.. _tutorial_dev_get_started_config_output:

7.1.4 [output]
++++++++++++++

The ``[output]`` section provides a path where the output from the analysis run
will be written, the option to output the results web pages to another
location, and a list of analysis to be generated (or explicitly skipped).

.. code-block:: ini

    [output]
    ## options related to writing out plots, intermediate cached data sets, logs,
    ## etc.

    # The subdirectory for the analysis and output on the web portal
    subdir = ${runs:mainRunName}/clim_${climatology:startYear}-${climatology:endYear}_ts_${timeSeries:startYear}-${timeSeries:endYear}

    # directory where analysis should be written
    # NOTE: This directory path must be specific to each test case.
    baseDirectory = /lcrc/group/e3sm/${web_portal:username}/analysis/${output:subdir}

    # provide an absolute path to put HTML in an alternative location (e.g. a web
    # portal)
    htmlSubdirectory = ${web_portal:base_path}/${web_portal:username}/analysis/${output:subdir}

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
    generate = ['all', 'no_BGC', 'no_icebergs', 'no_index', 'no_eke',
                'no_landIceCavities']

In this example, I have made liberal use of
`extended interpolation <https://docs.python.org/3/library/configparser.html#configparser.ExtendedInterpolation>`_
in the config file to make use of config options in other config options.

``subdir`` is the subdirectory for both the analysis and the output on the
web portal.  It typically indicates the run being used and the years covered
by the climatology (and sometimes the time series as in this example).  See
:ref:`tutorial_dev_get_started_config_clim` for more info on these config
options.

``baseDirectory`` is any convenient location for the output.  In this example,
I have used a typical path on Anvil or Chrysalis, including the
``${web_portal:username}`` that will be populated automatically on a supported
machine and ``${output:subdir}``, the subdirectory from above.

``htmlSubdirectory`` is set using the location of the web portal, which is
automatically determined on an E3SM machine, the user name, and the same
subdirectory used for analysis output.  You can modify the path as needed to
match your own preferred workflow.

.. note::

    On some E3SM supported machines like Chicoma, there is no
    web portal so you will want to just manually replace the part of the
    ``basePath`` given by ``/lcrc/group/e3sm/${web_portal:username}`` in the
    example above.

    You will need to just put the web output in an ``html`` subdirectory within
    the analysis output:

    .. code-block:: ini

        htmlSubdirectory = html

    and copy this from the supercomputer to your laptop to view it in your
    browser.

Finally, the ``generate`` option provides a python list of flags that can be
used to determine which analysis will be generated.  In this case, we are
turning off some analysis that will not work because some features
(biogeochemistry, icebergs, eddy kinetic energy and land-ice cavities) are not
available in this run and one (the El Niño climate index) is not useful.

.. _tutorial_dev_get_started_config_clim:

7.1.5. [climatology], [timeSeries] and [index]
++++++++++++++++++++++++++++++++++++++++++++++

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

    # start and end years for timeseries analysis. Use endYear = end to indicate
    # that the full range of the data should be used.  If errorOnMissing = False,
    # the start and end year will be clipped to the valid range.  Otherwise, out
    # of bounds values will lead to an error.  In a "control" config file used in
    # a "main vs. control" analysis run, the range of years must be valid and
    # cannot include "end" because the original data may not be available.
    startYear = 1
    endYear = 5

    [index]
    ## options related to producing nino index.

    # start and end years for El Nino 3.4 analysis. Use endYear = end to indicate
    # that the full range of the data should be used.  If errorOnMissing = False,
    # the start and end year will be clipped to the valid range.  Otherwise, out
    # of bounds values will lead to an error.  In a "control" config file used in
    # a "main vs. control" analysis run, the range of years must be valid and
    # cannot include "end" because the original data may not be available.
    startYear = 1
    endYear = 5

For each of these, options a full year of data must exist for that year to
be included in the analysis.

For the example QU480 simulation, only 5 years of output are available, so we
are doing a climatology over the last 3 years (3 to 5) and displaying time
series over the full 5 years.  (If the El Niño index weren't disabled, it would
also be displayed over the full 5 years.)

7.2 Running MPAS-Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

The hard work is done.  Now that we have a config file, we are ready to run.

To run MPAS-Analysis, you should either create a job script or log into
an interactive session on a compute node.  Then, activate the ``mpas_dev``
conda environment as in :ref:`tutorial_dev_get_started_activ_env`.

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

    $ mpas_analysis -m <machine> myrun.cfg

where ``<machine>`` is the name of the machine (all lowercase).  On Cori, we
only support the Haswell nodes (so the machine name is ``cori-haswell``).  For
now, we only support CPU nodes on Perlmutter (``pm-cpu``) and Chicoma
(``chicoma-cpu``).

Typical output is the analysis is running correctly looks something like:

.. code-block:: none

    $ mpas_analysis -m anvil myrun.cfg
    Detected E3SM supported machine: anvil
    Using the following config files:
       /gpfs/fs1/home/ac.xylar/code/mpas-analysis/add_my_fancy_task/mpas_analysis/default.cfg
       /gpfs/fs1/home/ac.xylar/anvil/mambaforge/envs/mpas_dev/lib/python3.10/site-packages/mache/machines/anvil.cfg
       /gpfs/fs1/home/ac.xylar/code/mpas-analysis/add_my_fancy_task/mpas_analysis/configuration/anvil.cfg
       /gpfs/fs1/home/ac.xylar/code/mpas-analysis/add_my_fancy_task/mpas_analysis/__main__.py
       /gpfs/fs1/home/ac.xylar/code/mpas-analysis/add_my_fancy_task/myrun.cfg
    copying /gpfs/fs1/home/ac.xylar/code/mpas-analysis/add_my_fancy_task/myrun.cfg to HTML dir.

    running: /gpfs/fs1/home/ac.xylar/anvil/mambaforge/envs/mpas_dev/bin/ESMF_RegridWeightGen --source /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/tmp76l7of28/src_mesh.nc --destination /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/tmp76l7of28/dst_mesh.nc --weight /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/map_oQU480_to_0.5x0.5degree_bilinear.nc --method bilinear --netcdf4 --no_log --src_loc center --src_regional --ignore_unmapped
    running: /gpfs/fs1/home/ac.xylar/anvil/mambaforge/envs/mpas_dev/bin/ESMF_RegridWeightGen --source /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/tmpj94wpf9y/src_mesh.nc --destination /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/tmpj94wpf9y/dst_mesh.nc --weight /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/map_oQU480_to_6000.0x6000.0km_10.0km_Antarctic_stereo_bilinear.nc --method bilinear --netcdf4 --no_log --src_loc center --src_regional --dst_regional --ignore_unmapped
    running: /gpfs/fs1/home/ac.xylar/anvil/mambaforge/envs/mpas_dev/bin/ESMF_RegridWeightGen --source /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/tmp6zm13a0s/src_mesh.nc --destination /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/tmp6zm13a0s/dst_mesh.nc --weight /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/map_oQU480_to_WOCE_transects_5km_bilinear.nc --method bilinear --netcdf4 --no_log --src_loc center --src_regional --dst_regional --ignore_unmapped
    Preprocessing SOSE transect data...
      temperature
      salinity
      potentialDensity
      zonalVelocity
      meridionalVelocity
      velMag
      Done.
    running: /gpfs/fs1/home/ac.xylar/anvil/mambaforge/envs/mpas_dev/bin/ESMF_RegridWeightGen --source /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/tmpe2a9yblb/src_mesh.nc --destination /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/tmpe2a9yblb/dst_mesh.nc --weight /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/mapping/map_oQU480_to_SOSE_transects_5km_bilinear.nc --method bilinear --netcdf4 --no_log --src_loc center --src_regional --dst_regional --ignore_unmapped

    Running tasks: 100% |##########################################| Time:  0:06:42

    Log files for executed tasks can be found in /lcrc/group/e3sm/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5/logs
    Total setup time: 0:02:13.78
    Total run time: 0:08:55.86
    Generating webpage for viewing results...
    Web page: https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5


The first part of the output, before the progress bar, is the "setup" phase
where MPAS-Analysis is checking if the requested analysis can be run on the
simulation results.  The specific output shown here is related to 1)
listing the config files used to determine the final set of config options
used in the analysis, and 2) creating mapping files that are used to
interpolate between the ``oQU480`` mesh and the various grids MPAS-Analysis
uses to compare with observations. Since MPAS-Analysis didn't know about that
``oQU480`` mesh ahead of time, it is creating mapping files and regions masks
for this mesh on the fly.

The ``mpas_analysis`` command-line tool has several more options you can
explore with

.. code-block:: bash

    $ mpas_analysis --help

These include listing the available analysis tasks and their tags, purging a
previous analysis run before running the analysis again, plotting all available
color maps, and outputting verbose python error messages when the analysis
fails during the setup phase (before a progress bar appears).

7.3 Viewing the Output
~~~~~~~~~~~~~~~~~~~~~~

The primary output from MPAS-Analysis is a set of web pages, each containing
galleries of figures.  The output can be found in the directory you provided in
:ref:`tutorial_dev_get_started_config_output` and given in the last line of
the analysis output (if you are on a supported machine with a web portal),
e.g.:

.. code-block:: none

    Web page: https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.xylar/analysis/A_WCYCL1850.ne4_oQU480.anvil/clim_3-5_ts_1-5

.. note::

    On Cori and Perlmutter, you will need to change the permissions so you can
    see the webpage online:

    .. code-block:: bash

        $ chmod -R ugo+rX /global/cfs/cdirs/e3sm/www/<username>

    where ``<username>`` is your NERSC username.

If the web page is incomplete, it presumably means there was an error during
the analysis run, since the web page is generated as the final step.  Check
the analysis output and then the log files for individual analysis tasks to
see what when wrong.  See :ref:`tutorial_getting_started_trouble` or ask for
help if you run into trouble.

The main web page has links to the ocean and sea-ice web pages as well as some
"provenance" information about which version of MPAS-Analysis you were using
and how it was configured.

The web page generated by this tutorial should look something like this
(somewhat outdated)
`example output <https://mpas-dev.github.io/MPAS-Analysis/examples/QU480>`_.
