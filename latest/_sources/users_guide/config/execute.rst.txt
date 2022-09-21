.. _config_execute:

Execute
=======

The ``[execute]`` section of a configuration file contains options used to
control how tasks are executed within an MPAS-Analysis run::

    [execute]
    ## options related to executing parallel tasks

    # the number of parallel tasks (1 means tasks run in serial, the default)
    parallelTaskCount = 1

    # the parallelism mode in ncclimo ("serial", "bck" or "mpi")
    # Set this to "bck" (background parallelism) in most cases.  The default number
    # of threads (see below) is 12, one for each monthly climatology. Set to "mpi"
    # to run one MPI task on each node and however many threads per node to reach
    # 12 total threads.
    ncclimoParallelMode = serial

    # the number of total threads to use when ncclimo runs in "bck" or "mpi" mode.
    # Reduce this number if ncclimo is crashing (maybe because it is out of memory).
    # The number of threads must be a factor of 12 (1, 2, 3, 4, 6 or 12).
    ncclimoThreads = 12

    # the number of MPI tasks to use in creating mapping files (1 means tasks run in
    # serial, the default)
    mapMpiTasks = 1

    # "None" if ESMF should perform mapping file generation in serial without a
    # command, or one of "srun" or "mpirun" if it should be run in parallel (or ins
    # serial but with a command)
    mapParallelExec = None

    # "None" if ncremap should perform remapping without a command, or "srun"
    # possibly with some flags if it should be run with that command
    ncremapParallelExec = None

    # Multiprocessing method used in python mask creation ("forkserver", "fork" or
    # Multiprocessing method used in python mask creation ("forkserver", "fork" or
    # "spawn").  We have found that "spawn" is the only one that works in python
    # 3.7 on Anvil so this is the default
    multiprocessingMethod = spawn

Parallel Tasks
--------------

By default, MPAS-Analysis runs one task at a time, displaying any logging
output directly to the screen, rather than storing it in log files.  However,
the analysis will run much more quickly if parallel processes are used. For
example, to run 12 tasks in parallel, simply set::

  parallelTaskCount = 12

MPAS-Analysis currently supports parallelism through spawning processes, rather
than with MPI, meaning that MPAS-Analysis should be run on a single node of a
multi-node machine such as a cluster or supercomputer.  Given that some tasks
themselves spawn multiple threads and that some tasks are memory intensive, it
may not be desirable to launch one task per core on a node with limited memory.

Because MPAS-Analysis does not use MPI parallelism, it can typically be run on
the login nodes of supercomputing facilities.  Check with the policies of your
center to see if this is permitted and make sure not to run with a large number
of parallel tasks so as to overwhelm the shared resource.

Parallelism in NCO
------------------

The ``ncclimo`` command from the `NetCDF Operators (NCO) package`_ is used
internally in MPAS-Analysis. This command supports three options for
parallelism: ``serial``, ``bck`` or ``mpi``.  If set to ``serial``, the
default, any MPAS-Analysis tasks that use ``ncclimo`` will compute
climatologies one month and then one season at a time.  If ``bck`` mode is
used, ``ncclimoThreads`` threads are spawned (default is 12, one for each
month), and then separate threads are used to compute each season.  Given that
computing climatologies takes up a significant portion of the runtime in
MPAS-Analysis, the speed-up of nearly a factor of ``ncclimoThreads`` in these
computations can be quite noticeable.  For very big data sets, it may be
necessary to run ``ncclimo`` either with fewer threads (reducing
``ncclimoThreads``, noting that it must be a factor of 12) or on multiple nodes
to prevent running out of memory.  To run an MPI job, spawn a job with between
2 and 12 nodes, and set ``ncclimoParallelMode = mpi`` to run the 12 ``ncclimo``
threads on multiple nodes.

Again, when running MPAS-Analysis on login nodes of supercomputing facilities,
it is important to be aware of the policies regarding using shared resources.
On login nodes, ``bck`` may only be appropriate with ``ncclimoThreads`` set to a
small number and ``mpi`` mode may not work at all.

Parallel Mapping File Creation
------------------------------

If mapping files from the MPAS mesh to the comparison grids aren't already
available in the diagnostics directory, they will be created before any other
MPAS-Analysis tasks are run.  If you are running MPAS-Analysis out of
E3SM-Unified on a compute node, on many systems (see below), ESMF has been
built with the system version of MPI and you must run mapping-file generation
with ``srun``.  If you are running with ``parallelTaskCount > 1``, the mapping
files will be generated in parallel.

.. code-block:: cfg

    mapParallelExec = srun

Similarly, some systems (Anvil and Chrysalis) require a parallel executable
for calls to ``ncremap`` from E3SM-Unified on compute nodes:

.. code-block:: cfg

    ncremapParallelExec = srun

E3SM supported machines with system MPI support in E3SM-Unified 1.5.0:

* Anvil

* Badger

* Chrysalis

* Compy

* Cori-Haswell

* Grizzly

These machines do **not** have MPI support in E3SM-Unified:

* Cori-KNL

* Cooley

* Andes

* Acme1

In the very near future, we hope to add a capability to MPAS-Analysis so that
it will automatically recognize which machine it is on (or you can specify if
need be), allowing these and other config options to be set automatically.

Parallel Mask Creation
----------------------

Tasks that involve :ref:`config_region_groups` can generate the masks for each
region in the group on the fly.  This is done with the mask generation
command-line tools form MPAS-Tools (see 
`Mask Creation with Python Multiprocessing <http://mpas-dev.github.io/MPAS-Tools/stable/mesh_conversion.html#mask-creation-with-pthon-multiprocessing>`_),
which support 3 modes of parallelism: "spawn", "fork" and "forkserver". For
technical details on these modes, see
`Contexts and start methods <https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods>`_.
We have found that "spawn" seems to be the most reliable option on Anvil under
python 3.7 and 3.8.  Any of these methods works well under python 3.8 but only
"spawn" was reliable under python 3.7.  Therefore, we use "spawn" as the
default.

As we gain more experience with this setting, we may update config files for
specific machines to have different defaults.



.. _`NetCDF Operators (NCO) package`: http://nco.sourceforge.net/nco.html
