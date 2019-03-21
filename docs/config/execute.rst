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
  # Set this to "bck" (background parallelism) if running on a machine that can
  # handle 12 simultaneous processes, one for each monthly climatology.
  # Set to "mpi" to run one MPI task on each node and however many threads per
  # node to reach 12 total threads.
  ncclimoParallelMode = serial

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
used, 12 threads are spawned, one for each month, and then separate threads
are used to compute each season.  Given that computing climatologies takes up
a significant portion of the runtime in MPAS-Analysis, the speed-up of nearly
a factor of 12 in these computations can be quite noticeable.  For very big
data sets, it may be necessary to run ``ncclimo`` on multiple nodes to prevent
running out of memory.  Spawn a job with between 2 and 12 nodes, and set
``ncclimoParallelMode = mpi`` to run the 12 ``ncclimo`` threads on multiple
nodes.

Again, when running MPAS-Analysis on login nodes of supercomputing facilities,
it is important to be aware of the policies regarding using shared resources,
and ``bck`` may only be appropriate when running jobs on the compute nodes.
``mpi`` mode may not work at all on login nodes.

.. _`NetCDF Operators (NCO) package`: http://nco.sourceforge.net/nco.html