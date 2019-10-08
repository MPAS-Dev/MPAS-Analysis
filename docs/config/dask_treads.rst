.. _dask_treads:

Dask treads and subprocess count
================================

Several tasks and subtasks have config options ``daskThreads`` and
``subprocessCount`` used to control threading within a subtask::

  # The number of threads dask is allowed to spawn for each task/subtask.
  # Decrease this number if tasks/subtasks are running out of available threads
  daskThreads = 2

  # The number of subprocesses that each task/subtask gets counted as
  # occupying. Increase this number if tasks/subtasks are running out of
  # memory, so that fewer tasks will be allowed to run at once
  subprocessCount = 1

Dask treads
-----------

Dask and xarray support thread-parallel operations on data sets.  They also
support chunk-wise operation on data sets that can't fit in memory.  These
capabilities are very powerful but also difficult to configure for general
cases.  Dask is also not desigend by default with the idea that multiple tasks,
each with multiple dask treads, might operate simultaneously.  As a result,
it is possible to spawn huge numbers of dask treads in MPAS-Analysis that both
slow down analysis and lead to errors when the node runs out of threads
completely.

To prevent this, many tasks or subtasks that use dask treading take the number
of execution threads from a config option, typically in the config section for
the parent task.  Typically, the number of ``daskThreads`` should be around
the same as the number of cores on a node divided by the number of tasks
that will run simultaneiously.  Since the number of running tasks is controlled
by ``subprocessCount``, see below, this number might differ from
``parallelTaskCount``.

Subprocess count
----------------

Tasks or subtasks that use dask treading may consume too much memory or use
too many threads to "count" as a single task.  That is, it might not be safe to
run with ``parallelTaskCount`` simultaneious instances of the task/subtask and
it would be better if it occupied the slot of multiple tasks in the pool of
tasks.  MPAS-Analysis will treat a dask-based task or subtask as occupying
the number of task slots given by the ``subprocessCount`` option.  For example,
if ``parallelTaskCount = 8`` and ``subprocessCount = 2``, up to 4 tasks or
subtasks would be allowed to run simultaneiously.
