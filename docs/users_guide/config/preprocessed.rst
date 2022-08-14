.. _config_preprocessed:

Preprocessed Reference Runs
===========================

The ``[oceanPreprocessedReference]`` and ``[seaIcePreprocessedReference]``
sections of a configuration file contain options used to point to preprocessed
data from E3SM v0 reference runs::

  [oceanPreprocessedReference]
  ## options related to preprocessed ocean reference run with which the results
  ## will be compared (e.g. a POP, CESM or ACME v0 run)

  # directory where ocean reference simulation results are stored
  baseDirectory = /dir/to/ocean/reference

  ...

  [seaIcePreprocessedReference]
  ## options related to preprocessed sea ice reference run with which the results
  ## will be compared (e.g. a CICE, CESM or ACME v0 run)

  # directory where ocean reference simulation results are stored
  baseDirectory = /dir/to/seaice/reference

If such a preprocessed reference run is available, the name of the reference
run should be specified (see :ref:`config_runs`) and the base directories
should be specified here.
