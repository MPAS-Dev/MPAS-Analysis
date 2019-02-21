.. _config_input:

Input
=====

The ``[input]`` section of a configuration file contains options used to
specify the location of the "main" simulation and some settings for how
these data will be read in::

  [input]
  ## options related to reading in the results to be analyzed

  # directory containing model results
  baseDirectory = /dir/for/model/output

  # Note: an absolute path can be supplied for any of these subdirectories.
  # A relative path is assumed to be relative to baseDirectory.
  # By default, results are assumed to be directly in baseDirectory,
  # i.e. <baseDirecory>/./

  # subdirectory containing restart files
  runSubdirectory = .
  # subdirectory for ocean history files
  oceanHistorySubdirectory = .
  # subdirectory for sea ice history files
  seaIceHistorySubdirectory = .

  # names of namelist and streams files, either a path relative to baseDirectory
  # or an absolute path.
  oceanNamelistFileName = mpaso_in
  oceanStreamsFileName = streams.ocean
  seaIceNamelistFileName = mpassi_in
  seaIceStreamsFileName = streams.seaice

  # names of ocean and sea ice meshes (e.g. oEC60to30, oQU240, oRRS30to10, etc.)
  mpasMeshName = mesh

  # Large datasets can encounter a memory error.  Specification of a maximum
  # chunk size `maxChunkSize` can be helpful to prevent the memory error.  The
  # current maximum chunk size assumes approximately 64GB of ram and large files
  # with a single time slice.
  maxChunkSize = 10000

  # Directory for mapping files (if they have been generated already). If mapping
  # files needed by the analysis are not found here, they will be generated and
  # placed in the output mappingSubdirectory
  # mappingDirectory = /dir/for/mapping/files

Input Directories
-----------------

The input directories are specified through a base directory and optionally
separate subdirectories for the run (continaing namelist, streams and and one
or more restart files for each MPAS component) and for each component's
simulation output ("history" files).  You will always need to set
``baseDirectory`` in a custom config file.

By default, all files are assumed to be located in the base directory.
However, E3SM supports short-term archiving of component output to separate
directories.  If this feature was enabled for the E3SM run to be analyzed, the
configuration should be::

  baseDirectory = /dir/for/model/output
  runSubdirectory = run
  oceanHistorySubdirectory = archive/ocn/hist
  seaIceHistorySubdirectory = archive/ice/hist

where ``/dir/for/model/output`` is replaced with the directory containing
the ``run`` and ``archive`` directories.

Namelist Files
--------------

MPAS components are configured with a namelist file containing a very large
number of parameters and other configuration settings.  Part of the strength
of MPAS-Analysis is that it is aware of these namelists and can automatically
disable analysis tasks that are not supported under a given configuration.
By default, the namelist files for ocean and sea ice components are::

  oceanNamelistFileName = mpaso_in
  seaIceNamelistFileName = mpassi_in

For older E3SM (v1 alpha and beta) runs, a different naming convention was used
and these options will need to be updated to::

  oceanNamelistFileName = mpas-o_in
  seaIceNamelistFileName = mpas-cice_in

The location of the namelist files is relative to the baseDirectory, so if
they are located within a run directory inside the base directory, they should
instead be specified as::

  oceanNamelistFileName = run/mpaso_in
  seaIceNamelistFileName = run/mpassi_in

Streams Files
-------------

Streams files are used to configure input and output from MPAS components.
MPAS-Analysis can parse these files to determine the locations of output files
(using the ``*HistorySubdirectory`` to find them if they have been moved by
short-term archiving).  Similarly to namelist files, there are default names::

  oceanStreamsFileName = streams.ocean
  seaIceStreamsFileName = streams.seaice

alterations appropriate for E3SM v1 alpha and beta runs::


  oceanStreamsFileName = streams.ocean
  seaIceStreamsFileName = streams.cice

and the addition of the ``run/`` subdirectory if analyzing a run that used
short-term archiving::

  oceanStreamsFileName = run/streams.ocean
  seaIceStreamsFileName = run/streams.seaice

Mesh Name
---------

The MPAS-Ocean and MPAS-Seaice run on the same mesh.  There are a number of
standard E3SM ocean/sea ice meshes at various resolutions.  The meshes
currently supported by the public release of MPAS-Analysis include:

  * ``oEC60to30v3``: An Eddy-Closure (EC) mesh with 30-km resolution at the
    poles and equator and 60-km resolution at mid latitudes,
  * ``oRRS30to10v3``: A Rossby-Radius-Scaled (RRS) mesh with 10-km resolution
    at the poles and 30-km resolution at the equator,
  * ``oRRS18to6``: An RRS mesh with 6-km resolution at the poles and 18-km
    resolution at the equator.

Mapping files (see :ref:`config_mapping_files` below) and region mask files
(see :ref:`config_region_mask_files`) are provided from the
`E3SM public data repository`_ for these meshes.  For assistance with other
mesh resolutions, please contact the MPAS-Analysis developers.

Xarray and Dask
---------------

MPAS-Analysis makes extensive use of the `xarray package`_, which uses the
`dask package`_ internally to perform operations that are too large to fit
in memory.  While most tasks in MPAS-Analysis have moved away from opening
multi-file data sets using xarray in favor of concatinating these data sets
together using NCO tools, there are some legacy options that users can modify
if they experience errors related to dask::

  maxChunkSize = 10000

If an out of memory error occurs, it may first be worth reducing the number
of parallel tasks running (see :ref:`config_execute`) but if the error is
clearly related to dask (which might be the case, for example, if the error
occrus in the ``streamfunctionMOC`` task), you may wish to reduce the
``maxChunkSize``.  This will make tasks using dask slower but will reduce their
memory usage.

.. _`E3SM public data repository`: https://web.lcrc.anl.gov/public/e3sm/diagnostics/
.. _`xarray package`: https://xarray.pydata.org/en/stable/
.. _`dask package`: https://dask.pydata.org/en/latest/
