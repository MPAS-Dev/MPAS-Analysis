.. _config_observations:

Ocean and Sea Ice Observations
==============================

The ``[oceanObservations]`` and ``[seaIceObservations]`` sections of a
configuration file contain options used to point to the observations files and
folders::

  [oceanObservations]
  ## options related to ocean observations with which the results will be compared

  # directory where ocean observations are stored
  baseDirectory = /dir/to/ocean/observations
  sstSubdirectory = SST
  sssSubdirectory = SSS
  mldSubdirectory = MLD
  ninoSubdirectory = Nino
  mhtSubdirectory = MHT
  meltSubdirectory = Melt
  soseSubdirectory = SOSE
  sshSubdirectory = SSH
  argoSubdirectory = ARGO
  schmidtkoSubdirectory = Schmidtko

  # interpolation order for observations. Likely values are
  #   'bilinear', 'neareststod' (nearest neighbor) or 'conserve'
  interpolationMethod = bilinear

  # The directories where observation climatologies will be stored if they need
  # to be computed.  If a relative path is supplied, it is relative to the output
  # base directory.  If an absolute path is supplied, this should point to
  # cached climatology files on the desired comparison grid.  If cached remapped
  # files are supplied, there is no need to provide cached files before
  # remapping.
  climatologySubdirectory = clim/obs
  remappedClimSubdirectory = clim/obs/remapped

  ...

  [seaIceObservations]
  ## options related to sea ice observations with which the results will be
  ## compared

  # directory where sea ice observations are stored
  baseDirectory = /dir/to/seaice/observations

  # interpolation order for observations. Likely values are
  #   'bilinear', 'neareststod' (nearest neighbor) or 'conserve'
  interpolationMethod = bilinear

  # The directories where observation climatologies will be stored if they need
  # to be computed.  If a relative path is supplied, it is relative to the output
  # base directory.  If an absolute path is supplied, this should point to
  # cached climatology files on the desired comparison grid.  If cached remapped
  # files are supplied, there is no need to provide cached files before
  # remapping.
  climatologySubdirectory = clim/obs
  remappedClimSubdirectory = clim/obs/remapped

Files and Directories
---------------------

The input directories are specified through a base directory and either
subdirectories or file names for each set of observations.  You will always
need to set ``baseDirectory`` for both ocean and sea ice observation to
a directory ending in ``observations/Ocean`` and ``observations/SeaIce``,
respectively.  The rest of these paths is determined by where you downloaded
the observationsal data from the `E3SM public data repository`_ (see the
:ref:`quick_start` for details).  Typically, the subdirectories and individual
files will not need to be changed except for debugging purposes.

The directories for storing cached datasets before and afer remapping
(specified in ``climatologySubdirectory`` and ``remappedClimSubdirectory``)
may be given any relative or absoltue path, but should typically be left as the
default values.

Remapping
---------

Observational climatologies are remapped from the native grid (typically
global latitude/longitude or Antarctic stereographic) to common
comparison grids.  The remapping can be performed with any of three methods:
``bilinear``, ``neareststod`` (nearest neighbor) or ``conserve``.  Mapping
files are created with the `ESMF_RegridWeightGen tool`_.  The default method
is ``bilinear`` and these are the mapping files distributed from the
`E3SM public data repository`_.  The ``conseve`` method is know to be *much*
slower to compute and should only be used if it is necessary (e.g. because
remapped data will be checked for conservation).

.. _`ESMF_RegridWeightGen tool`: http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/ESMF_refdoc/node3.html#SECTION03020000000000000000
.. _`E3SM public data repository`: https://web.lcrc.anl.gov/public/e3sm/diagnostics/
