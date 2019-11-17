.. _config_climatology:

Climatology
===========

The ``[climatology]`` section of a configuration file contains options used to
determine the start and end years of climatologies, the reference years for
anomalies and to control remapping of climatologies to comparions grids::

  [climatology]
  ## options related to producing climatologies, typically to compare against
  ## observations and previous runs

  # the year from which to compute anomalies if not the start year of the
  # simulation.  This might be useful if a long spin-up cycle is performed and
  # only the anomaly over a later span of years is of interest.
  # anomalyRefYear = 249

  # the first year over which to average climatalogies
  startYear = 11
  # the last year over which to average climatalogies
  endYear = 20

  # The comparison lat/lon grid resolution in degrees
  comparisonLatResolution = 0.5
  comparisonLonResolution = 0.5

  # The comparison Antarctic polar stereographic grid size and resolution in km
  comparisonAntarcticStereoWidth = 6000.
  comparisonAntarcticStereoResolution = 10.

  # interpolation order for model and observation results. Likely values are
  #   'bilinear', 'neareststod' (nearest neighbor) or 'conserve'
  mpasInterpolationMethod = bilinear

  # should climatologies be performed with ncclimo or with xarray/dask
  useNcclimo = True

  # should remapping be performed with ncremap or with the Remapper class
  # directly in MPAS-Analysis
  useNcremap = True

  # The minimum weight of a destination cell after remapping. Any cell with
  # weights lower than this threshold will therefore be masked out.
  renormalizationThreshold = 0.01

  # if useNcclimo = False, the number of threads dask is allowed to spawn for
  # each process computing a climatology for a given month or season
  # Decrease this number if mpasClimatology* subtasks are running
  # out of available threads
  daskThreads = 2

  # if useNcclimo = False, the number of subprocesses that each climatology
  # subtask gets counted as occupying.
  # Increase this number if mpasClimatology* subtasks are running
  # out of memory, and fewer tasks will be allowed to run at once
  subprocessCount = 1


Start and End Year
------------------

A custom config file should specify a start and end year for climatologies.
Simulation data must exist for all 12 months of each year in this range.
Otherwise, the range will be reduced to those years with complete data and
a warning message will be displayed.


Anomaly Reference Year
----------------------

Anomalies between a climatology and the associated field from a reference year
is used in a few analysis tasks.  By default, the reference year is not
specified in the configuration file and it taken to be the start of the
simulation (determined from the contents of a restart file).  Under certain
circumstances (e.g. repetition of forcing data for several cycles, as in
the `Common Ocean Reference Experiments, CORE`_), it may be desirable to
specify a different year to use for computing anomalies::

  anomalyRefYear = 249

.. _config_remapping:

Remapping Options
-----------------

Climatologies are remapped from MPAS meshes and observations grids to common
comparison grids.  The remapping can be performed with any of three methods:
``bilinear``, ``neareststod`` (nearest neighbor) or ``conserve``.  Mapping
files are created with the `ESMF_RegridWeightGen tool`_.  The default method
is ``bilinear`` and these are the mapping files distributed from the
`E3SM public data repository`_.  The ``conseve`` method is know to be *much*
slower to compute and should only be used if it is necessary (e.g. because
remapped data will be checked for conservation).

MPAS-Analysis typically uses the `NCO`_ tool ``ncremap`` to perform remapping.
However, a python remapping capability is also available.  The user can force
remapping to use the python-based remapping by specifying::

  useNcremap = False

This capability is available largely for debugging purposes.

Remapped data typically only makes sense if it is renormalized after remapping.
For remapping of conserved quatntities like fluxes, renormalization would not
be desirable but for quantities like potential temperature, salinity and
potential density commonly used in MPAS-Anlaysis tasks, values become
physically meaningless near land boundaries and regions without data unless
renormalization is performed.  A threshold is needed to determine how much of a
cell's area on the output grid must contain valid data from the input grid or
mesh, below which that cell is considered invalid and is masked out of the
destination data set.  This threshold is specified as a fraction::

  renormalizationThreshold = 0.01

If noisy or unphysical values occur near maked regions on the comparison grid,
it might be necessary to increase this threshold.  If too much data appears to
be being masked out unnecessarily on the comparison grid, perhaps this value
should be made smaller.


Computing climatologies
-----------------------

MPAS-Analysis typically uses the `NCO`_ tool ``ncclimo`` to compute
climatologies.  For some large data sets on a single node, ``ncclimo`` runs
out of memory in ``bck`` mode but is painfully slow in ``serial`` mode and
wastes extra nodes in ``mpi`` mode.  (See :ref: `config_execute` for more on
configuring ``ncclimo``.)  For such cases, there is also an xarray/dask method
of computing climatologies::

  useNcremap = False


Other Options
-------------

* :ref:`config_comparison_grids`
* :ref:`dask_treads`

.. _`Common Ocean Reference Experiments, CORE`: http://data1.gfdl.noaa.gov/nomads/forms/mom4/CORE.html
.. _`ESMF_RegridWeightGen tool`: http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/ESMF_refdoc/node3.html#SECTION03020000000000000000
.. _`E3SM public data repository`: https://web.lcrc.anl.gov/public/e3sm/diagnostics/
.. _`NCO`: http://nco.sourceforge.net/nco.html
