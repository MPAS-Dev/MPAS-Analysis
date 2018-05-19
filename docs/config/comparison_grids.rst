.. _config_comparison_grids:

Comparison Grids
================

Currently, MPAS-Analysis supports two comparison grids to which both model
results and observations are remapped for comparison and plotting.  These are
``latlon``, a global latitude/longitude grid, and ``antarctic``, a
stereographic grid centered at the south pole.  Many analysis tasks support
both of these grids, in which case either or both can be specified in a list::

  # comparison grid(s) ('latlon', 'antarctic') on which to plot analysis
  comparisonGrids = ['latlon']

The resolution of each of these grids are set through a set of configuration
options that affect all tasks using these grids::

  [climatology]
  ## options related to producing climatologies, typically to compare against
  ## observations and previous runs

  ...

  # The comparison lat/lon grid resolution in degrees
  comparisonLatResolution = 0.5
  comparisonLonResolution = 0.5

  # The comparison Antarctic polar stereographic grid size and resolution in km
  comparisonAntarcticStereoWidth = 6000.
  comparisonAntarcticStereoResolution = 10.

These options can be altered by the user, in which case new weights for
remapping model results and observations to these grids will be generated on
the fly.  Mapping files for interpolating to the standard comparison grids from
both observations and standard MPAS meshes are provided when you download the
observations files.
