.. _config_transects:

Output Grids for Transects
==========================

Several tasks for producing transects use a common methodology for producing
the comparison grid for each transect::

    # The approximate horizontal resolution (in km) of each transect.  Latitude/
    # longitude between observation points will be subsampled at this interval.
    # Use 'obs' to indicate no subsampling. Use 'mpas' to indicate plotting of
    # model data on the native grid, in which case comparison with observations
    # will take place on the observation grid.
    #horizontalResolution = mpas
    #horizontalResolution = obs
    horizontalResolution = 5

    # The name of the vertical comparison grid.  Valid values are 'mpas' for the
    # MPAS vertical grid, 'obs' to use the locations of observations or
    # any other name if the vertical grid is defined by 'verticalComparisonGrid'.
    # If horizontalResolution is 'mpas', model data (both main and control) will be
    # plotted on the MPAS vertical grid, regardless of the comparison grid.
    #verticalComparisonGridName = mpas
    #verticalComparisonGridName = obs
    verticalComparisonGridName = uniform_0_to_4000m_at_10m

    # The vertical comparison grid if 'verticalComparisonGridName' is not 'mpas' or
    # 'obs'.  This should be numpy array of (typically negative) elevations (in m).
    verticalComparisonGrid = numpy.linspace(0, -4000, 401)

    # A range for the y axis (if any)
    verticalBounds = []

The ``horizontalResolution`` of all transects can be ``obs``, ``mpas`` or a
number of kilometers.  If ``obs``, model data are sampled at latitude and
longitude points corresponding to the observations.  If the horizontal grid
is ``mpas``, then the native MPAS-Ocean mesh is used for both the horizontal and
vertical grids.  If a number of kilometers is given, linear interpolation
between observation points is performed with approximately the requested
resolution.  The distance between observation points is always divided into an
integer number of segments of equal length so the resolution may be slightly
above or below ``horizontalResolution``.

The vertical grid is determined by two parameters,
``verticalComparisonGridName`` and ``verticalComparisonGrid``.  If
``verticalComparisonGridName = mpas``, but ``horizontalResoltuion`` is not
``mpas``, the MPAS-Ocean vertical coordinate will be interpolated horizontally
from grid cell centers to the latitude and longitude of each point along the
transect, and the observations will be interpolated vertically to the resulting
grid.  If ``verticalComparisonGridName = obs``, the vertical grid of the
observations is used instead.  If ``verticalComparisonGridName`` is anything
else, it is taken to be the name of a user-defined vertical grid (best to make
it descriptive and unique, e.g. ``uniform_0_to_4000m_at_10m``) and
``verticalComparisonGrid`` should be assigned a valid array of positive-up
depth values (in the form of a python list or numpy array), e.g.::

    verticalComparisonGrid = numpy.linspace(0, -4000, 401)

produces points between 0 and -4000 m sampled every 10 m.

``verticalBounds`` is a list of minimum and maximum limits for the vertical axis
of the transect.  The default is an empty list, which means ``matplotlib``
selects the axis limits to encompass the full range of the vertical grid.

.. note::

  Some types of transects (e.g. those produce with geojson files) do not have
  a vertical grid associated with them (just horizontal latitude/longitude
  points), meaning that ``verticalComparisonGridName = obs`` is not a valid
  option for tasks with these transects.
