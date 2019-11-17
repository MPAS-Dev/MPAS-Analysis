.. _config_transects:

Output Grids for Transects
==========================

Several tasks for producing transects use a common methodology for producing
the comparison grid for each transect::

    # The approximate horizontal resolution (in km) of each transect.  Latitude/
    # longitude between observation points will be subsampled at this interval.
    # Use 'obs' to indicate no subsampling.
    horizontalResolution = obs

    # The name of the vertical comparison grid.  Valid values are 'mpas' for the
    # MPAS vertical grid, 'obs' to use the locations of observations or
    # any other name if the vertical grid is defined by 'verticalComparisonGrid'
    # verticalComparisonGridName = obs
    verticalComparisonGridName = uniform_0_to_4000m_at_10m
    #verticalComparisonGridName = mpas

    # The vertical comparison grid if 'verticalComparisonGridName' is not 'mpas' or
    # 'obs'.  This should be numpy array of (typically negative) elevations (in m).
    verticalComparisonGrid = numpy.linspace(0, -4000, 401)

The ``horizontalResolution`` of all transects can be either ``obs`` or a number
of kilometers.  If ``obs``, model data are sampled at latitute and longitude
points corresponding to WOCE stations.  It a number of kilometers is given,
linear interpolation between observation points is performed with approximately
the requested resolution.  The distance between stations is always divided into
an integer number of segments of equal length so the resolution may be slightly
above or below ``horizontalResolution``.


The vertical grid is determined by two parameters,
``verticalComparisonGridName`` and ``verticalComparisonGrid``.  If
``verticalComparisonGridName = mpas``, the MPAS-Ocean vertical coordinate will
be interpolated horitontally from grid cell centers to the latitude and
longitude of each point along the transect, and the observations will be
interpolated vertically to the resulting grid.  If
``verticalComparisonGridName = obs``, the vertical grid of the observations
is used instead.  If ``verticalComparisonGridName`` is anything else, it is
taken to be the name of a user-defined vertical grid (best to make it
descriptive and unique, e.g. ``uniform_0_to_4000m_at_10m``) and
``verticalComparisonGrid`` should be assigned a valid array of positive-up
depth values (in the form of a python list or numpy array), e.g.::

    verticalComparisonGrid = numpy.linspace(0, -4000, 401)

produces points between 0 and -4000 m sampled every 10 m.

.. note::

  Some types of transects (e.g. those produce with geojson files) do not have
  a vertical grid associated with them (just horizontal latitude/longitude
  points), meaning that ``verticalComparisonGridName = obs`` is not a valid
  option for tasks with these transects.
