.. _task_geojsonTransects:

geojsonTransects
================

An analysis task for interpolating MPAS fields to transects specified by files
in ``geojson`` format..

Component and Tags::

    component: ocean
    tags: climatology, transect, geojson

Configuration Options
---------------------

The following configuration options are available for this task::

    [geojsonTransects]
    ## options related to plotting model transects at points determined by a
    ## geojson file.  To generate your own geojson file, go to:
    ##   http://geojson.io/
    ## and draw one or more polylines, then add a name to each:
    ##
    ##      "properties": {
    ##        "name": "My Favorite Name"
    ##      },
    ## and save the file as GeoJSON (say transects.geojson).  Finally, set the
    ## option:
    ##      geojsonFiles = ['transects.geojson']
    ## (giving an absolute path if necessary) in your custom config file.

    # a list of geojson files containing lat/lon points in LineStrings to be
    # plotted.  If relative paths are given, they are relative to the current
    # working directory.  The files must be listed in quotes, e.g.:
    # geojsonFiles = ['file1.geojson', '/path/to/file2.geojson']
    geojsonFiles = []

    # a list of dictionaries for each field to plot.  The dictionary includes
    # prefix (used for file names, task names and sections) as well as the mpas
    # name of the field, units for colorbars and a the name as it should appear
    # in figure titles and captions.
    fields =
        [{'prefix': 'temperature',
          'mpas': 'timeMonthly_avg_activeTracers_temperature',
          'units': r'$\degree$C',
          'titleName': 'Potential Temperature'},
         {'prefix': 'salinity',
          'mpas': 'timeMonthly_avg_activeTracers_salinity',
          'units': r'PSU',
          'titleName': 'Salinity'},
         {'prefix': 'potentialDensity',
          'mpas': 'timeMonthly_avg_potentialDensity',
          'units': r'kg m$^{-3}$',
          'titleName': 'Potential Density'},
         {'prefix': 'zonalVelocity',
          'mpas': 'timeMonthly_avg_velocityZonal',
          'units': r'm s$^{-1}$',
          'titleName': 'Zonal Velocity'},
         {'prefix': 'meridionalVelocity',
          'mpas': 'timeMonthly_avg_velocityMeridional',
          'units': r'm s$^{-1}$',
          'titleName': 'Meridional Velocity'}]

    # Times for comparison times (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct,
    # Nov, Dec, JFM, AMJ, JAS, OND, ANN)
    seasons =  ['ANN']

    # The approximate horizontal resolution (in km) of each transect.  Latitude/
    # longitude between observation points will be subsampled at this interval.
    # Use 'obs' to indicate no subsampling.
    horizontalResolution = 5

    # The name of the vertical comparison grid.  Valid values are 'mpas' for the
    # MPAS vertical grid or any other name if the vertical grid is defined by
    # 'verticalComparisonGrid'
    #verticalComparisonGridName = mpas
    verticalComparisonGridName = uniform_0_to_4000m_at_10m

    # The vertical comparison grid if 'verticalComparisonGridName' is not 'mpas'.
    # This should be numpy array of (typically negative) elevations (in m).
    verticalComparisonGrid = numpy.linspace(0, -4000, 401)

    # The minimum weight of a destination cell after remapping. Any cell with
    # weights lower than this threshold will therefore be masked out.
    renormalizationThreshold = 0.01


    [geojsonTemperatureTransects]
    ## options related to plotting geojson transects of potential temperature

    # colormap for model/observations
    colormapNameResult = RdYlBu_r
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the SemiLogNorm
    normArgsResult = {'vmin': -2., 'vmax': 30.}
    # place the ticks automatically by default
    # colorbarTicksResult = numpy.linspace(-2., 2., 9)

    # colormap for differences
    colormapNameDifference = balance
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the SemiLogNorm
    normArgsDifference = {'vmin': -2., 'vmax': 2.}
    # place the ticks automatically by default
    # colorbarTicksDifference = numpy.linspace(-2., 2., 9)


    [geojsonSalinityTransects]
    ## options related to plotting geojson transects of salinity

    # colormap for model/observations
    colormapNameResult = haline
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the SemiLogNorm
    normArgsResult = {'vmin': 30, 'vmax': 39.0}
    # place the ticks automatically by default
    # colorbarTicksResult = numpy.linspace(34.2, 35.2, 9)

    # colormap for differences
    colormapNameDifference = balance
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the SemiLogNorm
    normArgsDifference = {'vmin': -0.5, 'vmax': 0.5}
    # place the ticks automatically by default
    # colorbarTicksDifference = numpy.linspace(-0.5, 0.5, 9)


    [geojsonPotentialDensityTransects]
    ## options related to plotting geojson transects of potential density

    # colormap for model/observations
    colormapNameResult = Spectral_r
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the norm
    normArgsResult = {'vmin': 1026.5, 'vmax': 1028.}
    # place the ticks automatically by default
    # colorbarTicksResult = numpy.linspace(1026., 1028., 9)

    # colormap for differences
    colormapNameDifference = balance
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the norm
    normArgsDifference = {'vmin': -0.3, 'vmax': 0.3}
    # place the ticks automatically by default
    # colorbarTicksDifference = numpy.linspace(-0.3, 0.3, 9)


    [geojsonZonalVelocityTransects]
    ## options related to plotting geojson transects of zonal velocity

    # colormap for model/observations
    colormapNameResult = delta
    # color indices into colormapName for filled contours
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the norm
    normArgsResult = {'vmin': -0.2, 'vmax': 0.2}
    # determine the ticks automatically by default, uncomment to specify
    # colorbarTicksResult = numpy.linspace(-0.2, 0.2, 9)

    # colormap for differences
    colormapNameDifference = balance
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the norm
    normArgsDifference = {'vmin': -0.2, 'vmax': 0.2}
    # determine the ticks automatically by default, uncomment to specify
    # colorbarTicksDifference = numpy.linspace(-0.2, 0.2, 9)


    [geojsonMeridionalVelocityTransects]
    ## options related to plotting geojson transects of meridional velocity

    # colormap for model/observations
    colormapNameResult = delta
    # color indices into colormapName for filled contours
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the norm
    normArgsResult = {'vmin': -0.2, 'vmax': 0.2}
    # determine the ticks automatically by default, uncomment to specify
    # colorbarTicksResult = numpy.linspace(-0.2, 0.2, 9)

    # colormap for differences
    colormapNameDifference = balance
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the norm
    normArgsDifference = {'vmin': -0.2, 'vmax': 0.2}
    # determine the ticks automatically by default, uncomment to specify
    # colorbarTicksDifference = numpy.linspace(-0.2, 0.2, 9)

Geojson Files
-------------

This task takes a list of geojson file names (supplied as a python list of
``str``::

  geojsonFiles = ['file1.geojson', '/path/to/file2.geojson']

Transects are specified by ``LineString`` objects in the files.  Some examples
are provided in the `MPAS geometric features repository`_.  You can also
generate your own very easily at To generate your own geojson file, go to
`geojson.io`_ and draw one or more polylines, then add a name to each::

  ...
  "properties": {
    "name": "My Favorite Name"
  },
  ...

and save the file as GeoJSON (say transects.geojson).  Finally, set the
option::

    geojsonFiles = ['transects.geojson']

(giving an absolute path if necessary) in your custom config file.


Fields
------

Since there are no observations associated with geojson transects, you are
free to choose which MPAS fields you would like to plot.  These fields are
provided as a python dictionary.  The keys are names for the fields (anything
you would like use as a prefix on files and subtask names, best if it does
not contain spaces).  The values are python dictionaries.  The values
associate with the ``mpas`` key are the names of the 3D fields where transects
are desired.  The ``units`` entry indicates the units to display on the
colorbar.  The ``titleName`` entry specifies the name of the field to include
in plot titles and captions.

Each field must have a corresponding section in the config file defining its
color maps.  For example, ``temperature`` has an associated
``[geojsonTemperatureTransect]`` section.

Other Options
-------------

For details on the remaining configuration options, see:
 * :ref:`config_transects`
 * :ref:`config_remapping`
 * :ref:`config_colormaps`
 * :ref:`config_seasons`

Example Result
--------------

.. image:: examples/geojson_transect.png
   :width: 500 px
   :align: center

.. _`MPAS geometric features repository`: https://github.com/MPAS-Dev/geometric_features
.. _`geojson.io`: http://geojson.io/
