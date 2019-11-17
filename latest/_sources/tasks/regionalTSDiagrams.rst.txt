.. _task_regionalTSDiagrams:

regionalTSDiagrams
==================

An analysis task for plotting T-S (potential temperature vs. salinity)
diagrams of climatologis in ocean regions.

Component and Tags::

  component: ocean
  tags: climatology, regions, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

    [regionalTSDiagrams]
    ## options related to plotting T/S diagrams of ocean regions

    # the names of region groups to plot, each with its own section below
    regionGroups = ['Antarctic Regions', 'Ocean Basins']

    # a list of seasons to compute climatologies over
    seasons = ['ANN']

    # the number of threads dask is allowed to spawn for each process computing
    # the observational climatologies
    # Decrease this number if regionalTSDiagrams subtasks are running
    # out of available threads
    daskThreads = 8

    # the number of subprocesses that each observational climatology subtask gets
    # counted as occupying
    # Increase this number if regionalTSDiagrams subtasks are running
    # out of memory, and fewer tasks will be allowed to run at once
    subprocessCount = 4

    [TSDiagramsForAntarcticRegions]
    ## options related to plotting T/S diagrams of Antarctic regions

    # A geojson file containing Antarctic ocean regions.  Each region must have
    # 'zmin' and 'zmax' properites in addition to the usual properies for a region
    # in geometric_features
    regionMask = 'antarcticRegions.geojson'

    # list of regions to plot or ['all'] for all regions in the masks file.
    # See "regionNames" in the antarcticRegions masks file in
    # regionMaskSubdirectory for details.
    regionNames = ['all']

    # diagram type, either 'volumetric' or 'scatter', depending on if the points
    # should be binned the plot should show the volume fraction in each bin or
    # scattered points colored by their depth
    diagramType = volumetric

    # if diagramType == 'volumetric', the bin boundaries for T and S
    # if diagramType == 'scatter', only the min and max are important (and the
    #   bins are only used for computing neutral density contours)
    Tbins = numpy.linspace(-2.5, 4, 66)
    Sbins = numpy.linspace(33.8, 34.8, 51)

    # density contour interval
    rhoInterval = 0.1

    # The color map for depth or volume
    colormap = white_cmo_deep
    # The following is more appropriate if diagramType == 'scatter'
    # colormap = cmo.deep_r

    # The minimum and maximum depth over which fields are plotted, default is
    # to take these values from the geojson feature's zmin and zmax properties.
    # Add these to a custom config file to override the defaults.
    # zmin = -1000
    # zmax = -400

    # Obserational data sets to compare against
    obs = ['SOSE', 'WOA18']

    [TSDiagramsForOceanBasins]
    ## options related to plotting T/S diagrams of major ocean basins

    # A geojson file containing Antarctic ocean regions.  Each region must have
    # 'zmin' and 'zmax' properites in addition to the usual properies for a region
    # in geometric_features
    regionMask = 'oceanBasins.geojson'

    # list of regions to plot or ['all'] for all regions in the masks file.
    # See "regionNames" in the antarcticRegions masks file in
    # regionMaskSubdirectory for details.
    regionNames = ['all']

    # diagram type, either 'volumetric' or 'scatter', depending on if the points
    # should be binned the plot should show the volume fraction in each bin or
    # scattered points colored by their depth
    diagramType = volumetric

    # if diagramType == 'volumetric', the bin boundaries for T and S
    # if diagramType == 'scatter', only the min and max are important (and the
    #   bins are only used for computing neutral density contours)
    Tbins = numpy.linspace(-2.5, 16, 75)
    Sbins = numpy.linspace(33.8, 35.8, 101)

    # density contour interval
    rhoInterval = 0.2

    # The color map for depth or volume
    colormap = white_cmo_deep
    # The following is more appropriate if diagramType == 'scatter'
    # colormap = cmo.deep_r

    # The minimum and maximum depth over which fields are plotted.
    zmin = -1000
    zmax = 0

    # Obserational data sets to compare against
    obs = ['WOA18']

Region Groups
-------------

A list of groups of regions, each of which will get its own gallery group on
the resulting analysis webpage.  These can be any name without punctuation.
For each region group, there should be a corresponding
``TSDiagramsFor<RegionGroup>`` section of the config file, with any spaces
removed from the name of the region group.  By default, the only region group
is "Antarctic Regions".

Region Mask
-----------

The ``regionMask`` is a geojson file produce from the ``geometric_features``
packge.  It should include any number of ocean regions, each of which includes
properties ``zmin`` and ``zmax``.  Examples of how to create such a set of
features can be found in `antarctic_ocean_regions`_.

Region Names
------------

The ``regionNames`` can be set to ``['all']`` (the default) to plot all of the
regions in the ``regionMask`` file.  In the case of "Antarctic Regions", these
are::

  ["Southern Ocean", "Southern Ocean 60S", "Eastern Weddell Sea Shelf",
   "Eastern Weddell Sea Deep", "Western Weddell Sea Shelf",
   "Western Weddell Sea Deep", "Weddell Sea Shelf", "Weddell Sea Deep",
   "Bellingshausen Sea Shelf", "Bellingshausen Sea Deep", "Amundsen Sea Shelf",
   "Amundsen Sea Deep", "Eastern Ross Sea Shelf", "Eastern Ross Sea Deep",
   "Western Ross Sea Shelf", "Western Ross Sea Deep",
   "East Antarctic Seas Shelf", "East Antarctic Seas Deep"]

For "Ocean Basins", they are::

  ["Atlantic_Basin", "Pacific_Basin", "Indian_Basin", "Arctic_Basin",
   "Southern_Ocean_Basin", "Mediterranean_Basin", "Global Ocean",
   "Global Ocean 65N to 65S", "Global Ocean 15S to 15N"]

Diagram Type
------------

By default, a "volumetric" diagram is produced, where the volume of ocean water
in a reigon is binned in T and S space, and the volume of each bin is plotted.
This allows for more quantitative comparison with observation- and model-based
climatologies.

If ``diagramType`` is set to ``scatter``, a point cloud is plotted instead,
shaded by depth.  We anticipate that this will useful for plotting data sets
that are spatially scattered (e.g. the MEOP seal data), because each sample
does not correspond to a volume.  This type of dyagram may also be helpful for
comparison with publications that use scatter plots.

Bins and Contour Intervals
--------------------------
If ``diagramType = volumetric``, the number and spacing of potential
temperature and salinity bins is set in ``Tbins`` and ``Sbins``.  For
``diagramType = scatter``, ``Tbins`` and ``Sbins`` are used to make contour
plots of neutral density and are used to determine the bounds of the figure
in T/S space.  ``rhoInterval`` is the interval between contours of neutral
density.  ``zmin`` and ``zmax`` are the minimum and maximum depths (positive
up) of the ocean region.  If available (e.g. for "Antarctic Regions"), the
default is to read them from ``regionMask``.

Observations
------------
The ``obs`` option contains a list of the names of observational data sets.
Currently, "SOSE" and "WOA18" are the only data sets available, but we
anticipate adding several additional data sets in the near future.

:ref:`sose`
:ref:`woa18_t_s`

Other Config Options
--------------------

For more details on the remaining config options, see
 * :ref:`config_regions`
 * :ref:`config_seasons`
 * :ref:`config_colormaps`
 * :ref:`dask_treads`

Example Result
--------------

.. image:: examples/so_ts_diag.png
   :width: 500 px
   :align: center

.. _`antarctic_ocean_regions`: https://github.com/MPAS-Dev/geometric_features/tree/master/feature_creation_scripts/antarctic_ocean_regions
