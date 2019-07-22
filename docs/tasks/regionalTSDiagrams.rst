.. _task_oceanRegionalProfiles:

oceanRegionalProfiles
=====================

An analysis task for plotting T-S (potential temperature vs. salinity)
diagrams of climatologis in ocean regions.

Component and Tags::

  component: ocean
  tags: climatology, regions, antarctic, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

  [regionalTSDiagrams]
  ## options related to plotting time series of groups of ocean regions

  # the names of region groups to plot, each with its own section below
  regionGroups = ['Antarctic Regions']

  # a list of seasons to compute climatologies over
  seasons = ['ANN']


  [TSDiagramsForAntarcticRegions]
  ## options related to plotting time series of Antarctic regions

  # A geojson file containing Antarctic ocean regions.  Each region must have
  # 'zmin' and 'zmax' properites in addition to the usual properies for a region
  # in geometric_features
  regionMask = 'antarcticRegions.geojson'

  # list of regions to plot or ['all'] for all regions in the masks file.
  # See "regionNames" in the antarcticRegions masks file in
  # regionMaskSubdirectory for details.
  regionNames = ['all']

  # min and max for T and S
  TBounds = [-2.5, 4]
  SBounds = [33.8, 34.8]

  # density contour interval
  rhoInterval = 0.1

  # The color map for depth
  colormap = cmo.deep_r

  # The minimum and maximum depth over which fields are plotted, default is
  # to take these values from the geojson feature's zmin and zmax properties.
  # Add these to a custom config file to override the defaults.
  # zmin = -1000
  # zmax = -400

  # Obserational data sets to compare against
  obs = ['SOSE']


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


Bounds
------
The range of potential temperature and salinity in each diagram are set by
``TBounds`` and ``SBounds``, respectively.  ``rhoInterval`` is not currently
used by will specify the interval between contours of density in a future
update.  ``zmin`` and ``zmax`` can be used to override the lower and upper
limits on depth, which are otherwise taken from ``regionMask``.

Observations
------------
The ``obs`` option contains a list of the names of observational data sets.
Currently, "SOSE" is the only data set available, but we anticipate adding
several additional data sets in the near future.

:ref:`sose`

Other Config Options
--------------------

For more details on the remaining config options, see
 * :ref:`config_regions`
 * :ref:`config_seasons`
 * :ref:`config_colormaps`


Example Result
--------------

.. image:: examples/so_ts_diag.png
   :width: 500 px
   :align: center

.. _`antarctic_ocean_regions`: https://github.com/MPAS-Dev/geometric_features/tree/master/feature_creation_scripts/antarctic_ocean_regions
