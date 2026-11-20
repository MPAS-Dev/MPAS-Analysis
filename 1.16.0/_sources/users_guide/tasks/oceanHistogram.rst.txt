.. _task_oceanHistogram:

oceanHistogram
==============

An analysis task for plotting histograms of 2-d variables of climatologies
in ocean regions.

Component and Tags::

  component: ocean
  tags: climatology, histogram, regions, publicObs

Configuration Options
---------------------

The following configuration options are available for this task:

.. code-block:: cfg

    [oceanHistogram]
    ## options related to plotting histograms of climatologies of 2-d ocean
    ## variables

    # list of variables to plot
    variableList = []

    # list of observations to compare against
    obsList = ['AVISO']

    # list of ocean variables by which to weight variables in variable list
    weightList = []

    # list of regions to plot from the region list in [regions] below
    regionGroups = ['Ocean Basins']

    # list of region names within the region group listed above
    regionNames = []

    # Seasons to conduct analysis over
    # Note: start and end year will be inherited from climatology section
    seasons =  ['ANN']

    # Number of histogram bins
    bins = 40

Region Groups
-------------

A list of groups of regions, each of which will get its own gallery on
the resulting analysis web page.  See :ref:`config_region_groups` for
more information on the available region groups.  By default,
the only region group is ``'Ocean Basins'``.

Region Names
------------

The ``regionNames`` can be set to ``['all']`` to plot all of the
regions in the geojson file.  In the case of ``Antarctic Regions``, these
are:

.. code-block:: cfg

    ["Southern Ocean", "Southern Ocean 60S", "Eastern Weddell Sea Shelf",
     "Eastern Weddell Sea Deep", "Western Weddell Sea Shelf",
     "Western Weddell Sea Deep", "Weddell Sea Shelf", "Weddell Sea Deep",
     "Bellingshausen Sea Shelf", "Bellingshausen Sea Deep", "Amundsen Sea Shelf",
     "Amundsen Sea Deep", "Eastern Ross Sea Shelf", "Eastern Ross Sea Deep",
     "Western Ross Sea Shelf", "Western Ross Sea Deep",
     "East Antarctic Seas Shelf", "East Antarctic Seas Deep"]

For ``Ocean Basins``, they are:

.. code-block:: cfg

    ["Atlantic_Basin", "Pacific_Basin", "Indian_Basin", "Arctic_Basin",
     "Southern_Ocean_Basin", "Mediterranean_Basin", "Global Ocean",
     "Global Ocean 65N to 65S", "Global Ocean 15S to 15N"]

Variable List
-------------
The ``variableList`` option determines the variables to plot. Currently, only
2-d variables are supported but we envision extending to 3-d variables in the
near future. The variables are to be listed according to their registry names, with ``timeMontly_avg_`` prepended in the code.

Observations
------------
The ``obsList`` option contains a list of the names of observational data sets.
Currently, "AVISO" is the only data set available, but we anticipate adding
several additional data sets in the near future.

:ref:`aviso_ssh`

Weight List
-----------
The ``weightList`` option determines the variables to normalize by in plotting
the histogram. This is an optional feature. If used, ``weightList`` should be
of the same length as ``variableList``. We currently only support for 2-d
variables that are the same size as their corresponding variables in
``variableList`` and that are present in the restart file. For example, we use
``areaCell`` to weight cell-centered variables such as ``ssh``.

Bins
----
``bins`` sets the number of bins to include in the histogram.

Other Config Options
--------------------

Other config options include ``lineWidth``, ``mainColor``, ``obsColor``,
``controlColor``, ``titleFontSize``, ``defaultFontSize``. For more details on
the remaining config options, see :ref:`config_seasons`.

Example Result
--------------

.. image:: examples/histogram_ssh_aviso_atl.png
   :width: 500 px
   :align: center
