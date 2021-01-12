.. _config_regions:

Regions
=======

Region Config Section
---------------------

The ``[regions]`` section of a configuration file contains options related
to regions either defined in MPAS components' online analysis or determined
within MPAS-Analysis using region mask files::

  [regions]
  ## options related to ocean regions used in several analysis modules

  # list of region names (needs to be in the same order as region indices in
  # time-series stats)
  regions = ['arctic', 'equatorial', 'so', 'nino3', 'nino4', 'nino3.4', 'global']
  # list of plot titles (needs to be in the same order as region indices in
  # time-series stats)
  plotTitles = ['Arctic', 'Equatorial (15S-15N)', 'Southern Ocean', 'Nino 3',
                'Nino 4', 'Nino 3.4', 'Global Ocean']


MPAS-Ocean currently has hard coded into its online analysis 7 regions in a
fixed order, as given in the ``regions`` option.  This should not be altered
unless corresponding changes to the MPAS-Ocean code have been made.

The corresponding ``plotTitles`` can be modified as desired to update how
these regions are named in plot titles and in gallery names on the resulting
web page.

.. _config_region_groups:

Region Groups
-------------

Currently, seven analysis tasks (:ref:`task_climatologyMapAntarcticMelt`,
:ref:`task_hovmollerOceanRegions`, :ref:`task_oceanRegionalProfiles`,
:ref:`task_regionalTSDiagrams`, :ref:`task_streamfunctionMOC`,
:ref:`task_timeSeriesAntarcticMelt`, and :ref:`task_timeSeriesOceanRegions`)
use masks that define regions in an MPAS mesh as part of their analysis.  Most
of these region group are defined in
:py:func:`mpas_analysis.shared.regions.compute_region_masks_subtask.get_region_info()`.
Several tasks (:ref:`task_hovmollerOceanRegions`,
:ref:`task_oceanRegionalProfiles`, :ref:`task_regionalTSDiagrams`, and
:ref:`task_timeSeriesOceanRegions`) can use any of the defined region groups.
Currently, available region groups are: ``Antarctic Regions``,
``Ocean Basins``, ``Ice Shelves``, and ``Ocean Subbasins``.

The option ``regionMaskSubdirectory`` in the ``[diagnostics]`` section specifies
the path to cached mask files for these region groups, typically
``diagnostics/mpas_analysis/region_masks``. Region masks for common MPAS Ocean
and Seaice meshes are supplied as part of the data from the
`E3SM public data repository <https://web.lcrc.anl.gov/public/e3sm/diagnostics/>`_
(see the :ref:`quick_start`).

If masks for for a given grid don't already exist in the cached mask location,
they will be generated automatically from the aggregation function from the
``geometric_features`` package, see
`Aggregate Existing Features <http://mpas-dev.github.io/geometric_features/stable/aggregation.html>`_.
The mask data will be stored in a geojson file with the region group's prefix
and date stamp (e.g. ``iceShelves20200621.geojson``).  Then, masks on the MPAS
Ocean and Seaice mesh will be computed, a process that can be time consuming for
large meshes.  To generate the masks in advance (using threading to speed up the
process), see the example utility script ``utility_scripts/make_region_mask.py``.

