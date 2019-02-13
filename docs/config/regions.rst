.. _config_regions:

Regions
=======

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

Region Names
------------

MPAS-Ocean currently has hard coded into its online analysis 7 regions in a
fixed order, as given in the ``regions`` option.  This should not be altered
unless corresponding changes to the MPAS-Ocean code have been made.

The corresponding ``plotTitles`` can be modified as desired to update how
these regions are named in plot titles and in gallery names on the resulting
web page.


.. _config_region_mask_files:

Region Mask Files
-----------------

Currently, three analysis tasks (``streamfunctionMOC``,
``timeSeriesAntarcticMelt`` and ``oceanRegionalProfiles``) use masks that
define regions in an MPAS mesh as part of their analysis.  The option
``regionMaskSubdirectory`` specifies the path to these mask files, typically
``diagnostics/mpas_analysis/region_masks``.  Region  masks for common MPAS
Ocean and Seaice grids are supplied as part of the data from the
`E3SM public data repository`_ (see the :ref:`quick_start`).

If ice shelf or ocean basin masks for a given grid don't already exist in the
cached mask location, they will be generated automatically from the associated
geojson files (``iceShelves.geojson`` or ``oceanBasins.geojson``), a process
that can be very time consuming.  To generate them in advance (possibly using
threading to speed up the process), see the example utility script
``utility/make_region_mask.py``

.. _`E3SM public data repository`: https://web.lcrc.anl.gov/public/e3sm/diagnostics/
