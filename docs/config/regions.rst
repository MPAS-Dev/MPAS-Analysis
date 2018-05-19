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

  # Directory for region mask files
  regionMaskDirectory = /path/to/masks/

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

Currently, two analysis tasks (``streamfunctionMOC`` and
``timeSeriesAntarcticMelt``) use masks that define regions in an MPAS mesh
as part of their analysis.  ``regionMaskDirectory`` specifies the path to
these mask files.  Region masks for common MPAS Ocean and Seaice grids are
supplied as part of the data from the `E3SM public data repository`_ (see
the :ref:`quick_start`).  The path should typically end with
``mpas_analysis/maps``.

.. _`E3SM public data repository`: https://web.lcrc.anl.gov/public/e3sm/diagnostics/
