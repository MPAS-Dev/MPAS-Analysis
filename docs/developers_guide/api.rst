#############
API reference
#############

This page provides an auto-generated summary of the MPAS-Analysis API. For
more details and examples, refer to the relevant chapters in the main part of
the documentation.

Top-level script: mpas_analysis
===============================

.. currentmodule:: mpas_analysis.__main__

.. autosummary::
   :toctree: generated/

   main
   build_analysis_list
   determine_analyses_to_generate
   add_task_and_subtasks
   update_generate
   run_analysis
   wait_for_task

Downloading data
================

.. currentmodule:: mpas_analysis.download_data

.. autosummary::
   :toctree: generated/

   download_analysis_data

Analysis tasks
==============

Base Class
----------

.. currentmodule:: mpas_analysis.shared

.. autosummary::
   :toctree: generated/

   AnalysisTask
   AnalysisTask.setup_and_check
   AnalysisTask.run_task
   AnalysisTask.run_after
   AnalysisTask.add_subtask
   AnalysisTask.run
   AnalysisTask.check_generate
   AnalysisTask.check_analysis_enabled
   AnalysisTask.set_start_end_date

Ocean tasks
-----------

.. currentmodule:: mpas_analysis.ocean

.. autosummary::
   :toctree: generated/

   ConservationTask
   ClimatologyMapSST
   ClimatologyMapSSS
   ClimatologyMapMLD
   ClimatologyMapMLDMinMax
   ClimatologyMapSSH
   ClimatologyMapEKE
   ClimatologyMapFluxes
   ClimatologyMapOHCAnomaly
   ClimatologyMapAntarcticMelt
   ClimatologyMapSose
   ClimatologyMapVel
   ClimatologyMapWoa
   ClimatologyMapArgoTemperature
   ClimatologyMapArgoSalinity
   ClimatologyMapWaves
   IndexNino34
   MeridionalHeatTransport
   OceanHistogram
   OceanRegionalProfiles
   RegionalTSDiagrams
   SoseTransects
   StreamfunctionMOC
   TimeSeriesOHCAnomaly
   TimeSeriesTemperatureAnomaly
   TimeSeriesSalinityAnomaly
   TimeSeriesSST
   TimeSeriesAntarcticMelt
   TimeSeriesOceanRegions
   TimeSeriesTransport
   WoaTransects
   WoceTransects

Ocean subtasks
--------------

.. currentmodule:: mpas_analysis.ocean.compute_anomaly_subtask

.. autosummary::
   :toctree: generated/

   ComputeAnomalySubtask

.. currentmodule:: mpas_analysis.ocean.plot_depth_integrated_time_series_subtask

.. autosummary::
   :toctree: generated/

   PlotDepthIntegratedTimeSeriesSubtask

.. currentmodule:: mpas_analysis.ocean.plot_hovmoller_subtask

.. autosummary::
   :toctree: generated/

   PlotHovmollerSubtask


Ocean utilities
---------------

.. currentmodule:: mpas_analysis.ocean.utility

.. autosummary::
   :toctree: generated/

   add_standard_regions_and_subset
   get_standard_region_names
   compute_zmid


Sea ice tasks
-------------

.. currentmodule:: mpas_analysis.sea_ice

.. autosummary::
   :toctree: generated/

   ClimatologyMapSeaIceConc
   ClimatologyMapSeaIceThick
   ClimatologyMapSeaIceSnowDepth
   ClimatologyMapSeaIceSnowiceFormation
   ClimatologyMapSeaIceSnowMelt
   ClimatologyMapSeaIcePondArea
   ClimatologyMapSeaIceAreaFractionRidge
   ClimatologyMapSeaIceVolumeRidge
   ClimatologyMapSeaIceAlbedo
   ClimatologyMapSeaIceProduction
   ClimatologyMapSeaIceMelting
   ClimatologyMapSeaIceAreaTendencyThermo
   ClimatologyMapSeaIceAreaTendencyTransp
   ClimatologyMapSeaIceVolumeTendencyThermo
   ClimatologyMapSeaIceVolumeTendencyTransp
   TimeSeriesSeaIce
   ClimatologyMapIcebergConc


Shared modules
==============

Reading MPAS Datasets
---------------------
.. currentmodule:: mpas_analysis.shared.io

.. autosummary::
   :toctree: generated/

   open_mpas_dataset

.. currentmodule:: mpas_analysis.shared.mpas_xarray

.. autosummary::
   :toctree: generated/

   mpas_xarray.open_multifile_dataset
   mpas_xarray.preprocess
   mpas_xarray.remove_repeated_time_index
   mpas_xarray.subset_variables

.. currentmodule:: mpas_analysis.shared.generalized_reader

.. autosummary::
   :toctree: generated/

   generalized_reader.open_multifile_dataset


Climatology
-----------
.. currentmodule:: mpas_analysis.shared.climatology

.. autosummary::
   :toctree: generated/

   get_comparison_descriptor
   get_remapper
   compute_monthly_climatology
   compute_climatology
   add_years_months_days_in_month
   get_unmasked_mpas_climatology_directory
   get_unmasked_mpas_climatology_file_name
   get_masked_mpas_climatology_file_name
   get_remapped_mpas_climatology_file_name

   MpasClimatologyTask
   MpasClimatologyTask.add_variables
   MpasClimatologyTask.get_file_name

   RemapMpasClimatologySubtask
   RemapMpasClimatologySubtask.setup_and_check
   RemapMpasClimatologySubtask.run_task
   RemapMpasClimatologySubtask.add_comparison_grid_descriptor
   RemapMpasClimatologySubtask.get_masked_file_name
   RemapMpasClimatologySubtask.get_remapped_file_name
   RemapMpasClimatologySubtask.customize_masked_climatology
   RemapMpasClimatologySubtask.customize_remapped_climatology

   RemapObservedClimatologySubtask
   RemapObservedClimatologySubtask.get_observation_descriptor
   RemapObservedClimatologySubtask.build_observational_dataset
   RemapObservedClimatologySubtask.get_file_name

Time Series
-----------
.. currentmodule:: mpas_analysis.shared.time_series

.. autosummary::
   :toctree: generated/

    cache_time_series
    compute_moving_avg_anomaly_from_start
    compute_moving_avg

    MpasTimeSeriesTask

Namelist and Streams Files
--------------------------
.. currentmodule:: mpas_analysis.shared.io.namelist_streams_interface

.. autosummary::
   :toctree: generated/

   convert_namelist_to_dict
   NameList.__init__
   NameList.__getattr__
   NameList.__getitem__
   NameList.get
   NameList.getint
   NameList.getfloat
   NameList.getbool

   StreamsFile.__init__
   StreamsFile.read
   StreamsFile.readpath
   StreamsFile.has_stream
   StreamsFile.find_stream

I/O Utilities
-------------
.. currentmodule:: mpas_analysis.shared.io

.. autosummary::
   :toctree: generated/

   utility.paths
   utility.make_directories
   utility.build_config_full_path
   utility.check_path_exists
   write_netcdf


Plotting
--------
.. currentmodule:: mpas_analysis.shared.plot

.. autosummary::
   :toctree: generated/

   timeseries_analysis_plot
   timeseries_analysis_plot_polar
   plot_polar_comparison
   plot_global_comparison
   plot_1D
   plot_vertical_section_comparison
   plot_vertical_section
   colormap.setup_colormap
   ticks.plot_xtick_format
   add_inset

   PlotClimatologyMapSubtask
   PlotClimatologyMapSubtask.set_plot_info



Projection
----------
.. currentmodule:: mpas_analysis.shared.projection

.. autosummary::
   :toctree: generated/

   get_pyproj_projection
   get_cartopy_projection


Regions
-------
.. currentmodule:: mpas_analysis.shared.regions

.. autosummary::
   :toctree: generated/

   compute_region_masks.ComputeRegionMasks
   compute_region_masks_subtask.ComputeRegionMasksSubtask
   compute_region_masks_subtask.get_feature_list

Timekeeping
-----------
.. currentmodule:: mpas_analysis.shared.timekeeping

.. autosummary::
   :toctree: generated/

   utility.get_simulation_start_time
   utility.string_to_datetime
   utility.string_to_relative_delta
   utility.string_to_days_since_date
   utility.days_to_datetime
   utility.datetime_to_days
   utility.date_to_days
   MpasRelativeDelta.MpasRelativeDelta

Transects
---------
.. currentmodule:: mpas_analysis.shared.transects

.. autosummary::
   :toctree: generated/

   compute_transect_masks_subtask.compute_mpas_transect_masks
   compute_transect_masks_subtask.ComputeTransectMasksSubtask

