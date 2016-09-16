import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib.colors as cols

import numpy as np
import numpy.ma as ma
import xarray as xr
import datetime
#import calendar
from netCDF4 import Dataset as netcdf_dataset

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, remove_repeated_time_index
from ..shared.plot.plotting import plot_polar_comparison

def seaice_modelvsobs(config):

    indir = config.get('paths','archive_dir_ocn')
    plots_dir = config.get('paths','plots_dir')
    obsdir = config.get('paths','obs_seaicedir')

    casename = config.get('case','casename')

    remapfile = config.get('data','mpas_remapfile')
    climodir = config.get('data','mpas_climodir')

    climo_yr1 = config.getint('time','climo_yr1')
    climo_yr2 = config.getint('time','climo_yr2')
    yr_offset = config.getint('time','yr_offset')

    #climodir = "%s/%s" % (climodir,casename)
    climodir_regridded = "%s/mpas_regridded" % climodir
    if not os.path.isdir("%s" % climodir):
        print "\nClimatology directory does not exist. Create it...\n"
        os.mkdir("%s" % climodir)
    if not os.path.isdir("%s" % climodir_regridded):
        print "\nRegridded directory does not exist. Create it...\n"
        os.mkdir("%s" % climodir_regridded)

    print indir
    print climodir

    # Model climo (output) filenames
    climofiles = {}
    climofiles['winNH'] = "mpas-cice_climo.years%04d-%04d.jfm.nc" % (climo_yr1,climo_yr2)
    climofiles['sumNH'] = "mpas-cice_climo.years%04d-%04d.jas.nc" % (climo_yr1,climo_yr2)
    climofiles['winSH'] = "mpas-cice_climo.years%04d-%04d.djf.nc" % (climo_yr1,climo_yr2)
    climofiles['sumSH'] = "mpas-cice_climo.years%04d-%04d.jja.nc" % (climo_yr1,climo_yr2)
    climofiles['on']  = "mpas-cice_climo.years%04d-%04d.on.nc" % (climo_yr1,climo_yr2)
    climofiles['fm']  = "mpas-cice_climo.years%04d-%04d.fm.nc" % (climo_yr1,climo_yr2)

    # make a dictionary of the months in each climotology
    monthsInClim = {}
    monthsInClim['winNH'] = [1,2,3]
    monthsInClim['sumNH'] = [7,8,9]
    monthsInClim['winSH'] = [12,1,2]
    monthsInClim['sumSH'] = [6,7,8]
    monthsInClim['on'] = [10,11]
    monthsInClim['fm'] = [2,3]


    # Obs filenames
    obs_iceconc_filenames = {}
    obs_iceconc_filenames['winNH_NASATeam'] = "%s/SSMI/NASATeam_NSIDC0051/SSMI_NASATeam_gridded_concentration_NH_jfm.interp0.5x0.5.nc" % obsdir
    obs_iceconc_filenames['sumNH_NASATeam'] = "%s/SSMI/NASATeam_NSIDC0051/SSMI_NASATeam_gridded_concentration_NH_jas.interp0.5x0.5.nc" % obsdir
    obs_iceconc_filenames['winSH_NASATeam'] = "%s/SSMI/NASATeam_NSIDC0051/SSMI_NASATeam_gridded_concentration_SH_djf.interp0.5x0.5.nc" % obsdir
    obs_iceconc_filenames['sumSH_NASATeam'] = "%s/SSMI/NASATeam_NSIDC0051/SSMI_NASATeam_gridded_concentration_SH_jja.interp0.5x0.5.nc" % obsdir
    obs_iceconc_filenames['winNH_Bootstrap'] = "%s/SSMI/Bootstrap_NSIDC0079/SSMI_Bootstrap_gridded_concentration_NH_jfm.interp0.5x0.5.nc" % obsdir
    obs_iceconc_filenames['sumNH_Bootstrap'] = "%s/SSMI/Bootstrap_NSIDC0079/SSMI_Bootstrap_gridded_concentration_NH_jas.interp0.5x0.5.nc" % obsdir
    obs_iceconc_filenames['winSH_Bootstrap'] = "%s/SSMI/Bootstrap_NSIDC0079/SSMI_Bootstrap_gridded_concentration_SH_djf.interp0.5x0.5.nc" % obsdir
    obs_iceconc_filenames['sumSH_Bootstrap'] = "%s/SSMI/Bootstrap_NSIDC0079/SSMI_Bootstrap_gridded_concentration_SH_jja.interp0.5x0.5.nc" % obsdir
    obs_icethick_filenames = {}
    obs_icethick_filenames['onNH'] = "%s/ICESat/ICESat_gridded_mean_thickness_NH_on.interp0.5x0.5.nc" % obsdir
    obs_icethick_filenames['fmNH'] = "%s/ICESat/ICESat_gridded_mean_thickness_NH_fm.interp0.5x0.5.nc" % obsdir
    obs_icethick_filenames['onSH'] = "%s/ICESat/ICESat_gridded_mean_thickness_SH_on.interp0.5x0.5.nc" % obsdir
    obs_icethick_filenames['fmSH'] = "%s/ICESat/ICESat_gridded_mean_thickness_SH_fm.interp0.5x0.5.nc" % obsdir

    # Checks on directory/files existence:
    for climName in obs_iceconc_filenames:
        obs_filename = obs_iceconc_filenames[climName]
        if os.path.isfile("%s" % obs_filename) != True:
            raise SystemExit("Obs file %s not found. Exiting..." % obs_filename)
    for climName in obs_icethick_filenames:
        obs_filename = obs_icethick_filenames[climName]
        if os.path.isfile("%s" % obs_filename) != True:
            raise SystemExit("Obs file %s not found. Exiting..." % obs_filename)


    # Load data
    print "  Load sea-ice data..."
    print indir

    infiles = "".join([indir,"/am.mpas-cice.timeSeriesStatsMonthly.????-*.nc"])
    #infiles = "".join([indir,"/am.mpas-cice.timeSeriesStatsMonthly.001[5-6]-*.nc"])
    ds = xr.open_mfdataset(infiles, preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset,
                                timeSeriesStats=True,
                                timestr='timeSeriesStatsMonthly_avg_daysSinceStartOfSim_1',
                                onlyvars=['timeSeriesStatsMonthly_avg_iceAreaCell_1',
                                          'timeSeriesStatsMonthly_avg_iceVolumeCell_1']))
    ds = remove_repeated_time_index(ds)

    # Compute climatologies (first motnhly and then seasonally)
    print "  Compute seasonal climatologies..."
    time_start = datetime.datetime(yr_offset+climo_yr1,1,1)
    time_end = datetime.datetime(yr_offset+climo_yr2,12,31)
    ds_tslice = ds.sel(Time=slice(time_start,time_end))
    # check that each year has 24 months (?)
    monthly_clim = ds_tslice.groupby('Time.month').mean('Time')
    daysInMonth = [31,28,31,30,31,30,31,31,30,31,30,31]
    monthLetters = ['J','F','M','A','M','J','J','A','S','O','N','D']


    clims = {}
    for climName in monthsInClim:
        months = monthsInClim[climName]
        month = months[0]
        days = daysInMonth[month-1]
        climatology = days*monthly_clim.sel(month=month)
        totalDays = days
        for month in months[1:]:
            days = daysInMonth[month-1]
            climatology += days*monthly_clim.sel(month=month)
            totalDays += days
        climatology /= totalDays

        clims[climName] = climatology

    print "  Regrid fields to regular grid..."
    for climName in clims:
        # Save to netcdf files
        outFileName = "%s/%s" % (climodir,climofiles[climName])
        clims[climName].to_netcdf(outFileName)
        args = ["ncremap", "-P", "mpas", "-i", outFileName, "-m", remapfile, "-O", climodir_regridded]
        try:
            subprocess.check_call(args)
        except subprocess.CalledProcessError, e:
            print 'Error with call ', ' '.join(args)
            print e
            raise e

    print "  Make ice concentration plots..."
    suptitle = "Ice concentration"

    # interate over observations of sea-ice concentration
    first = True
    for climName in ['winNH','winSH','sumNH','sumSH']:
        hemisphere = climName[-2:]
        season = climName[:-2]

        if hemisphere == 'NH':
            plotProjection = 'npstere'
        else:
            plotProjection = 'spstere'

        clevsModelObs = config.getlist('seaice_modelvsobs','clevsModelObs_conc_%s'%season,listType=float)
        cmap = plt.get_cmap(config.get('seaice_modelvsobs','cmapModelObs'))
        cmapIndices = config.getlist('seaice_modelvsobs','cmapIndicesModelObs',listType=int)
        cmapModelObs = cols.ListedColormap(cmap(cmapIndices),"cmapModelObs")

        clevsDiff = config.getlist('seaice_modelvsobs','clevsDiff_conc_%s'%season,listType=float)
        cmap = plt.get_cmap(config.get('seaice_modelvsobs','cmapDiff'))
        cmapIndices = config.getlist('seaice_modelvsobs','cmapIndicesDiff',listType=int)
        cmapDiff = cols.ListedColormap(cmap(cmapIndices),"cmapDiff")

        lon0 = config.getfloat('seaice_modelvsobs','lon0_%s'%hemisphere)
        latmin = config.getfloat('seaice_modelvsobs','latmin_%s'%hemisphere)

        # Load in sea-ice data
        #  Model...
        # ice concentrations
        f = netcdf_dataset("%s/%s" % (climodir_regridded,climofiles[climName]),mode='r')
        iceconc = f.variables["timeSeriesStatsMonthly_avg_iceAreaCell_1"][:]
        if(first):
            lons = f.variables["lon"][:]
            lats = f.variables["lat"][:]
            print "Min lon: ", np.amin(lons), "Max lon: ", np.amax(lons)
            print "Min lat: ", np.amin(lats), "Max lat: ", np.amax(lats)
            Lons, Lats = np.meshgrid(lons, lats)
            first = False
        f.close()

        #  ...and observations
        # ice concentrations from NASATeam (or Bootstrap) algorithm
        for obsName in ['NASATeam','Bootstrap']:

            f = netcdf_dataset(obs_iceconc_filenames['%s_%s'%(climName, obsName)],mode='r')
            obs_iceconc = f.variables["AICE"][:]
            f.close()

            diff = iceconc - obs_iceconc

            monthsName = []
            for month in monthsInClim[climName]:
                monthsName.append(monthLetters[month-1])
            monthsName = ''.join(monthsName)

            plot_polar_comparison(
                config,
                Lons,
                Lats,
                iceconc,
                obs_iceconc,
                diff,
                cmapModelObs,
                clevsModelObs,
                cmapDiff,
                clevsDiff,
                title = "%s (%s, years %04d-%04d)" % (suptitle, monthsName, climo_yr1, climo_yr2),
                fileout = "%s/iceconc%s%s_%s_%s_years%04d-%04d.png" % (plots_dir, obsName, hemisphere,
                                                                       casename, monthsName, climo_yr1, climo_yr2),
                plotProjection = plotProjection,
                latmin = latmin,
                lon0 = lon0,
                modelTitle = "%s" % casename,
                obsTitle = "Observations (SSM/I %s)" % (obsName),
                diffTitle = "Model-Observations",
                cbarlabel = "%")


    print "  Make ice thickness plots..."
    # Plot Northern Hemisphere FM sea-ice thickness
    suptitle = "Ice thickness"
    # interate over observations of sea-ice thickness
    for climName in ['fm','on']:

        # Load in sea-ice data
        #  Model...
        # ice concentrations
        f = netcdf_dataset("%s/%s" % (climodir_regridded,climofiles[climName]),mode='r')
        icethick = f.variables["timeSeriesStatsMonthly_avg_iceVolumeCell_1"][:]
        f.close()

        monthsName = []
        for month in monthsInClim[climName]:
            monthsName.append(monthLetters[month-1])
        monthsName = ''.join(monthsName)

        for hemisphere in ['NH','SH']:
            #  ...and observations
            # ice concentrations from NASATeam (or Bootstrap) algorithm

            clevsModelObs = config.getlist('seaice_modelvsobs','clevsModelObs_thick_%s'%hemisphere,listType=float)
            cmap = plt.get_cmap(config.get('seaice_modelvsobs','cmapModelObs'))
            cmapIndices = config.getlist('seaice_modelvsobs','cmapIndicesModelObs',listType=int)
            cmapModelObs = cols.ListedColormap(cmap(cmapIndices),"cmapModelObs")

            clevsDiff = config.getlist('seaice_modelvsobs','clevsDiff_thick_%s'%hemisphere,listType=float)
            cmap = plt.get_cmap(config.get('seaice_modelvsobs','cmapDiff'))
            cmapIndices = config.getlist('seaice_modelvsobs','cmapIndicesDiff',listType=int)
            cmapDiff = cols.ListedColormap(cmap(cmapIndices),"cmapDiff")

            lon0 = config.getfloat('seaice_modelvsobs','lon0_%s'%hemisphere)
            latmin = config.getfloat('seaice_modelvsobs','latmin_%s'%hemisphere)

            f = netcdf_dataset(obs_icethick_filenames['%s%s'%(climName,hemisphere)],mode='r')
            obs_icethick = f.variables["HI"][:]
            f.close()
            # Mask thickness fields
            icethick[ icethick == 0 ] = ma.masked
            obs_icethick = ma.masked_values(obs_icethick,0)
            if hemisphere == 'NH':
                # Obs thickness should be nan above 86 (ICESat data)
                obs_icethick[ Lats > 86 ] = ma.masked
                plotProjection = 'npstere'
            else:
                plotProjection = 'spstere'

            diff = icethick - obs_icethick



            plot_polar_comparison(
                config,
                Lons,
                Lats,
                icethick,
                obs_icethick,
                diff,
                cmapModelObs,
                clevsModelObs,
                cmapDiff,
                clevsDiff,
                title = "%s (%s, years %04d-%04d)" % (suptitle, monthsName, climo_yr1, climo_yr2),
                fileout = "%s/icethick%s_%s_%s_years%04d-%04d.png" % (plots_dir, hemisphere,
                                                                      casename, monthsName, climo_yr1, climo_yr2),
                plotProjection = plotProjection,
                latmin = latmin,
                lon0 = lon0,
                modelTitle = "%s" % casename,
                obsTitle = "Observations (ICESat)",
                diffTitle = "Model-Observations",
                cbarlabel = "m")

