

import matplotlib.pyplot as plt
import matplotlib.colors as cols

import numpy as np
import xarray as xr
import datetime
from netCDF4 import Dataset as netcdf_dataset

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, remove_repeated_time_index

from ..shared.plot.plotting import plot_global_comparison

from ..shared.interpolation.interpolate import interp_fields, init_tree

def mld_modelvsobs(config):
    
    indir = config.get('paths','archive_dir_ocn')
    plots_dir = config.get('paths','plots_dir')
    obsdir = config.get('paths','obs_mlddir')
    
    casename = config.get('case','casename')
    meshfile = config.get('data','mpas_meshfile')

    climo_yr1 = config.getint('time','climo_yr1')
    climo_yr2 = config.getint('time','climo_yr2')
    yr_offset = config.getint('time','yr_offset')

    infiles = "".join([indir,"/am.mpas-o.timeSeriesStats.????-*.nc"])

    monthdictionary={'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,'Oct':10,
                         'Nov':11,'Dec':12,'JFM':np.array([1,2,3]),'AMJ':np.array([4,5,6]),'JAS':np.array([7,8,9]),
                         'OND':np.array([10,11,12]),'ANN':np.arange(1,13)}
                         
    outputTimes = config.getlist('mld_modelvsobs','comparisonTimes')

    obs_filename = "%s/holtetalley_mld_climatology.nc" % obsdir
    
      # Load in model data
    ds = xr.open_mfdataset(infiles, preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset,
                                timeSeriesStats=True,
                                timestr='time_avg_daysSinceStartOfSim',
                                onlyvars=['time_avg_dThreshMLD']))
    ds = remove_repeated_time_index(ds)

    f = netcdf_dataset(meshfile,mode='r')
    lonCell = f.variables["lonCell"][:]
    latCell = f.variables["latCell"][:]
    
# Load observation MLD data
    dsData = xr.open_mfdataset(obs_filename)
    dsData.iMONTH.values += 1
    
    latDATA, lonDATA = np.meshgrid(dsData.lat.values,dsData.lon.values)
    latDATA = latDATA.flatten()
    lonDATA = lonDATA.flatten()

    dinmonth = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    
    daysarray = np.ones((12,dsData.mld_dt_mean.values.shape[1],dsData.mld_dt_mean.values.shape[0]))

    for i, dval in enumerate(dinmonth):
        daysarray[i,:,:] = dval
        inds = np.where(np.isnan(dsData.mld_dt_mean[:,:,i].values.T))
        daysarray[i,inds[0],inds[1]] = np.NaN

    # Compute climatologies (first motnhly and then seasonally)
    time_start = datetime.datetime(yr_offset+climo_yr1,1,1)
    time_end = datetime.datetime(yr_offset+climo_yr2,12,31)
    ds_tslice = ds.sel(Time=slice(time_start,time_end))
    monthly_clim = ds_tslice.groupby('Time.month').mean('Time')

# initialize interpolation variables
    
    dLongitude = 1.   # Target longitude spacing
    dLatitude = 1.   # Target latitude spacing
    
    dtr = 180/np.pi
    d2, inds2, lonTarg, latTarg = init_tree(lonCell*dtr,latCell*dtr, dLongitude, dLatitude)
    d, inds, lonTargD, latTargD = init_tree(lonDATA, latDATA, dLongitude, dLatitude)

    nLon = lonTarg.shape[0]
    nLat = lonTarg.shape[1]
    
    MLD_MODEL = np.zeros((len(outputTimes),nLon,nLat))
    MLD_DATA = np.zeros((len(outputTimes),nLon,nLat))
    MLD_BIAS = np.zeros((len(outputTimes),nLon,nLat))

# Interpolate and compute biases

    for i, timestring in enumerate(outputTimes):
        
        monthsvalue = monthdictionary[timestring]
        
        if isinstance(monthsvalue, int):
            modeldata = monthly_clim.sel(month=monthsvalue).time_avg_dThreshMLD.values
            obsdata = dsData.sel(iMONTH=monthsvalue).mld_dt_mean.values
        
        else:
            modeldata = np.sum(dinmonth[monthsvalue-1]*monthly_clim.sel(month=monthsvalue).time_avg_dThreshMLD.values.T,axis=1) / np.sum(dinmonth[monthsvalue-1])
            obsdata = np.nansum(daysarray[monthsvalue-1,:,:]*dsData.sel(iMONTH=monthsvalue).mld_dt_mean.values.T,axis=0) / np.nansum(daysarray[monthsvalue-1,:,:],axis=0)

        MLD_MODEL[i,:,:] = interp_fields(modeldata, d2, inds2, lonTarg)

        MLD_DATA[i,:,:] = interp_fields(obsdata.flatten(), d, inds, lonTargD)
        

    for i in range(len(outputTimes)):
        MLD_BIAS[i,:,:] = MLD_MODEL[i,:,:] - MLD_DATA[i,:,:]

    clevsModelObs = config.getlist('mld_modelvsobs','clevsModelObs',listType=float)
    cmap = plt.get_cmap(config.get('mld_modelvsobs','cmapModelObs'))
    cmapIndices = config.getlist('mld_modelvsobs','cmapIndicesModelObs',listType=int)
    cmapModelObs = cols.ListedColormap(cmap(cmapIndices),"cmapModelObs")

    clevsDiff = config.getlist('mld_modelvsobs','clevsDiff',listType=float)
    cmap = plt.get_cmap(config.get('mld_modelvsobs','cmapDiff'))
    cmapIndices = config.getlist('mld_modelvsobs','cmapIndicesDiff',listType=int)
    cmapDiff = cols.ListedColormap(cmap(cmapIndices),"cmapDiff")

    for i in range(len(outputTimes)):
        plot_global_comparison(
            config,
            lonTarg,
            latTarg,
            MLD_MODEL[i,:,:],
            MLD_DATA[i,:,:],
            MLD_BIAS[i,:,:],
            cmapModelObs,
            clevsModelObs,
            cmapDiff,
            clevsDiff,
            fileout = "%s/mldHolteTalleyARGO_%s_%s_years%04d-%04d.png" % (plots_dir, casename, outputTimes[i], climo_yr1, climo_yr2),
            title = "MLD (%s, years %04d-%04d)" % (outputTimes[i], climo_yr1, climo_yr2),
            modelTitle = "%s" % casename,
            obsTitle = "Observations (HolteTalley density threshold MLD)",
            diffTitle = "Model-Observations",
            cbarlabel = "m")    
    
def sst_modelvsobs(config):

    indir = config.get('paths','archive_dir_ocn')
    plots_dir = config.get('paths','plots_dir')
    obsdir = config.get('paths','obs_sstdir')

    casename = config.get('case','casename')

    #remapfile = config.get('data','mpas_remapfile')
    #climodir = config.get('data','mpas_climodir')
    meshfile = config.get('data','mpas_meshfile')

    climo_yr1 = config.getint('time','climo_yr1')
    climo_yr2 = config.getint('time','climo_yr2')
    yr_offset = config.getint('time','yr_offset')
    #figure_tag = '.png' #type of figure to create

    infiles = "".join([indir,"/am.mpas-o.timeSeriesStats.????-*.nc"])
    obs_filename = "%s/MODEL.SST.HAD187001-198110.OI198111-201203.nc" % obsdir
    
    monthdictionary={'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,'Oct':10,'Nov':11,'Dec':12,'JFM':np.array([1,2,3]), \
                       'AMJ':np.array([4,5,6]),'JAS':np.array([7,8,9]),'OND':np.array([10,11,12]),'ANN':np.arange(1,13)}
    outputTimes = config.getlist('sst_modelvsobs','comparisonTimes',listType=str)

    if yr_offset < 1900:
        preIndustrial = True
    else:
        preIndustrial = False

    # Load in model data
    ds = xr.open_mfdataset(infiles, preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset,
                                timeSeriesStats=True,
                                timestr='time_avg_daysSinceStartOfSim',
                                onlyvars=['time_avg_activeTracers_temperature'],
                                selvals={'nVertLevels':1}))
    ds = remove_repeated_time_index(ds)

    f = netcdf_dataset(meshfile,mode='r')
    lonCell = f.variables["lonCell"][:]
    latCell = f.variables["latCell"][:]

    # Load in the observational data

    dsData = xr.open_mfdataset(obs_filename)
    if preIndustrial:
        time_start = datetime.datetime(1870,1,1)
        time_end = datetime.datetime(1900,12,31)
        preIndustrial_txt = "pre-industrial 1870-1900"
    else:
        time_start = datetime.datetime(1990,1,1)
        time_end = datetime.datetime(2011,12,31)
        preIndustrial_txt = "present-day 1990-2011"

    latDATA, lonDATA = np.meshgrid(dsData.lat.values,dsData.lon.values)
    latDATA = latDATA.flatten()
    lonDATA = lonDATA.flatten()
    
    ds_tslice = dsData.sel(time=slice(time_start,time_end))
    monthly_clim_data = ds_tslice.groupby('time.month').mean('time')
    dinmonth = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    daysarray = np.ones((12,monthly_clim_data.SST.values.shape[1],monthly_clim_data.SST.values.shape[2]))

    for i, dval in enumerate(dinmonth):
        daysarray[i,:,:] = dval
        inds = np.where(np.isnan(monthly_clim_data.SST[i,:,:].values))
        daysarray[i,inds[0],inds[1]] = np.NaN
        
    # Compute climatologies (first motnhly and then seasonally)
    time_start = datetime.datetime(yr_offset+climo_yr1,1,1)
    time_end = datetime.datetime(yr_offset+climo_yr2,12,31)
    ds_tslice = ds.sel(Time=slice(time_start,time_end))
    monthly_clim = ds_tslice.groupby('Time.month').mean('Time')
    
    dLongitude = 1.   # Target longitude spacing
    dLatitude = 1.   # Target latitude spacing
    
    dtr = 180/np.pi
    d2, inds2, lonTarg, latTarg = init_tree(lonCell*dtr,latCell*dtr, dLongitude, dLatitude)
    d, inds, lonTargD, latTargD = init_tree(lonDATA, latDATA, dLongitude, dLatitude)
    
    nLon = lonTarg.shape[0]
    nLat = lonTarg.shape[1]
    SST_MODEL = np.zeros((len(outputTimes),nLon,nLat))
    SST_DATA = np.zeros((len(outputTimes),nLon,nLat))
    SST_BIAS = np.zeros((len(outputTimes),nLon,nLat))

    for i, timestring in enumerate(outputTimes):
        
        monthsvalue = monthdictionary[timestring]
        
        if isinstance(monthsvalue, int):
            modeldata = monthly_clim.sel(month=monthsvalue).time_avg_dThreshMLD.values
            obsdata = monthly_clim_data.sel(month=monthsvalue).mld_dt_mean.values
        
        else:
            modeldata = np.sum(dinmonth[monthsvalue-1]*monthly_clim.sel(month=monthsvalue).time_avg_activeTracers_temperature.values.T,axis=1) / np.sum(dinmonth[monthsvalue-1])
            obsdata = np.nansum(daysarray[monthsvalue-1,:,:]*monthly_clim_data.sel(month=monthsvalue).SST.values,axis=0) / np.nansum(daysarray[monthsvalue-1,:,:],axis=0)
        
        SST_MODEL[i,:,:] = interp_fields(modeldata, d2, inds2, lonTarg)
        SST_DATA[i,:,:] = interp_fields(obsdata.T.flatten(), d, inds, lonTargD)
        
    for i in range(len(outputTimes)):
        SST_BIAS[i,:,:] = SST_MODEL[i,:,:] - SST_DATA[i,:,:]

    clevsModelObs = config.getlist('sst_modelvsobs','clevsModelObs',listType=float)
    cmap = plt.get_cmap(config.get('sst_modelvsobs','cmapModelObs'))
    cmapIndices = config.getlist('sst_modelvsobs','cmapIndicesModelObs',listType=int)
    cmapModelObs = cols.ListedColormap(cmap(cmapIndices),"cmapModelObs")

    clevsDiff = config.getlist('sst_modelvsobs','clevsDiff',listType=float)
    cmap = plt.get_cmap(config.get('sst_modelvsobs','cmapDiff'))
    cmapIndices = config.getlist('sst_modelvsobs','cmapIndicesDiff',listType=int)
    cmapDiff = cols.ListedColormap(cmap(cmapIndices),"cmapDiff")

    for i in range(len(outputTimes)):
        plot_global_comparison(
            config,
            lonTarg,
            latTarg,
            SST_MODEL[i,:,:],
            SST_DATA[i,:,:],
            SST_BIAS[i,:,:],
            cmapModelObs,
            clevsModelObs,
            cmapDiff,
            clevsDiff,
            fileout = "%s/sstHADOI_%s_%s_years%04d-%04d.png" % (plots_dir, casename, outputTimes[i], climo_yr1, climo_yr2),
            title = "SST (%s, years %04d-%04d)" % (outputTimes[i], climo_yr1, climo_yr2),
            modelTitle = "%s" % casename,
            obsTitle = "Observations (Hadley/OI, %s)" % preIndustrial_txt,
            diffTitle = "Model-Observations",
            cbarlabel = r"$^o$C")    
            
    
