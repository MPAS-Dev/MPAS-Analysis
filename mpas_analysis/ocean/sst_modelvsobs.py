from scipy.spatial import cKDTree

import matplotlib.pyplot as plt
import matplotlib.colors as cols

import numpy as np
import numpy.ma as ma
import xarray as xr
import datetime
from netCDF4 import Dataset as netcdf_dataset

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, remove_repeated_time_index

from ..shared.plot.plotting import plot_global_comparison

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
    
    figure_labels=['January','February','March','April','May','June','July','August','September','October','November','December','JFM','AMJ','JAS','OND','Annual']
    #outputTimes = [1,6,12,14,16] #times to output (0-11 = Jan - Dec, 12 - JFM, 13 - AMJ, 14-JAS, 15-OND, 16 - Annual)
    outputTimes = [12,14,16]

    if yr_offset < 1900:
        preIndustrial = True
    else:
        preIndustrial = False

    dLongitude = 1.   # Target longitude spacing
    dLatitude = 1.   # Target latitude spacing

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
    nLon = 360
    nLat = 180

    dsData = xr.open_mfdataset(obs_filename)
    if preIndustrial:
        time_start = datetime.datetime(1870,1,1)
        time_end = datetime.datetime(1900,12,31)
        preIndustrial_txt = "pre-industrial 1870-1900"
    else:
        time_start = datetime.datetime(1990,1,1)
        time_end = datetime.datetime(2011,12,31)
        preIndustrial_txt = "present-day 1990-2011"

    ds_tslice = dsData.sel(time=slice(time_start,time_end))
    monthly_clim_data = ds_tslice.groupby('time.month').mean('time')
    dinmonth = [31,28,31,30,31,30,31,31,30,31,30,31]
    daysarray = np.ones((12,monthly_clim_data.SST.values.shape[1],monthly_clim_data.SST.values.shape[2]))
    # The NaN values (last line in here) don't seem to be recognized by python:
    #for i, dval in enumerate(dinmonth):
    #    daysarray[i,:,:] = dval
    #    inds = np.where(np.isnan(monthly_clim_data.SST[i,:,:].values))
    #    daysarray[i,inds[0],inds[1]] = NaN

    JFMData = (dinmonth[0]*monthly_clim_data.sel(month=1) + dinmonth[1]*monthly_clim_data.sel(month=2) + dinmonth[2]*monthly_clim_data.sel(month=3)) / np.sum(dinmonth[:3])
    AMJData = (dinmonth[3]*monthly_clim_data.sel(month=4) + dinmonth[4]*monthly_clim_data.sel(month=5) + dinmonth[5]*monthly_clim_data.sel(month=6)) / np.sum(dinmonth[3:6])
    JASData = (dinmonth[6]*monthly_clim_data.sel(month=7) + dinmonth[7]*monthly_clim_data.sel(month=8) + dinmonth[8]*monthly_clim_data.sel(month=9)) / np.sum(dinmonth[6:9])
    ONDData = (dinmonth[9]*monthly_clim_data.sel(month=10) + dinmonth[10]*monthly_clim_data.sel(month=11) + dinmonth[11]*monthly_clim_data.sel(month=12)) / np.sum(dinmonth[9:12])
    ANNData = monthly_clim_data.mean(dim='month')

    latDATA, lonDATA = np.meshgrid(dsData.lat.values,dsData.lon.values)
    latDATA = latDATA.flatten()
    lonDATA = lonDATA.flatten()

    # Compute climatologies (first motnhly and then seasonally)
    time_start = datetime.datetime(yr_offset+climo_yr1,1,1)
    time_end = datetime.datetime(yr_offset+climo_yr2,12,31)
    ds_tslice = ds.sel(Time=slice(time_start,time_end))
    monthly_clim = ds_tslice.groupby('Time.month').mean('Time')

    JFM = (dinmonth[0]*monthly_clim.sel(month=1) + dinmonth[1]*monthly_clim.sel(month=2) + dinmonth[2]*monthly_clim.sel(month=3)) / np.sum(dinmonth[:3])
    AMJ = (dinmonth[3]*monthly_clim.sel(month=4) + dinmonth[4]*monthly_clim.sel(month=5) + dinmonth[5]*monthly_clim.sel(month=6)) / np.sum(dinmonth[3:6])
    JAS = (dinmonth[6]*monthly_clim.sel(month=7) + dinmonth[7]*monthly_clim.sel(month=8) + dinmonth[8]*monthly_clim.sel(month=9)) / np.sum(dinmonth[6:9])
    OND = (dinmonth[9]*monthly_clim.sel(month=10) + dinmonth[10]*monthly_clim.sel(month=11) + dinmonth[11]*monthly_clim.sel(month=12)) / np.sum(dinmonth[9:12])
    ANN = monthly_clim.mean(dim='month')

    def lon_lat_to_cartesian(lon, lat, R = 6371222):
        """
        calculates lon, lat coordinates of a point on a sphere with
        radius R
        """
        if max(abs(lon)) < 3.0*np.pi:
            lon_r = lon
            lat_r = lat
        else:
            lon_r = np.radians(lon)
            lat_r = np.radians(lat)

        x = R * np.cos(lat_r) * np.cos(lon_r)
        y = R * np.cos(lat_r) * np.sin(lon_r)
        z = R * np.sin(lat_r)
        return x,y,z
                
    def init_tree(lon_input, lat_input, dLon, dLat):                    
        lon_input = lon_input.flatten()
        lat_input = lat_input.flatten()

        if max(lon_input) > 180:
            inds = np.where(lon_input > 180)
            lon_input[inds] -= 360

        xs, ys, zs = lon_lat_to_cartesian(lon_input,lat_input)
        tree = cKDTree(zip(xs, ys, zs))

        lonVals = np.arange(-180. + dLon/2., 181 - dLon/2., dLon)
        latVals = np.arange(-90. + dLat/2., 91 - dLat/2., dLat)
        
        latTarg, lonTarg = np.meshgrid(latVals,lonVals)
        xt, yt, zt = lon_lat_to_cartesian(lonTarg.flatten(),latTarg.flatten())

        d, inds = tree.query(zip(xt, yt, zt), k = 1)

        return d, inds, lonTarg, latTarg

    def interp_fields(field, d, inds, lonTarg):
        #d, inds = tree.query(zip(xt, yt, zt), k = 10)
        #w = 1.0 / d**2
        #print d.shape
        #interpFld = np.sum(w * field.flatten()[inds],axis=1) / np.sum(w,axis=1)
        #interpFld.shape = lonTarg.shape

        #return interpFld
        return field.flatten()[inds].reshape(lonTarg.shape)

    dtr = 180/np.pi
    d2, inds2, lonTarg, latTarg = init_tree(lonCell*dtr,latCell*dtr, dLongitude, dLatitude)

    d, inds, lonTargD, latTargD = init_tree(lonDATA, latDATA, dLongitude, dLatitude)

    nLon = lonTarg.shape[0]
    nLat = lonTarg.shape[1]
    SST_MPAS = np.zeros((17,nLon,nLat))
    SST_DATA = np.zeros((17,nLon,nLat))
    SST_BIAS = np.zeros((17,nLon,nLat))

    for i in range(12):
        cur_mon = monthly_clim.sel(month=i+1)
        SST_MPAS[i,:,:] = interp_fields(cur_mon.time_avg_activeTracers_temperature[:].values, d2, inds2, lonTarg)

        cur_mon_data = monthly_clim_data.sel(month=i+1)
        SST_DATA[i,:,:] = interp_fields(cur_mon_data.SST.values.T.flatten(), d, inds, lonTargD)

    SST_MPAS[12,:,:] = interp_fields(JFM.time_avg_activeTracers_temperature[:].values, d2, inds2, lonTarg)
    SST_DATA[12,:,:] = interp_fields(JFMData.SST.values.T, d, inds, lonTargD)

    SST_MPAS[13,:,:] = interp_fields(AMJ.time_avg_activeTracers_temperature[:].values, d2, inds2, lonTarg)
    SST_DATA[13,:,:] = interp_fields(AMJData.SST.values.T, d, inds, lonTargD)

    SST_MPAS[14,:,:] = interp_fields(JAS.time_avg_activeTracers_temperature[:].values, d2, inds2, lonTarg)
    SST_DATA[14,:,:] = interp_fields(JASData.SST.values.T, d, inds, lonTargD)

    SST_MPAS[15,:,:] = interp_fields(OND.time_avg_activeTracers_temperature[:].values, d2, inds2, lonTarg)
    SST_DATA[15,:,:] = interp_fields(ONDData.SST.values.T, d, inds, lonTargD)

    SST_MPAS[16,:,:] = interp_fields(ANN.time_avg_activeTracers_temperature[:].values, d2, inds2, lonTarg)
    SST_DATA[16,:,:] = interp_fields(ANNData.SST.values.T, d, inds, lonTargD)

    for i in range(17):
        SST_BIAS[i,:,:] = SST_MPAS[i,:,:] - SST_DATA[i,:,:]
    SST_BIAS[16,:,:] = SST_MPAS[16,:,:] - SST_DATA[16,:,:]

    clevsModelObs = config.getlist('sst_modelvsobs','clevsModelObs',listType=float)
    cmap = plt.get_cmap(config.get('sst_modelvsobs','cmapModelObs'))
    cmapIndices = config.getlist('sst_modelvsobs','cmapIndicesModelObs',listType=int)
    cmapModelObs = cols.ListedColormap(cmap(cmapIndices),"cmapModelObs")

    clevsDiff = config.getlist('sst_modelvsobs','clevsDiff',listType=float)
    cmap = plt.get_cmap(config.get('sst_modelvsobs','cmapDiff'))
    cmapIndices = config.getlist('sst_modelvsobs','cmapIndicesDiff',listType=int)
    cmapDiff = cols.ListedColormap(cmap(cmapIndices),"cmapDiff")

    for i in outputTimes:
        plot_global_comparison(
            config,
            lonTarg,
            latTarg,
            SST_MPAS[i,:,:],
            SST_DATA[i,:,:],
            SST_BIAS[i,:,:],
            cmapModelObs,
            clevsModelObs,
            cmapDiff,
            clevsDiff,
            fileout = "%s/sstHADOI_%s_%s_years%04d-%04d.png" % (plots_dir, casename, figure_labels[i], climo_yr1, climo_yr2),
            title = "SST (%s, years %04d-%04d)" % (figure_labels[i], climo_yr1, climo_yr2),
            modelTitle = "%s" % casename,
            obsTitle = "Observations (Hadley/OI, %s)" % preIndustrial_txt,
            diffTitle = "Model-Observations",
            cbarlabel = r"$^o$C")    
   
#    nice_cmap = plt.get_cmap('viridis')
#    lev_cmap = nice_cmap([20,80,110,140,170,200,230,235,240])
#    new_cmap = cols.ListedColormap(lev_cmap,"mv_cmap")
#    norm = mpl.colors.BoundaryNorm(MLDclevs, new_cmap.N)
#    cs = m.contourf(x,y,SST_MPAS[i,:,:],cmap='viridis',spacing='uniform',levels= np.linspace(-1,30,15))
#    cbar = m.colorbar(cs,location='right',pad="15%",spacing='uniform',extendfrac='auto',
#                      extendrect='True',ticks=np.linspace(-1,30,15), boundaries=np.linspace(-1,30,15))
#
#    nice_cmap = plt.get_cmap('RdBu_r')
#    lev_cmap = nice_cmap([20,50,80,110,140,170,200,230,235,240])
#    new_cmap = cols.ListedColormap(lev_cmap,"mv_cmap")
#    norm = mpl.colors.BoundaryNorm(BIASclevs, new_cmap.N)
#    cs = m.contourf(x, y, SST_BIAS[i,:,:], cmap = 'RdBu_r', spacing = 'uniform', levels = np.linspace(-5,5,14))   
#    cbar = m.colorbar(cs,location='right',pad="15%",spacing='uniform',extendfrac='auto',
#                      extendrect='True',ticks=np.linspace(-5,5,14), boundaries=np.linspace(-5,5,14))
