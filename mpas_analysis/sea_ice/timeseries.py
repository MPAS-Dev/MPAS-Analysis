import numpy as np
import xarray as xr
import pandas as pd
import datetime

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, remove_repeated_time_index

from ..shared.plot.plotting  import timeseries_analysis_plot

def seaice_timeseries(config):

    varnames=['iceAreaCell','iceVolumeCell']

    plot_titles = {'iceAreaCell':'Sea-ice area',
                   'iceVolumeCell':'Sea-ice volume',
                   'iceThickness':'Sea-ice thickness'}

    units_dict = {'iceAreaCell':'[km$^2$]',
                  'iceVolumeCell':'[10$^3$ km$^3$]',
                   'iceThickness':'[m]'}

    obs_filenames = {'iceAreaCell':[config.get('seaIceData','obs_iceareaNH'),
                                    config.get('seaIceData','obs_iceareaSH')],
                     'iceVolumeCell':[config.get('seaIceData','obs_icevolNH'),
                                      config.get('seaIceData','obs_icevolSH')]}


    # Some plotting rules
    title_font_size = config.get('seaice_timeseries', 'title_font_size')

    indir = config.get('paths','archive_dir_ocn')
    meshfile = config.get('data','mpas_meshfile')

    casename = config.get('case','casename')
    ref_casename_v0 = config.get('case','ref_casename_v0')
    indir_v0data = config.get('paths','ref_archive_v0_seaicedir')

    compare_with_obs = config.getboolean('seaice_timeseries','compare_with_obs')

    plots_dir = config.get('paths','plots_dir')

    yr_offset = config.getint('time','yr_offset')

    N_movavg = config.getint('seaice_timeseries','N_movavg')

    print "  Load sea-ice data..."
    # Load mesh
    dsmesh = xr.open_dataset(meshfile)

    # Load data
    infiles = "".join([indir,'/am.mpas-cice.timeSeriesStatsMonthly.????-??-??.nc'])
    ds = xr.open_mfdataset(infiles, preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset,
                               timeSeriesStats=True,
                               timestr='timeSeriesStatsMonthly_avg_daysSinceStartOfSim_1',
                               onlyvars=['timeSeriesStatsMonthly_avg_iceAreaCell_1',
                                         'timeSeriesStatsMonthly_avg_iceVolumeCell_1']))
    ds = remove_repeated_time_index(ds)

    ds = ds.merge(dsmesh)

    year_start = (pd.to_datetime(ds.Time.min().values)).year
    year_end   = (pd.to_datetime(ds.Time.max().values)).year
    time_start = datetime.datetime(year_start,1,1)
    time_end   = datetime.datetime(year_end,12,31)

    if ref_casename_v0 != "None":
        infiles_v0data = "".join([indir_v0data,'/icevol.',ref_casename_v0,'.year*.nc'])
        ds_v0 = xr.open_mfdataset(infiles_v0data,
            preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset))
        ds_v0_tslice = ds_v0.sel(Time=slice(time_start,time_end))

    # Make Northern and Southern Hemisphere partition:
    areaCell = ds.areaCell
    ind_nh = ds.latCell > 0
    ind_sh = ds.latCell < 0
    areaCell_nh = areaCell.where(ind_nh)
    areaCell_sh = areaCell.where(ind_sh)

    for varname in varnames:
        obs_filenameNH = obs_filenames[varname][0]
        obs_filenameSH = obs_filenames[varname][1]
        plot_title = plot_titles[varname]
        units = units_dict[varname]

        print "  Compute NH and SH time series of %s..."%(varname)
        if varname == "iceThickCell":
            varnamefull = "timeSeriesStatsMonthly_avg_iceVolumeCell_1"
        else:
            varnamefull = "".join(["timeSeriesStatsMonthly_avg_",varname,"_1"])
        var = ds[varnamefull]

        var_nh = var.where(ind_nh)*areaCell_nh
        var_sh = var.where(ind_sh)*areaCell_sh

        ind_iceext = var > 0.15
        var_nh_iceext = var_nh.where(ind_iceext)
        var_sh_iceext = var_sh.where(ind_iceext)

        if varname == "iceAreaCell":
            var_nh = var_nh.sum('nCells')
            var_sh = var_sh.sum('nCells')
            var_nh = 1e-6*var_nh # m^2 to km^2
            var_sh = 1e-6*var_sh # m^2 to km^2
            var_nh_iceext = 1e-6*var_nh_iceext.sum('nCells')
            var_sh_iceext = 1e-6*var_sh_iceext.sum('nCells')
        elif varname == "iceVolumeCell":
            var_nh = var_nh.sum('nCells')
            var_sh = var_sh.sum('nCells')
            var_nh = 1e-3*1e-9*var_nh # m^3 to 10^3 km^3
            var_sh = 1e-3*1e-9*var_sh # m^3 to 10^3 km^3
        else:
            var_nh = var_nh.mean('nCells')/areaCell_nh.mean('nCells')
            var_sh = var_sh.mean('nCells')/areaCell_sh.mean('nCells')

        print "  Make plots..."

        xlabel = "Time [years]"

        if ref_casename_v0 != "None":
            figname_nh = "%s/%sNH_%s_%s.png" % (plots_dir,varname,casename,ref_casename_v0)
            figname_sh = "%s/%sSH_%s_%s.png" % (plots_dir,varname,casename,ref_casename_v0)
        else:
            figname_nh = "%s/%sNH_%s.png" % (plots_dir,varname,casename)
            figname_sh = "%s/%sSH_%s.png" % (plots_dir,varname,casename)

        title_nh = "%s (NH), %s (r)" % (plot_title,casename)
        title_sh = "%s (SH), %s (r)" % (plot_title,casename)

        if compare_with_obs:
            if varname == "iceAreaCell":
                title_nh = "%s\nSSM/I observations, annual cycle (k)" % title_nh
                title_sh = "%s\nSSM/I observations, annual cycle (k)" % title_sh
            elif varname == "iceVolumeCell":
                title_nh = "%s\nPIOMAS, annual cycle (k)" % title_nh
                title_sh = "%s\n" % title_sh

        if ref_casename_v0 != "None":
            title_nh = "%s\n %s (b)" % (title_nh,ref_casename_v0)
            title_sh = "%s\n %s (b)" % (title_sh,ref_casename_v0)


        if varname == "iceAreaCell":

            if compare_with_obs:
                ds_obs = xr.open_mfdataset(obs_filenameNH,
                    preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset))
                ds_obs = remove_repeated_time_index(ds_obs)
                var_nh_obs = ds_obs.IceArea
                var_nh_obs = replicate_cycle(var_nh,var_nh_obs)

                ds_obs = xr.open_mfdataset(obs_filenameSH,
                    preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset))
                ds_obs = remove_repeated_time_index(ds_obs)
                var_sh_obs = ds_obs.IceArea
                var_sh_obs = replicate_cycle(var_sh,var_sh_obs)

            if ref_casename_v0 != "None":
                infiles_v0data = "".join([indir_v0data,'/icearea.',ref_casename_v0,'.year*.nc'])
                ds_v0 = xr.open_mfdataset(infiles_v0data,
                    preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset))
                ds_v0_tslice = ds_v0.sel(Time=slice(time_start,time_end))
                var_nh_v0 = ds_v0_tslice.icearea_nh
                var_sh_v0 = ds_v0_tslice.icearea_sh

        elif varname == "iceVolumeCell":

            if compare_with_obs:
                ds_obs = xr.open_mfdataset(obs_filenameNH,
                    preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset))
                ds_obs = remove_repeated_time_index(ds_obs)
                var_nh_obs = ds_obs.IceVol
                var_nh_obs = replicate_cycle(var_nh,var_nh_obs)

                var_sh_obs = None

            if ref_casename_v0 != "None":
                infiles_v0data = "".join([indir_v0data,'/icevol.',ref_casename_v0,'.year*.nc'])
                ds_v0 = xr.open_mfdataset(infiles_v0data,
                    preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset))
                ds_v0_tslice = ds_v0.sel(Time=slice(time_start,time_end))
                var_nh_v0 = ds_v0_tslice.icevolume_nh
                var_sh_v0 = ds_v0_tslice.icevolume_sh

        if varname in ["iceAreaCell", "iceVolumeCell"]:
            if compare_with_obs:
                if ref_casename_v0 != "None":
                    vars_nh = [var_nh, var_nh_obs, var_nh_v0]
                    vars_sh = [var_sh, var_sh_obs, var_sh_v0]
                    lineStyles = ['r-', 'k-', 'b-']
                    lineWidths = [1.2, 1.2, 1.2]
                else:
                    # just v1 model and obs
                    vars_nh = [var_nh, var_nh_obs]
                    vars_sh = [var_sh, var_sh_obs]
                    lineStyles = ['r-', 'k-']
                    lineWidths = [1.2, 1.2]
            elif ref_casename_v0 != "None":
                # just v1 and v0 models
                vars_nh = [var_nh, var_nh_v0]
                vars_sh = [var_sh, var_sh_v0]
                lineStyles = ['r-', 'b-']
                lineWidths = [1.2, 1.2]

            if compare_with_obs or ref_casename_v0 != "None":
                # separate plots for nothern and southern hemispheres
                timeseries_analysis_plot(config, vars_nh, N_movavg, title_nh,
                                         xlabel, units, figname_nh,
                                         lineStyles=lineStyles,
                                         lineWidths=lineWidths,
                                         title_font_size=title_font_size)
                timeseries_analysis_plot(config, vars_sh, N_movavg, title_sh,
                                         xlabel, units, figname_sh,
                                         lineStyles=lineStyles,
                                         lineWidths=lineWidths,
                                         title_font_size=title_font_size)
            else:
                # we will combine north and south onto a single graph
                figname = "%s/%s.%s.png" % (plots_dir,casename,varname)
                title = "%s, NH (r), SH (k)\n%s" % (plot_title,casename)
                timeseries_analysis_plot(config, [var_nh, var_sh], N_movavg,
                                         title, xlabel, units, figname,
                                         lineStyles=['r-','k-'],
                                         lineWidths=[1.2, 1.2],
                                         title_font_size=title_font_size)

        elif varname == "iceThickCell":

            figname = "%s/%s.%s.png" % (plots_dir,casename,varname)
            title = "%s NH (r), SH (k)\n%s" % (plot_title,casename)
            timeseries_analysis_plot(config, [var_nh,var_sh], N_movavg, title,
                                      xlabel, units, figname,
                                      lineStyles=['r-','k-'],
                                      lineWidths=[1.2, 1.2],
                                      title_font_size=title_font_size)

        else:
            raise SystemExit("varname variable %s not supported for plotting" % varname)


def replicate_cycle(ds,ds_toreplicate):
    dsshift = ds_toreplicate.copy()
    shiftT = ((dsshift.Time.max() - dsshift.Time.min())
            + (dsshift.Time.isel(Time=1) - dsshift.Time.isel(Time=0)))
    nT = np.ceil((ds.Time.max() - ds.Time.min())/shiftT)

    # replicate cycle:
    for i in np.arange(nT):
        dsnew = ds_toreplicate.copy()
        dsnew['Time'] = dsnew.Time + (i+1)*shiftT
        dsshift = xr.concat([dsshift,dsnew], dim='Time')
    # constrict replicated ds_short to same time dimension as ds_long:
    dsshift = dsshift.sel(Time=ds.Time.values, method='nearest')
    return dsshift
