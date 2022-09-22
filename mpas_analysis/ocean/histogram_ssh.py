# -*- coding: utf-8 -*-
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
#
import os
import xarray
import numpy
import matplotlib.pyplot as plt

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io import open_mpas_dataset
from mpas_analysis.shared.io.utility import build_config_full_path
from mpas_analysis.shared.climatology import compute_climatology, \
    get_unmasked_mpas_climatology_file_name

from mpas_analysis.shared.constants import constants
from mpas_analysis.shared.plot import histogram_analysis_plot, savefig
from mpas_analysis.shared.html import write_image_xml

class HistogramSSH(AnalysisTask):
    """
    Plots a histogram of the global mean sea surface height.

    Attributes
    ----------
    variableDict : dict
        A dictionary of variables from the time series stats monthly output
        (keys), together with shorter, more convenient names (values)

    histogramFileName : str
        The name of the file where the ssh histogram is stored

    controlConfig : mpas_tools.config.MpasConfigParser
        Configuration options for a control run (if one is provided)

    filePrefix : str
        The basename (without extension) of the PNG and XML files to write out
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, regionMasksTask, controlConfig=None):

        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasHistogram: ``MpasHistogramTask``
            The task that extracts the time series from MPAS monthly output

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        controlConfig : mpas_tools.config.MpasConfigParser
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super().__init__(
            config=config,
            taskName='histogramSSH',
            componentName='ocean',
            tags=['climatology', 'regions', 'histogram', 'ssh', 'publicObs'])

        self.run_after(mpasClimatologyTask)
        self.mpasClimatologyTask = mpasClimatologyTask

        #self.histogramFileName = ''
        self.controlConfig = controlConfig
        #TODO make the filePrefix reflect the regionGroup and regionName
        self.filePrefix = 'histogram_ssh'

        self.variableDict = {}
        for var in ['ssh']:
            key = 'timeMonthly_avg_{}'.format(var)
            self.variableDict[key] = var

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        OSError
            If files are not present
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super().setup_and_check()

        config = self.config

        mainRunName = config.get('runs', 'mainRunName')
        regionGroups = config.getexpression(self.taskName, 'regionGroups')

        self.seasons = config.getexpression(self.taskName, 'seasons')
        #TODO add to config file
        #self.variableList = config.getexpression(self.taskName, 'variableList')
        self.variableList = ['timeMonthly_avg_ssh']

        # Add xml file names for each season
        self.xmlFileNames = []
        self.filePrefix = f'ssh_histogram_{mainRunName}'
        for season in self.seasons:
            self.xmlFileNames.append(f'{self.plotsDirectory}/{self.filePrefix}_{season}.xml')

        # Specify variables and seasons to compute climology over
        self.mpasClimatologyTask.add_variables(variableList=self.variableList,
                                               seasons=self.seasons)

    def run_task(self):
        """
        Performs histogram analysis of the output of sea-surface height
        (SSH).
        """
        # Authors
        # -------
        # Carolyn Begeman, Adrian Turner, Xylar Asay-Davis

        self.logger.info("\nPlotting histogram of SSH...")

        config = self.config
        calendar = self.calendar
        seasons = self.seasons

        #startYear = self.startYear
        #endYear = self.endYear
        #TODO determine whether this is needed
        #startDate = '{:04d}-01-01_00:00:00'.format(self.startYear)
        #endDate = '{:04d}-12-31_23:59:59'.format(self.endYear)

        mainRunName = config.get('runs', 'mainRunName')

        variableList = ['ssh']

        baseDirectory = build_config_full_path(
            config, 'output', 'histogramSubdirectory')
        print(f'baseDirectory={baseDirectory}')
        print(f'plotsDirectory={self.plotsDirectory}')

        # the variable mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        try:
            restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one'
                          ' restart file to plot T-S diagrams')
        #dsRestart = xarray.open_dataset(restartFileName)
        #dsRestart = dsRestart.isel(Time=0)

        for season in seasons:
            #TODO get the filename of the climatology file from the climatology task
            #TODO determine whether this is actually histogramSSH
            #TODO make sure that the climatology spans the appropriate years
            inFileName = get_unmasked_mpas_climatology_file_name(
                config, season, self.componentName, op='avg')
            # Use xarray to open climatology dataset
            ds = xarray.open_dataset(inFileName)
            ds = self._multiply_ssh_by_area(ds)

            #TODO add region specification
            #ds.isel(nRegions=self.regionIndex))

            fields = [ds[variableList[0]]]
            #TODO add depth masking

            if config.has_option('histogramSSH', 'lineColors'):
                lineColors = [config.get('histogramSSH', 'mainColor')]
            else:
                lineColors = None
            lineWidths = [3]
            legendText = [mainRunName]
            #TODO add later
            #if plotControl:
            #    fields.append(refData.isel(nRegions=self.regionIndex))
            #    lineColors.append(config.get('histogram', 'controlColor'))
            #    lineWidths.append(1.2)
            #    legendText.append(controlRunName)
            #TODO make title more informative
            title = mainRunName
            if config.has_option('histogramSSH', 'titleFontSize'):
                titleFontSize = config.getint('histogramSSH',
                                              'titleFontSize')
            else:
                titleFontSize = None

            if config.has_option('histogramSSH', 'defaultFontSize'):
                defaultFontSize = config.getint('histogramSSH',
                                                'defaultFontSize')
            else:
                defaultFontSize = None

            yLabel = 'normalized Probability Density Function'
            xLabel = f"{ds.ssh.attrs['long_name']} ({ds.ssh.attrs['units']})"

            histogram_analysis_plot(config, fields, calendar=calendar,
                                    title=title, xlabel=xLabel, ylabel=yLabel,
                                    lineColors=lineColors, lineWidths=lineWidths,
                                    legendText=legendText,
                                    titleFontSize=titleFontSize, defaultFontSize=defaultFontSize)

            #TODO whether this should be plotsDirectory or baseDirectory
            outFileName = f'{self.plotsDirectory}/{self.filePrefix}_{season}.png'
            print(f'outFileName={outFileName}')
            savefig(outFileName, config)

            #TODO should this be in the outer loop instead?
            caption = 'Normalized probability density function for SSH climatologies in the {} Region'.format(title)
            write_image_xml(
                config=config,
                filePrefix=f'{self.filePrefix}_{season}',
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Histograms',
                groupLink='histogramSSH',
                gallery='SSH Histogram',
                thumbnailDescription=title,
                imageDescription=caption,
                imageCaption=caption)


    def _multiply_ssh_by_area(self, ds):

        """
        Compute a time series of the global mean water-column thickness.
        """

        restartFileName = self.runStreams.readpath('restart')[0]

        dsRestart = xarray.open_dataset(restartFileName)
        dsRestart = dsRestart.isel(Time=0)

        #TODO load seaIceArea for sea ice histograms
        #landIceFraction = dsRestart.landIceFraction.isel(Time=0)
        areaCell = dsRestart.areaCell

        # for convenience, rename the variables to simpler, shorter names
        ds = ds.rename(self.variableDict)

        ds['sshAreaCell'] = \
            ds['ssh'] / areaCell
        ds.sshAreaCell.attrs['units'] = 'm^2'
        ds.sshAreaCell.attrs['description'] = \
            'Sea-surface height multiplied by the cell area'

        return ds
