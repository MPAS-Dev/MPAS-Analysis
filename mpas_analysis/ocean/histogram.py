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
from mpas_analysis.shared.io.utility import build_config_full_path, build_obs_path
from mpas_analysis.shared.climatology import compute_climatology, \
    get_unmasked_mpas_climatology_file_name

from mpas_analysis.shared.constants import constants
from mpas_analysis.shared.plot import histogram_analysis_plot, savefig
from mpas_analysis.shared.html import write_image_xml

class OceanHistogram(AnalysisTask):
    """
    Plots a histogram of a 2-d ocean variable.

    Attributes
    ----------
    variableDict : dict
        A dictionary of variables from the time series stats monthly output
        (keys), together with shorter, more convenient names (values)

    histogramFileName : str
        The name of the file where the histogram is stored

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
            taskName='oceanHistogram',
            componentName='ocean',
            tags=['climatology', 'regions', 'histogram', 'publicObs'])

        self.run_after(mpasClimatologyTask)
        self.mpasClimatologyTask = mpasClimatologyTask

        self.controlConfig = controlConfig


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

        #TODO make the filePrefix reflect the regionGroup and regionName
        self.filePrefix = f'histogram_{mainRunName}'

        regionGroups = config.getexpression(self.taskName, 'regionGroups')

        self.seasons = config.getexpression(self.taskName, 'seasons')
        self.variableList = config.getexpression(self.taskName, 'variableList')

        variableList = []
        self.variableDict = {}
        for var in self.variableList:
            key = f'timeMonthly_avg_{var}'
            variableList.append(key)
            self.variableDict[key] = var

            # Add xml file names for each season
            self.xmlFileNames = []
            for season in self.seasons:
                self.xmlFileNames.append(f'{self.plotsDirectory}/{self.filePrefix}_{var}_{season}.xml')

        # Specify variables and seasons to compute climology over
        self.mpasClimatologyTask.add_variables(variableList=variableList,
                                               seasons=self.seasons)
        #TODO fixup
        if controlConfig is None:
            if 'ssh' in self.variableList:
                refTitleLabel = 'Observations (AVISO Dynamic ' \
                    'Topography, 1993-2010)'

                observationsDirectory = build_obs_path(
                    config, 'ocean', 'sshSubdirectory')

                obsFileName = \
                    "{}/zos_AVISO_L4_199210-201012_20180710.nc".format(
                        observationsDirectory)
                print(obsFileName)


    def run_task(self):
        """
        Performs histogram analysis of the output of variables in variableList.
        """
        # Authors
        # -------
        # Carolyn Begeman, Adrian Turner, Xylar Asay-Davis

        self.logger.info("\nPlotting histogram of ocean vars...")

        config = self.config
        calendar = self.calendar
        seasons = self.seasons

        mainRunName = config.get('runs', 'mainRunName')

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

        if config.has_option(self.taskName, 'areaVarName'):
            areaVarName = config.get(self.taskName, 'areaVarName')
        else:
            areaVarName = 'areaCell'

        for season in seasons:
            inFileName = get_unmasked_mpas_climatology_file_name(
                config, season, self.componentName, op='avg')
            # Use xarray to open climatology dataset
            ds = xarray.open_dataset(inFileName)
            ds = self._multiply_var_by_area(ds, self.variableList, areaVarName=areaVarName)

            #TODO add region specification
            #ds.isel(nRegions=self.regionIndex))
            if config.has_option(self.taskName, 'lineColors'):
                lineColors = [config.get(self.taskName, 'mainColor')]
            else:
                lineColors = None
            lineWidths = [3]
            legendText = [mainRunName]

            title = mainRunName
            if config.has_option(self.taskName, 'titleFontSize'):
                titleFontSize = config.getint(self.taskName,
                                              'titleFontSize')
            else:
                titleFontSize = None
            if config.has_option(self.taskName, 'titleFontSize'):
                axisFontSize = config.getint(self.taskName,
                                             'axisFontSize')
            else:
                axisFontSize = None

            if config.has_option(self.taskName, 'defaultFontSize'):
                defaultFontSize = config.getint(self.taskName,
                                                'defaultFontSize')
            else:
                defaultFontSize = None
            if config.has_option(self.taskName, 'bins'):
                bins = config.getint(self.taskName, 'bins')
            else:
                bins = None

            yLabel = 'normalized Probability Density Function'

            for var in self.variableList:

                fields = [ds[f'{var}_{areaVarName}']]
                #Note: if we want to support 3-d variable histograms, we need to add depth masking

                #TODO add later
                #if plotControl:
                #    fields.append(refData.isel(nRegions=self.regionIndex))
                #    lineColors.append(config.get('histogram', 'controlColor'))
                #    lineWidths.append(1.2)
                #    legendText.append(controlRunName)
                #TODO make title more informative
                xLabel = f"{ds.ssh.attrs['long_name']} ({ds.ssh.attrs['units']})"

                histogram_analysis_plot(config, fields, calendar=calendar,
                                        title=title, xlabel=xLabel, ylabel=yLabel, bins=bins,
                                        lineColors=lineColors, lineWidths=lineWidths,
                                        legendText=legendText,
                                        titleFontSize=titleFontSize, defaultFontSize=defaultFontSize)

                outFileName = f'{self.plotsDirectory}/{self.filePrefix}_{var}_{season}.png'
                savefig(outFileName, config)

                #TODO should this be in the outer loop instead?
                caption = 'Normalized probability density function for SSH climatologies in the {} Region'.format(title)
                write_image_xml(
                    config=config,
                    filePrefix=f'{self.filePrefix}_{var}_{season}',
                    componentName='Ocean',
                    componentSubdirectory='ocean',
                    galleryGroup='Histograms',
                    groupLink=f'histogram{var}',
                    gallery=f'{var} Histogram',
                    thumbnailDescription=title,
                    imageDescription=caption,
                    imageCaption=caption)


    def _multiply_var_by_area(self, ds, variableList, areaVarName='areaCell', fractionVarName=None):

        """
        Compute a time series of the global mean water-column thickness.
        """

        restartFileName = self.runStreams.readpath('restart')[0]

        dsRestart = xarray.open_dataset(restartFileName)
        dsRestart = dsRestart.isel(Time=0)

        areaCell = dsRestart[areaVarName]
        print(f'shape(areaCell) = {numpy.shape(areaCell.values)}')

        # for convenience, rename the variables to simpler, shorter names
        ds = ds.rename(self.variableDict)

        for varName in variableList:
            varAreaName = f'{varName}_{areaVarName}'
            ds[varAreaName] = ds[varName] / areaCell
            if fractionVarName is not None:
                frac = ds[fractionVarName]
                ds[varAreaName] = ds[varName] * frac
            ds[varAreaName].attrs['units'] = 'm^2'
            ds[varAreaName].attrs['description'] = \
                f'{varName} multiplied by the cell area'

        return ds
