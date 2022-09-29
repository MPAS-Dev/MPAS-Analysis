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

import xarray as xr
from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_obs_path


class ClimatologyMapSeaIceThick(AnalysisTask):
    """
    An analysis task for comparison of sea ice thickness against
    observations
    """
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

    def __init__(self, config, mpasClimatologyTask, hemisphere,
                 controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        hemisphere : {'NH', 'SH'}
            The hemisphere to plot

        controlConfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        taskName = f'climatologyMapSeaIceThick{hemisphere}'

        fieldName = 'seaIceThick'

        tags = ['climatology', 'horizontalMap', fieldName, 'publicObs']
        if hemisphere == 'NH':
            tags = tags + ['arctic']
        else:
            tags = tags + ['antarctic']

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceThick, self).__init__(
            config=config, taskName=taskName,
            componentName='seaIce',
            tags=tags)

        mpasFieldName = 'timeMonthly_avg_iceVolumeCell'
        iselValues = None

        sectionName = taskName

        if hemisphere == 'NH':
            hemisphereLong = 'Northern'
        else:
            hemisphereLong = 'Southern'

        # read in what seasons we want to plot
        seasons = config.getexpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError(f'config section {sectionName} does not contain '
                             f'a valid list of seasons')

        comparisonGridNames = config.getexpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError(f'config section {sectionName} does not contain '
                             f'a valid list of comparison grids')

        # the variable self.mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=f'{fieldName}{hemisphere}',
            variableList=[mpasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if controlConfig is None:
            refTitleLabel = 'Observations (ICESat)'
            galleryName = 'Observations: ICESat'
            diffTitleLabel = 'Model - Observations'
            refFieldName = 'seaIceThick'
        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = f'Control: {controlRunName}'
            refFieldName = mpasFieldName
            diffTitleLabel = 'Main - Control'

        for season in seasons:
            if controlConfig is None:
                obsFileName = build_obs_path(
                    config, 'seaIce',
                    relativePathOption=f'thickness{hemisphere}_{season}',
                    relativePathSection=sectionName)

                remapObservationsSubtask = RemapObservedThickClimatology(
                    parentTask=self, seasons=[season],
                    fileName=obsFileName,
                    outFilePrefix=f'{refFieldName}{hemisphere}_{season}',
                    comparisonGridNames=comparisonGridNames,
                    subtaskName=f'remapObservations{season}')
                self.add_subtask(remapObservationsSubtask)
            else:
                remapObservationsSubtask = None

            for comparisonGridName in comparisonGridNames:

                imageCaption = \
                    f'Climatology Map of {hemisphereLong}-Hemisphere ' \
                    f'Sea-Ice Thickness.'
                galleryGroup = \
                    f'{hemisphereLong}-Hemisphere Sea-Ice Thickness'
                # make a new subtask for this season and comparison grid
                subtaskName = f'plot{season}_{comparisonGridName}'

                subtask = PlotClimatologyMapSubtask(
                    parentTask=self, season=season,
                    comparisonGridName=comparisonGridName,
                    remapMpasClimatologySubtask=remapClimatologySubtask,
                    remapObsClimatologySubtask=remapObservationsSubtask,
                    subtaskName=subtaskName)

                subtask.set_plot_info(
                    outFileLabel=f'icethick{hemisphere}',
                    fieldNameInTitle='Sea ice thickness',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'm',
                    imageCaption=imageCaption,
                    galleryGroup=galleryGroup,
                    groupSubtitle=None,
                    groupLink=f'{hemisphere.lower()}_thick',
                    galleryName=galleryName,
                    maskValue=0)

                self.add_subtask(subtask)


class RemapObservedThickClimatology(RemapObservedClimatologySubtask):
    """
    A subtask for reading and remapping sea ice thickness observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def get_observation_descriptor(self, fileName):
        """
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(fileName=fileName,
                                                  latVarName='t_lat',
                                                  lonVarName='t_lon')
        return obsDescriptor

    def build_observational_dataset(self, fileName):
        """
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        Returns
        -------
        dsObs : ``xarray.Dataset``
            The observational dataset
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)
        dsObs = dsObs.rename({'HI': 'seaIceThick'})
        return dsObs
