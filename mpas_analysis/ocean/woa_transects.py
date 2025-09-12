# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
import os
from collections import OrderedDict

import numpy as np
import xarray as xr

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.compute_transects_subtask import \
    TransectsObservations

from mpas_analysis.ocean.compute_transects_with_vel_mag import \
    ComputeTransectsWithVelMag

from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, build_obs_path
from mpas_analysis.shared.io import write_netcdf_with_fill


class WoaTransects(AnalysisTask):
    """
    Plot model output at transects and various longitudes around Antarctica,
    compared against the WOA23 climatology
    """

    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):
        """
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        config : tranche.Tranche
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted
            as a transect

        controlconfig : tranche.Tranche, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        tags = ['climatology', 'transect', 'woa', 'publicObs', 'antarctic']

        # call the constructor from the base class (AnalysisTask)
        super().__init__(
            config=config, taskName='woaTransects',
            componentName='ocean',
            tags=tags)

        sectionName = self.taskName

        seasons = config.getexpression(sectionName, 'seasons')

        horizontalResolution = config.get(sectionName, 'horizontalResolution')

        verticalComparisonGridName = config.get(sectionName,
                                                'verticalComparisonGridName')

        if verticalComparisonGridName in ['mpas', 'obs']:
            verticalComparisonGrid = None
        else:
            verticalComparisonGrid = config.getexpression(
                sectionName, 'verticalComparisonGrid', allow_numpy=True)

        verticalBounds = config.getexpression(sectionName, 'verticalBounds')

        longitudes = _get_longitudes(config)

        fields = \
            [{'prefix': 'temperature',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'units': r'$\degree$C',
              'titleName': 'Potential Temperature',
              'obsFilePrefix': 'pot_temp',
              'obsFieldName': 'pt_an'},
             {'prefix': 'salinity',
              'mpas': 'timeMonthly_avg_activeTracers_salinity',
              'units': r'PSU',
              'titleName': 'Salinity',
              'obsFilePrefix': 'salinity',
              'obsFieldName': 's_an'}]

        fieldList = config.getexpression(sectionName, 'fieldList')
        fields = [field for field in fields if field['prefix'] in fieldList]

        variableList = [field['mpas'] for field in fields
                        if field['mpas'] != 'velMag']

        transectCollectionName = 'WOA_transects'
        if horizontalResolution not in ['obs', 'mpas']:
            transectCollectionName = \
                f'{transectCollectionName}_{horizontalResolution}km'

        transectsObservations = WoaTransectsObservations(
            config, horizontalResolution,
            transectCollectionName, fields)

        computeTransectsSubtask = ComputeTransectsWithVelMag(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='WOA_transects',
            transectCollectionName=transectCollectionName,
            variableList=variableList,
            seasons=seasons,
            obsDatasets=transectsObservations,
            verticalComparisonGridName=verticalComparisonGridName,
            verticalComparisonGrid=verticalComparisonGrid)

        plotObs = controlConfig is None
        if plotObs:

            refTitleLabel = 'Observations (WOA23 1991-2020)'

            diffTitleLabel = 'Model - State Estimate'

        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = f'Control: {controlRunName}'

            diffTitleLabel = 'Main - Control'

        for field in fields:
            fieldPrefix = field['prefix']
            fieldPrefixUpper = fieldPrefix[0].upper() + fieldPrefix[1:]
            for lon in longitudes:
                transectName = f'lon_{lon}'

                for season in seasons:
                    outFileLabel = f'WOA_{fieldPrefix}_'
                    if plotObs:
                        refFieldName = field['obsFieldName']
                    else:
                        refFieldName = field['mpas']

                    fieldNameInTytle = \
                        fr'{field["titleName"]} at {lon}$\degree$ Lon.'
                    # make a new subtask for this season and comparison grid
                    subtask = PlotTransectSubtask(self, season, transectName,
                                                  fieldPrefix,
                                                  computeTransectsSubtask,
                                                  plotObs, controlConfig)

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=fieldNameInTytle,
                        mpasFieldName=field['mpas'],
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=field['units'],
                        imageCaption=f'{fieldNameInTytle} {season}',
                        galleryGroup='WOA23 Transects',
                        groupSubtitle=None,
                        groupLink='woa_transects',
                        galleryName=field['titleName'],
                        configSectionName=f'woa{fieldPrefixUpper}Transects',
                        verticalBounds=verticalBounds)

                    self.add_subtask(subtask)


class WoaTransectsObservations(TransectsObservations):
    """
    A class for loading and manipulating WOA transect data

    Attributes
    ----------

    fields : list of dict
        dictionaries defining the fields with associated WOA data
    """

    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, horizontalResolution,
                 transectCollectionName, fields):

        """
        Construct the object, setting the observations dictionary to None.

        Parameters
        ----------
        config : tranche.Tranche
            Configuration options

        horizontalResolution : str
            'obs' for the obs as they are or a size in km if subdivision is
            desired.

        transectCollectionName : str
            A name that describes the collection of transects (e.g. the name
            of the collection of observations) used to name the
            destination "mesh" for regridding

        fields : list of dict
            dictionaries defining the fields with associated WOA data
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # this will be constructed later
        obsFileNames = None

        # call the constructor from the base class (TransectsObservations)
        super(WoaTransectsObservations, self).__init__(
            config, obsFileNames, horizontalResolution,
            transectCollectionName)

        self.fields = fields

    def get_observations(self):
        """
        Read in and set up the observations.

        Returns
        -------
        obsDatasets : OrderedDict
            The observational dataset
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, extract the WOA transects at the requested lat range and lon
        # values
        if self.obsFileNames is None:
            self.extract_observations()

        # then, call the method from the parent class
        return super().get_observations()

    def extract_observations(self):
        """
        Extract WOA23 observations over the range of lat and lon requested
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        longitudes = _get_longitudes(config)

        observationsDirectory = build_obs_path(
            config, 'ocean', 'woa23Subdirectory')

        outObsDirectory = build_config_full_path(
            config=config, section='output',
            relativePathOption='climatologySubdirectory',
            relativePathSection='oceanObservations')

        make_directories(outObsDirectory)

        transectsFileName = \
            f'{outObsDirectory}/{self.transectCollectionName}.nc'
        obsFileNames = OrderedDict()
        for lon in longitudes:
            transectName = f'lon_{lon}'
            obsFileNames[transectName] = transectsFileName

        self.obsFileNames = obsFileNames

        if os.path.exists(transectsFileName):
            return

        print('Preprocessing WOA transect data...')

        minLat = config.getfloat('woaTransects', 'minLat')
        maxLat = config.getfloat('woaTransects', 'maxLat')

        fileName = f'{observationsDirectory}/' \
                   f'woa23_decav_04_pt_s_mon_ann.20241101.nc'

        dsObs = xr.open_dataset(fileName)

        lat = dsObs.lat.values
        mask = np.logical_and(lat >= minLat, lat <= maxLat)
        indices = np.argwhere(mask)
        dsObs = dsObs.isel(lat=slice(indices[0][0],
                                     indices[-1][0]))

        dsObs = dsObs.sel(lon=longitudes, method=str('nearest'))

        dsObs = dsObs.rename({'depth': 'z'})
        attrs = dsObs.z.attrs
        dsObs['z'] = -dsObs.z
        dsObs['z'].attrs = attrs
        dsObs['z'].attrs['positive'] = 'up'

        write_netcdf_with_fill(dsObs, transectsFileName)

        print('  Done.')

    def build_observational_dataset(self, fileName, transectName):
        """
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        transectName : str
            transect name

        Returns
        -------
        dsObs : xarray.Dataset
            The observational dataset
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)

        lon = float(transectName.rsplit('_', 1)[-1])

        dsObs = dsObs.sel(method='nearest', lon=lon)
        lon = dsObs.lon.values

        # must have capital Time
        dsObs = dsObs.rename({'time': 'Time'})
        # make sure month is a coord
        dsObs = dsObs.set_coords('month')
        # add a dummy year to the dataset
        dsObs.coords['year'] = ('Time', np.ones(dsObs.sizes['Time'], int))

        # do some dropping and renaming so we end up with the right coordinates
        # and dimensions for transects
        dsObs = dsObs.rename({'lat': 'nPoints', 'z': 'nz'})
        dsObs['lat'] = dsObs.nPoints
        dsObs['z'] = dsObs.nz
        dsObs['lon'] = ('nPoints', lon * np.ones(dsObs.sizes['nPoints']))
        dsObs = dsObs.drop_vars(['nPoints', 'nz'])

        return dsObs


def _get_longitudes(config):
    longitudes = config.getexpression('woaTransects', 'longitudes',
                                      allow_numpy=True)
    longitudes = np.array(longitudes)
    # make sure longitudes are between -180 and 180
    longitudes = np.sort(np.mod(longitudes + 180., 360.) - 180.)
    return longitudes
