'''
Analysis tasks for comparing Antarctic climatology maps against observations
and reanalysis data.

Authors
-------
Xylar Asay-Davis
'''

import xarray as xr
import numpy

from .climatology_map import ClimatologyMapOcean

from ..shared.constants import constants

from ..shared.io.utility import build_config_full_path

from ..shared.grid import ProjectionGridDescriptor
from ..shared.climatology import get_antarctic_stereographic_projection, \
    MpasClimatology

from ..shared.mpas_xarray import mpas_xarray


class MpasClimatologyMeltAntarctic(MpasClimatology):  # {{{
    """
    A task for creating and remapping climatologies of Antarctic melt rates

    Attributes
    ----------
    landIceMask : xarray.DataArray
        A mask indicating where there is land ice on the ocean grid (thus,
        where melt rates are valid)

    Authors
    -------
    Xylar Asay-Davis
    """

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """

        # first, call setup_and_check from the base class
        # (ClimatologyMapOcean), which will perform some common setup
        super(MpasClimatologyMeltAntarctic, self).setup_and_check()

        landIceFluxMode = self.namelist.get('config_land_ice_flux_mode')
        if landIceFluxMode not in ['standalone', 'coupled']:
            raise ValueError('*** climatologyMapMeltAntarctic requires '
                             'config_land_ice_flux_mode \n'
                             '    to be standalone or coupled.  Otherwise, no '
                             'melt rates are available \n'
                             '    for plotting.')

        # }}}

    def run(self):  # {{{
        """
        Compute climatologies of melt rates from ACME/MPAS output

        This function has been overridden to load ``landIceMask`` from a
        restart file for later use in masking the melt rate.  It then simply
        calls the run function from MpasClimatology

        Authors
        -------
        Xylar Asay-Davis
        """

        # first, load the land-ice mask from the restart file
        dsLandIceMask = xr.open_dataset(self.restartFileName)
        dsLandIceMask = mpas_xarray.subset_variables(dsLandIceMask,
                                                     ['landIceMask'])
        dsLandIceMask = dsLandIceMask.isel(Time=0)
        self.landIceMask = dsLandIceMask.landIceMask > 0.

        # then, call run from the base class (MpasClimatology), which will
        # perform the main function of the task
        super(MpasClimatologyMeltAntarctic, self).run()

        # }}}

    def customize_climatology(self, climatology):  # {{{
        """
        First, calls the parent version of ``customize_climatology`` to set
        up the ``validMask``.  Then, masks the melt rate with the
        ``landIceMask`` and scales the field to be in m/yr instead of
        kg/m^2/s

        Parameters
        ----------
        climatology : ``xarray.Dataset`` object
            the climatology data set

        Returns
        -------
        climatology : ``xarray.Dataset`` object
            the modified climatology data set

        Authors
        -------
        Xylar Asay-Davis
        """

        # perform the custimization in MpasClimatology (adding a 'validMask')
        climatology = \
            super(MpasClimatologyMeltAntarctic, self).customize_climatology(
                    climatology)

        fieldName = self.variableList[0]

        # scale the field to m/yr from kg/m^2/s and mask out non-land-ice areas
        climatology[fieldName] = \
            constants.sec_per_year/constants.rho_fw * \
            climatology[fieldName].where(self.landIceMask)

        return climatology  # }}}

    # }}}


class ClimatologyMapMeltAntarctic(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of antarctic melt rates against
    observations

    Authors
    -------
    Xylar Asay-Davis
    """

    def create_mpas_climatology_task(self, comparisonGridNames, seasons):
        # {{{
        """
        Create an MpasClimatologyMeltAntarctic task to use as a prerequisite
        of this task

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed.

        Authors
        -------
        Xylar Asay-Davis
        """

        climatologyTask = MpasClimatologyMeltAntarctic(
                config=self.config,
                variableList=[self.mpasFieldName],
                taskSuffix='AntarcticMelt',
                componentName=self.componentName,
                comparisonGridNames=comparisonGridNames,
                seasons=seasons,
                tags=['climatology', 'landIceCavities'],
                iselValues=self.iselValues)

        return climatologyTask  # }}}

    def __init__(self, config, season, comparisonGridName):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        season : str
            A season (keys in ``shared.constants.monthDictionary``) to be
            plotted

        comparisonGridName : {'latlon', 'antarctic'}
            The name of the comparison grid that where model data an
            observations are compared

        Authors
        -------
        Xylar Asay-Davis
        """

        self.fieldName = 'melt'
        self.fieldNameInTitle = 'Melt Rate'
        self.season = season
        self.comparisonGridName = comparisonGridName

        self.mpasFieldName = 'timeMonthly_avg_landIceFreshwaterFlux'
        self.iselValues = None

        taskName, tags = self.get_task_name_and_tags()

        tags.extend([self.fieldName, 'landIceCavities'])

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapMeltAntarctic, self).__init__(
                config=config,
                taskName=taskName,
                componentName='ocean',
                tags=tags)
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """

        # first, call setup_and_check from the base class
        # (ClimatologyMapOcean), which will perform some common setup
        super(ClimatologyMapMeltAntarctic, self).setup_and_check()

        landIceFluxMode = self.namelist.get('config_land_ice_flux_mode')
        if landIceFluxMode not in ['standalone', 'coupled']:
            raise ValueError('*** climatologyMapMeltAntarctic requires '
                             'config_land_ice_flux_mode \n'
                             '    to be standalone or coupled.  Otherwise, no '
                             'melt rates are available \n'
                             '    for plotting.')

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            '{}/Rignot_2013_melt_rates_6000.0x6000.0km_10.0km_' \
            'Antarctic_stereo.nc'.format(observationsDirectory)

        self.obsFieldName = 'meltRate'

        # Set appropriate figure labels for melt rates
        self.observationTitleLabel = 'Observations (Rignot et al, 2013)'
        self.outFileLabel = 'meltRignot'
        self.unitsLabel = r'm a$^{-1}$'
        self.diffTitleLabel = 'Model - Observations'
        # }}}

    def build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        dsObs = xr.open_mfdataset(self.obsFileName)

        # create a descriptor of the observation grid using the x/y
        # coordinates
        projection = get_antarctic_stereographic_projection()
        obsDescriptor = ProjectionGridDescriptor.read(
            projection, fileName=self.obsFileName, xVarName='x', yVarName='y')

        return dsObs, obsDescriptor  # }}}

    # }}}


class MpasClimatologyDepthSlices(MpasClimatology):  # {{{
    """
    A task for creating and remapping climatologies of MPAS fields sliced
    at a given set of depths

    Attributes
    ----------
    depths : list of {float, 'bot'}
        A list of depths at which the climatology will be sliced in the
        vertical.

    maxLevelCell : xarray.DataArray
        The vertical index of the bottom cell in MPAS results

    verticalIndices : xarray.DataArray
        The vertical indices of slice to be plotted

    Authors
    -------
    Xylar Asay-Davis
    """

    def run(self):  # {{{
        """
        Compute climatologies of T or S  from ACME/MPAS output

        This function has been overridden to load ``maxLevelCell`` from a
        restart file for later use in indexing bottom T and S.
        ``verticalIndex`` is also computed for later indexing of
        the model level. It then simply calls the run function from
        ClimatologyMapOcean.

        Authors
        -------
        Xylar Asay-Davis
        """

        # first, load the land-ice mask from the restart file
        ds = xr.open_dataset(self.restartFileName)
        ds = mpas_xarray.subset_variables(ds, ['maxLevelCell',
                                               'bottomDepth',
                                               'layerThickness'])

        self.maxLevelCell = ds.maxLevelCell - 1

        verticalIndices = numpy.zeros((len(self.depths), ds.dims['nCells']),
                                      int)

        for depthIndex in range(len(self.depths)):
            depth = self.depths[depthIndex]
            if depth == 'bot':
                # switch to zero-based index
                verticalIndices[depthIndex, :] = ds.maxLevelCell.values - 1
            else:
                nCells = ds.dims['nCells']
                nVertLevels = ds.dims['nVertLevels']
                zInterface = numpy.zeros((nCells, nVertLevels+1))
                zInterface[:, nVertLevels] = -ds.bottomDepth.values
                for zIndex in range(nVertLevels-1, -1, -1):
                    zInterface[:, zIndex] = zInterface[:, zIndex+1] + \
                        ds.layerThickness.values[0, :, zIndex]
                zMid = 0.5*(zInterface[:, 0:-1] + zInterface[:, 1:])

                verticalIndex = numpy.argmin(numpy.abs(zMid-depth), axis=1)
                verticalIndices[depthIndex, :] = verticalIndex

        self.verticalIndices = \
            xr.DataArray.from_dict({'dims': ('depthSlice', 'nCells'),
                                    'data': verticalIndices})

        # then, call run from the base class (MpasClimatology), which will
        # perform the main function of the task
        super(MpasClimatologyDepthSlices, self).run()

    def customize_climatology(self, climatology):  # {{{
        """
        Uses ``verticalIndex`` to slice the 3D climatology field at each
        requested depth.  The resulting field has the depth appended to
        the variable name.

        Parameters
        ----------
        climatology : ``xarray.Dataset`` object
            the climatology data set

        Returns
        -------
        climatology : ``xarray.Dataset`` object
            the modified climatology data set

        Authors
        -------
        Xylar Asay-Davis
        """

        # perform the custimization in MpasClimatology (adding a 'validMask')
        climatology = \
            super(MpasClimatologyDepthSlices, self).customize_climatology(
                    climatology)

        climatology.coords['verticalIndex'] = \
            ('nVertLevels',
             numpy.arange(climatology.dims['nVertLevels']))

        depthNames = [str(depth) for depth in self.depths]

        climatology.coords['depthSlice'] = ('depthSlice', depthNames)

        for variableName in self.variableList:
            if 'nVertLevels' not in climatology[variableName].dims:
                continue

            # mask only the values with the right vertical index
            da = climatology[variableName].where(
                    climatology.verticalIndex == self.verticalIndices)

            # Each vertical layer has at most one non-NaN value so the "sum"
            # over the vertical is used to collapse the array in the vertical
            # dimension
            climatology[variableName] = \
                da.sum(dim='nVertLevels').where(
                    self.verticalIndices <= self.maxLevelCell)

        climatology = climatology.drop('verticalIndex')

        return climatology  # }}}

    # }}}


class ClimatologyMapSose(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of antarctic temperature, salinity, etc.
    against SOSE fields

    Attributes
    ----------

    Authors
    -------
    Xylar Asay-Davis
    """

    @classmethod
    def create_tasks(cls, config):  # {{{
        """
        For each comparison grid, construct one task for computing the
        climatologies and one plotting task for each season.  The climatology
        task is a prerequisite of the plotting tasks, but the plotting tasks
        can run in parallel with one another.

        Parameters
        ----------
        config : MpasAnalysisConfigParser object
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """

        # get the name of the class (something like ClimatologyMapSST)
        className = cls.__name__
        # convert the first letter to lowercase to get the corresponding
        # seciton name in the config file
        sectionName = className[0].lower() + className[1:]

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

        comparisonGridNames = config.getExpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(sectionName))

        # read in what seasons we want to plot
        depths = config.getExpression(sectionName, 'depths')

        if len(depths) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of depths'.format(sectionName))

        first = True
        tasks = []
        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                for depth in depths:
                    # one plotting task for each season, depth and comparison
                    # grid
                    climatologyMapTask = cls(
                            config=config,
                            depth=depth,
                            season=season,
                            comparisonGridName=comparisonGridName)

                    if first:
                        climatologyTask = \
                            climatologyMapTask.create_mpas_climatology_task(
                                    comparisonGridNames=comparisonGridNames,
                                    seasons=seasons,
                                    depths=depths)
                        tasks.append(climatologyTask)

                    climatologyMapTask.add_prerequisite_tasks(
                            [climatologyTask.taskName])

                    climatologyMapTask.climatologyTask = climatologyTask

                    tasks.append(climatologyMapTask)

                    first = False

        return tasks  # }}}

    def __init__(self, config, depth, season, comparisonGridName):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        depth : {float, 'bot'}
            Depth at which to perform the comparison, 'bot' for the ocean floor

        season : str
            A season (keys in ``shared.constants.monthDictionary``) to be
            plotted

        comparisonGridName : {'latlon', 'antarctic'}, optional
            The comparison grid on which to plot the analysis, required if
            mode == 'plot'

        Authors
        -------
        Xylar Asay-Davis
        """

        self.season = season
        self.comparisonGridName = comparisonGridName
        self.depth = depth

        upperField = self.fieldName[0].upper() + self.fieldName[1:]
        if depth == 'bot':
            self.fieldNameInTitle = 'Sea Floor {}'.format(upperField)
        else:
            self.fieldNameInTitle = '{} z={} m'.format(upperField, depth)

        taskName, tags = self.get_task_name_and_tags()
        tags.extend(['SOSE', self.fieldName])

        # append the depth on to the task name
        taskName = '{}_depth_{}'.format(taskName, depth)

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapSose, self).__init__(config=config,
                                                 taskName=taskName,
                                                 componentName='ocean',
                                                 tags=tags)
        # }}}

    def create_mpas_climatology_task(self, comparisonGridNames, seasons,
                                     depths):
        # {{{
        """
        Create an MpasClimatologyMeltAntarctic task to use as a prerequisite
        of this task

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        seasons : list of str
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed.

        depths : list of {float, 'bot'}
            A list of depths where the field is to be sliced

        Authors
        -------
        Xylar Asay-Davis
        """

        upperField = self.fieldName[0].upper() + self.fieldName[1:]
        climatologyTask = MpasClimatologyDepthSlices(
                config=self.config,
                variableList=[self.mpasFieldName],
                taskSuffix='{}DepthSlices'.format(upperField),
                componentName=self.componentName,
                comparisonGridNames=comparisonGridNames,
                seasons=seasons,
                tags=['climatology'],
                iselValues=self.iselValues)

        climatologyTask.depths = depths

        return climatologyTask  # }}}

    def post_process_remapped_mpas_climatology(self, remappedClimatology):
        # {{{
        '''
        Extract the appropriate slice of the remapped climatology for this
        depth.

        Parameters
        ----------
        remappedClimatology : ``xarray.Dataset`` object
            The remapped climatology to be processed

        Returns
        -------
        remappedClimatology : ``xarray.Dataset`` object
            The remapped climatology after processing

        Authors
        -------
        XylarAsay-Davis
        '''

        remappedClimatology = remappedClimatology.sel(
                depthSlice=str(self.depth), drop=True)

        return remappedClimatology

        # }}}

    def build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        dsObs = xr.open_dataset(self.obsFileName)
        dsObs = mpas_xarray.subset_variables(dsObs, [self.obsFieldName,
                                                     'month', 'year'])
        if self.depth != 'bot':
            dsObs = dsObs.sel(method='nearest', depth=self.depth)

        # create a descriptor of the observation grid using the x/y
        # coordinates
        projection = get_antarctic_stereographic_projection()
        obsDescriptor = ProjectionGridDescriptor.read(
            projection, fileName=self.obsFileName, xVarName='x', yVarName='y')

        return dsObs, obsDescriptor  # }}}

    # }}}


class ClimatologyMapSoseTemperature(ClimatologyMapSose):  # {{{
    """
    An analysis task for comparison of antarctic temperature against SOSE
    fields

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, config, depth, season, comparisonGridName):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        depth : {float, 'bot'}
            Depth at which to perform the comparison, 'bot' for the ocean floor

        season : str
            A season (keys in ``shared.constants.monthDictionary``) to be
            plotted

        comparisonGridName : {'latlon', 'antarctic'}, optional
            The comparison grid on which to plot the analysis, required if
            mode == 'plot'

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'temperature'
        self.mpasFieldName = 'timeMonthly_avg_activeTracers_temperature'
        self.iselValues = None

        # call the constructor from the base class (ClimatologyMapSose)
        super(ClimatologyMapSoseTemperature, self).__init__(
                config=config,
                depth=depth,
                season=season,
                comparisonGridName=comparisonGridName)
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class
        # (ClimatologyMapSose), which will perform some common setup
        super(ClimatologyMapSoseTemperature, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations', 'soseSubdirectory')

        self.obsFileName = '{}/SOSE_2005-2010_monthly_pot_temp_6000.0x' \
                           '6000.0km_10.0km_Antarctic_stereo.nc'.format(
                                   observationsDirectory)
        if self.depth == 'bot':
            self.obsFieldName = 'botTheta'
            self.outFileLabel = 'botTemp_SOSE'
        else:
            self.obsFieldName = 'theta'
            self.outFileLabel = 'temp_depth_{}_SOSE'.format(self.depth)
        self.unitsLabel = '$^\circ$C'

        # Set appropriate figure labels for melt rates
        self.observationTitleLabel = 'State Estimate (SOSE)'
        self.diffTitleLabel = 'Model - State Estimate'

        # }}}

    # }}}


class ClimatologyMapSoseSalinity(ClimatologyMapSose):  # {{{
    """
    An analysis task for comparison of antarctic salinity against SOSE
    fields

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, config, depth, season, comparisonGridName):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        depth : {float, 'bot'}
            Depth at which to perform the comparison, 'bot' for the ocean floor

        season : str
            A season (keys in ``shared.constants.monthDictionary``) to be
            plotted

        comparisonGridName : {'latlon', 'antarctic'}, optional
            The comparison grid on which to plot the analysis, required if
            mode == 'plot'

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'salinity'
        self.mpasFieldName = 'timeMonthly_avg_activeTracers_salinity'
        self.iselValues = None

        # call the constructor from the base class (ClimatologyMapSose)
        super(ClimatologyMapSoseSalinity, self).__init__(
                config=config,
                depth=depth,
                season=season,
                comparisonGridName=comparisonGridName)
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class
        # (ClimatologyMapSose), which will perform some common setup
        super(ClimatologyMapSoseSalinity, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations', 'soseSubdirectory')

        self.obsFileName = '{}/SOSE_2005-2010_monthly_salinity_6000.0x' \
                           '6000.0km_10.0km_Antarctic_stereo.nc'.format(
                                   observationsDirectory)
        if self.depth == 'bot':
            self.obsFieldName = 'botSalinity'
            self.outFileLabel = 'botSalin_SOSE'
        else:
            self.obsFieldName = 'salinity'
            self.outFileLabel = 'salin_depth_{}_SOSE'.format(self.depth)
        self.unitsLabel = 'PSU'

        # Set appropriate figure labels for melt rates
        self.observationTitleLabel = 'State Estimate (SOSE)'
        self.diffTitleLabel = 'Model - State Estimate'

        # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
