'''
Analysis tasks for comparing Antarctic climatology maps against observations
and reanalysis data.

Authors
-------
Xylar Asay-Davis
'''

import xarray as xr

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
        Compute climatologie sof melt rates from ACME/MPAS output

        This function has been overridden to load ``landIceMask`` from a
        restart file for later use in masking the melt rate.  It then simply
        calls the run function from

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

        # then, call run from the base class (ClimatologyMapOcean), which will
        # perform the main function of the task
        super(MpasClimatologyMeltAntarctic, self).run()

        # }}}

    def customize_climatology(self, climatology):  # {{{
        """
        A function that can be overridden by child classes for purposes of
        making custom changes to the climatology data set after slicing and
        masking and before remapping.  By default, a field 'validMask' is added
        to the climatology.  After remapping, this field indicates which cells
        on the remapped grid came from valid cells on the MPAS grid, useful for
        plotting a land mask (for example).

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

    def create_mpas_climatology(self, comparisonGridNames, seasons):  # {{{
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
        if landIceFluxMode not in  ['standalone', 'coupled']:
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

        # }}}

    def _build_observational_dataset(self):  # {{{
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

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
