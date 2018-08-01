# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
"""
An analysis subtasks for plotting comparison of a model transect against
observations or reference data.
"""
# Authors
# -------
# Xylar Asay-Davis, Greg Streletz

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy

from mpas_analysis.shared.plot.plotting import \
    plot_vertical_section, plot_vertical_section_comparison

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.climatology import compute_climatology, \
    get_remapped_mpas_climatology_file_name

from mpas_analysis.shared.constants import constants

from mpas_analysis.ocean.utility import nans_to_numpy_mask


class PlotTransectSubtask(AnalysisTask):  # {{{
    """
    An analysis task for plotting 2D model fields against observations.

    Attributes
    ----------
    season : str
        A season (key in ``shared.constants.monthDictionary``) to be
        plotted.

    transectName : str
        The name of the transect to plot

    remapMpasClimatologySubtask : ``RemapMpasClimatologySubtask``
        The subtask for remapping the MPAS climatology that this subtask
        will plot

    plotObs : bool, optional
        Whether to plot against observations.

    outFileLabel : str
        The prefix on each plot and associated XML file

    fieldNameInTitle : str
        The name of the field being plotted, as used in the plot title

    mpasFieldName : str
        The name of the variable in the MPAS timeSeriesStatsMonthly output

    obsFieldName : str
        The name of the variable to use from the observations file

    observationTitleLabel : str
        the title of the observations subplot

    diffTitleLabel : str, optional
        the title of the difference subplot

    unitsLabel : str
        the units of the plotted field, to be displayed on color bars

    imageCaption : str
        the caption when mousing over the plot or displaying it full
        screen

    galleryGroup : str
        the name of the group of galleries in which this plot belongs

    groupSubtitle : str
        the subtitle of the group in which this plot belongs (or blank
        if none)

    groupLink : str
        a short name (with no spaces) for the link to the gallery group

    galleryName : str
        the name of the gallery in which this plot belongs

    configSectionName : str
        the name of the section where the color map and range is defined
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Greg Streletz

    def __init__(self, parentTask, season, transectName, fieldName,
                 remapMpasClimatologySubtask, plotObs=True,
                 refConfig=None):
        # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (master) task for this subtask

        season : str
            A season (key in ``shared.constants.monthDictionary``) to be
            plotted.

        transectName : str
            The name of the transect to plot

        fieldName : str
            The name of the field to plot (for use in the subtask name only)

        remapMpasClimatologySubtask : ``RemapMpasClimatologySubtask``
            The subtask for remapping the MPAS climatology that this subtask
            will plot

        plotObs : bool, optional
            Whether to plot against observations.

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any), ignored if
            ``plotObs == True``
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.season = season
        self.transectName = transectName
        self.remapMpasClimatologySubtask = remapMpasClimatologySubtask
        self.plotObs = plotObs
        self.refConfig = refConfig
        subtaskName = 'plot{}_{}_{}'.format(season,
                                            transectName.replace(' ', '_'),
                                            fieldName)

        config = parentTask.config
        taskName = parentTask.taskName
        tags = parentTask.tags

        # call the constructor from the base class (AnalysisTask)
        super(PlotTransectSubtask, self).__init__(
                config=config, taskName=taskName, subtaskName=subtaskName,
                componentName='ocean', tags=tags)

        # this task should not run until the remapping subtasks are done, since
        # it relies on data from those subtasks
        self.run_after(remapMpasClimatologySubtask)
        # }}}

    def set_plot_info(self, outFileLabel, fieldNameInTitle, mpasFieldName,
                      refFieldName, refTitleLabel, unitsLabel,
                      imageCaption, galleryGroup, groupSubtitle, groupLink,
                      galleryName, configSectionName,
                      diffTitleLabel='Model - Observations'):
        # {{{
        '''
        Store attributes related to plots, plot file names and HTML output.

        Parameters
        ----------
        outFileLabel : str
            The prefix on each plot and associated XML file

        fieldNameInTitle : str
            The name of the field being plotted, as used in the plot title

        mpasFieldName : str
            The name of the variable in the MPAS timeSeriesStatsMonthly output

        refFieldName : str
            The name of the variable to use from the observations or reference
            file

        refTitleLabel : str
            the title of the observations or reference subplot

        unitsLabel : str
            the units of the plotted field, to be displayed on color bars

        imageCaption : str
            the caption when mousing over the plot or displaying it full
            screen

        galleryGroup : str
            the name of the group of galleries in which this plot belongs

        groupSubtitle : str
            the subtitle of the group in which this plot belongs (or blank
            if none)

        groupLink : str
            a short name (with no spaces) for the link to the gallery group

        galleryName : str
            the name of the gallery in which this plot belongs

        configSectionName : str
            the name of the section where the color map and range is defined

        diffTitleLabel : str, optional
            the title of the difference subplot
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.outFileLabel = outFileLabel
        self.fieldNameInTitle = fieldNameInTitle
        self.mpasFieldName = mpasFieldName
        self.refFieldName = refFieldName
        self.refTitleLabel = refTitleLabel
        self.diffTitleLabel = diffTitleLabel
        self.unitsLabel = unitsLabel

        # xml/html related variables
        self.imageCaption = imageCaption
        self.galleryGroup = galleryGroup
        self.groupSubtitle = groupSubtitle
        self.groupLink = groupLink
        self.galleryName = galleryName

        self.fieldNameInTitle = fieldNameInTitle
        self.thumbnailDescription = self.season

        self.configSectionName = configSectionName
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(PlotTransectSubtask, self).setup_and_check()

        config = self.config
        self.startYear = config.getint('climatology', 'startYear')
        self.endYear = config.getint('climatology', 'endYear')
        self.startDate = config.get('climatology', 'startDate')
        self.endDate = config.get('climatology', 'endDate')

        mainRunName = config.get('runs', 'mainRunName')

        self.xmlFileNames = []

        prefixPieces = []
        if self.outFileLabel is not '':
            prefixPieces.append(self.outFileLabel)
        prefixPieces.append(self.transectName.replace(' ', '_'))
        prefixPieces.append(mainRunName)
        years = 'years{:04d}-{:04d}'.format(self.startYear, self.endYear)
        prefixPieces.extend([self.season, years])

        self.filePrefix = '_'.join(prefixPieces)

        self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                    self.filePrefix))
        # }}}

    def run_task(self):  # {{{
        """
        Plots a comparison of ACME/MPAS output to SST, MLD or SSS observations
        or a reference run
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Greg Streletz

        season = self.season
        transectName = self.transectName
        self.logger.info("\nPlotting {} over season {}"
                         "...".format(self.fieldNameInTitle,
                                      season))

        # first read the model climatology
        remappedFileName = \
            self.remapMpasClimatologySubtask.get_remapped_file_name(
                    season=season, comparisonGridName=transectName)

        remappedModelClimatology = xr.open_dataset(remappedFileName)

        # now the observations or reference run
        if self.plotObs:
            verticalComparisonGridName = \
                self.remapMpasClimatologySubtask.verticalComparisonGridName
            remappedFileName = \
                self.remapMpasClimatologySubtask.obsDatasets.get_out_file_name(
                    transectName,
                    verticalComparisonGridName)
            remappedRefClimatology = xr.open_dataset(remappedFileName)

            # if Time is an axis, take the appropriate avarage to get the
            # climatology
            if 'Time' in remappedRefClimatology.dims:
                monthValues = constants.monthDictionary[season]
                remappedRefClimatology = compute_climatology(
                        remappedRefClimatology, monthValues, maskVaries=True)

        elif self.refConfig is not None:
            climatologyName = self.remapMpasClimatologySubtask.climatologyName
            remappedFileName = \
                get_remapped_mpas_climatology_file_name(
                    self.refConfig, season=season,
                    componentName=self.componentName,
                    climatologyName=climatologyName,
                    comparisonGridName=transectName)
            remappedRefClimatology = xr.open_dataset(remappedFileName)
            refStartYear = self.refConfig.getint('climatology', 'startYear')
            refEndYear = self.refConfig.getint('climatology', 'endYear')
            if refStartYear != self.startYear or refEndYear != self.endYear:
                self.refTitleLabel = '{}\n(years {:04d}-{:04d})'.format(
                        self.refTitleLabel, refStartYear, refEndYear)

        else:
            remappedRefClimatology = None

        self._plot_transect(remappedModelClimatology, remappedRefClimatology)

        remappedModelClimatology.close()

        if remappedRefClimatology is not None:
            remappedRefClimatology.close()

        # }}}

    def _plot_transect(self, remappedModelClimatology, remappedRefClimatology):
        # {{{
        """ plotting the transect """

        season = self.season
        config = self.config
        configSectionName = self.configSectionName

        mainRunName = config.get('runs', 'mainRunName')

        # broadcast x and z to have the same dimensions
        x, z = xr.broadcast(remappedModelClimatology.x,
                            remappedModelClimatology.z)

        # set lat and lon in case we want to plot versus these quantities
        lat = remappedModelClimatology.lat
        lon = remappedModelClimatology.lon

        # convert x, z, lat, and lon to numpy arrays; make a copy because
        # they are sometimes read-only (not sure why)
        x = x.values.copy().transpose()
        z = z.values.copy().transpose()
        lat = lat.values.copy().transpose()
        lon = lon.values.copy().transpose()
        self.lat = lat
        self.lon = lon

        # z is masked out with NaNs in some locations (where there is land) but
        # this makes pcolormesh unhappy so we'll zero out those locations
        z[numpy.isnan(z)] = 0.

        modelOutput = nans_to_numpy_mask(
            remappedModelClimatology[self.mpasFieldName].values)
        modelOutput = modelOutput.transpose()

        if remappedRefClimatology is None:
            refOutput = None
            bias = None
        else:
            refOutput = remappedRefClimatology[self.refFieldName]
            dims = refOutput.dims
            refOutput = nans_to_numpy_mask(refOutput.values)
            if dims[1] != 'nPoints':
                assert(dims[0] == 'nPoints')
                refOutput = refOutput.transpose()

            bias = modelOutput - refOutput

        filePrefix = self.filePrefix
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        title = '{}\n({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, season, self.startYear,
                self.endYear)


        xLabel = 'Distance [km]'
        yLabel = 'Depth [m]'


        # define the axis labels and the data to use for the upper
        # x axis or axes, if such additional axes have been requested

        upperXAxes = config.get('transects', 'upperXAxes')
        numUpperTicks = config.getint('transects', 'numUpperTicks')
        upperXAxisTickLabelPrecision = config.getint(
            'transects',
            'upperXAxisTickLabelPrecision')

        self._set_third_x_axis_to_none()

        if upperXAxes == 'neither':
            self._set_second_x_axis_to_none()
        elif upperXAxes == 'lat':
            self._set_second_x_axis_to_latitude()
        elif upperXAxes == 'lon':
            self._set_second_x_axis_to_longitude()
        elif upperXAxes == 'both':
            self._set_second_x_axis_to_longitude()
            self._set_third_x_axis_to_latitude()
        elif upperXAxes == 'greatestExtent':
            if self._greatest_extent(lat, lon):
                self._set_second_x_axis_to_latitude()
            else:
                self._set_second_x_axis_to_longitude()
        elif upperXAxes == 'strictlyMonotonic':
            if self._strictly_monotonic(lat, lon):
                self._set_second_x_axis_to_latitude()
            else:
                self._set_second_x_axis_to_longitude()
        elif upperXAxes == 'mostMonotonic':
            if self._most_monotonic(lat, lon):
                self._set_second_x_axis_to_latitude()
            else:
                self._set_second_x_axis_to_longitude()
        elif upperXAxes == 'mostStepsInSameDirection':
            if self._most_steps_in_same_direction(lat, lon):
                self._set_second_x_axis_to_latitude()
            else:
                self._set_second_x_axis_to_longitude()
        elif upperXAxes == 'fewestDirectionChanges':
            if self._fewest_direction_changes(lat, lon):
                self._set_second_x_axis_to_latitude()
            else:
                self._set_second_x_axis_to_longitude()
        else:
            raise ValueError('invalid option for upperXAxes')


        # get the parameters determining what type of plot to use,
        # what line styles and line colors to use, and whether and how
        # to label contours

        compareAsContours = config.getboolean('transects',
                                              'compareAsContoursOnSinglePlot')

        contourLineStyle = config.get('transects', 'contourLineStyle')
        contourLineColor = config.get('transects', 'contourLineColor')
        comparisonContourLineStyle = config.get('transects',
                                                'comparisonContourLineStyle')
        comparisonContourLineColor = config.get('transects',
                                                'comparisonContourLineColor')

        if compareAsContours:
            labelContours = config.getboolean('transects',
                                    'labelContoursOnContourComparisonPlots')
        else:
            labelContours = config.getboolean('transects',
                                              'labelContoursOnHeatmaps')
        
        contourLabelPrecision = config.getint('transects',
                                              'contourLabelPrecision')


        # construct a three-panel comparison plot for the transect, or a
        # single-panel contour comparison plot if compareAsContours is True
        
        plot_vertical_section_comparison(
            config,
            x,
            z,
            modelOutput,
            refOutput,
            bias,
            outFileName,
            configSectionName,
            cbarLabel=self.unitsLabel,
            xlabel=xLabel,
            ylabel=yLabel,
            title=title,
            modelTitle='{}'.format(mainRunName),
            refTitle=self.refTitleLabel,
            diffTitle=self.diffTitleLabel,
            secondXAxisData=self.secondXAxisData,
            secondXAxisLabel=self.secondXAxisLabel,
            thirdXAxisData=self.thirdXAxisData,
            thirdXAxisLabel=self.thirdXAxisLabel,
            numUpperTicks=numUpperTicks,
            upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
            invertYAxis=False,
            backgroundColor='#918167',
            compareAsContours=compareAsContours,
            lineStyle=contourLineStyle,
            lineColor=contourLineColor,
            comparisonContourLineStyle=comparisonContourLineStyle,
            comparisonContourLineColor=comparisonContourLineColor,
            labelContours=labelContours,
            contourLabelPrecision=contourLabelPrecision
            )

        caption = '{} {}'.format(season, self.imageCaption)
        write_image_xml(
            config,
            filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup=self.galleryGroup,
            groupSubtitle=self.groupSubtitle,
            groupLink=self.groupLink,
            gallery=self.galleryName,
            thumbnailDescription=self.thumbnailDescription,
            imageDescription=caption,
            imageCaption=caption)

        # }}}

    def _set_second_x_axis_to_latitude(self):
        self.secondXAxisData = self.lat
        self.secondXAxisLabel = 'Latitude'

    def _set_second_x_axis_to_longitude(self):
        self.secondXAxisData = self.lon
        self.secondXAxisLabel = 'Longitude'

    def _set_second_x_axis_to_none(self):
        self.secondXAxisData = None
        self.secondXAxisLabel = None

    def _set_third_x_axis_to_latitude(self):
        self.thirdXAxisData = self.lat
        self.thirdXAxisLabel = 'Latitude'

    def _set_third_x_axis_to_longitude(self):
        self.thirdXAxisData = self.lon
        self.thirdXAxisLabel = 'Longitude'

    def _set_third_x_axis_to_none(self):
        self.thirdXAxisData = None
        self.thirdXAxisLabel = None

    def _greatest_extent(self, lat, lon):
        # {{{
        """
        Returns True if lat has a greater extent (in degrees) than lon (or the
        same extent as lon), and False otherwise.

        Authors
        -------
        Greg Streletz
        """
        lat_extent = numpy.max(lat) - numpy.min(lat)

        lon_diff = numpy.diff(lon)
        lon_jump = numpy.where(
            numpy.logical_or(lon_diff > 180, lon_diff < -180), True, False)

        if numpy.all(lon_jump == False):  # transect does not cross
                                          # periodic boundary

            lon_extent = numpy.max(lon) - numpy.min(lon)

        else:  # transect crosses periodic boundary at least once, so min and
               # max cannot be used to calculate the longitudinal extent

            # calculate the indices at which the transect crosses the right-hand
            # periodic boundary and enters the leftmost region of the domain
            lon_r = numpy.sort(numpy.nonzero(lon_diff < -180)[0] + 1)

            # calculate the indices at which the transect crosses the left-hand
            # periodic boundary and enters the rightmost region of the domain
            lon_l = numpy.sort(numpy.nonzero(lon_diff > 180)[0] + 1)

            mins = []
            maxes = []
            last_idx = 0

            while(len(lon_r) > 0 and len(lon_l) > 0):
                if lon_r[0] < lon_l[0]:
                    mins.append(numpy.min(lon[last_idx:lon_r[0]]))
                    last_idx = lon_r[0]
                    lon_r = lon_r[1:]
                else:
                    maxes.append(numpy.max(lon[last_idx:lon_l[0]]))
                    last_idx = lon_l[0]
                    lon_l = lon_l[1:]

            min = numpy.min(mins)
            max = numpy.max(maxes)

            if min > max:

                lon_extent = 360.0 - (min - max)

                # The transect has crossed the periodic boundary,
                # so unless the transect has covered the entire
                # longitudinal range (360 degrees), the "leftmost"
                # or "min" coordinate of the transect will have a
                # greater value (in degrees longitude) that the
                # "rightmost" or "max" coordinate of the transect,
                # and (min - max) is the size (in degrees longitude)
                # of the region *not* accessed by the transect.  So,
                # the longitudinal extent of the transect is 360
                # degrees minus this value.

            else:   # if max >= min, the transect has covered the
                    # entire longitudinal range (360 degrees)

                lon_extent = 360.0

        if lat_extent >= lon_extent:
            return True
        else:
            return False
        # }}}

    def _strictly_monotonic(self, lat, lon):
        # {{{
        """
        Returns True if lat is strictly monotonic;  returns false if lon is
        strictly monotonic and lat is not strictly monotonic;  throws an error
        if neither are strictly monotonic.

        Authors
        -------
        Greg Streletz
        """
        lat_diff = numpy.diff(lat)
        if numpy.all(lat_diff > 0) or numpy.all(lat_diff < 0):
            return True
        lon_diff = numpy.diff(lon)
        lon_diff = numpy.where(lon_diff > 180, lon_diff - 360, lon_diff)
        lon_diff = numpy.where(lon_diff < -180, lon_diff + 360, lon_diff)
        if numpy.all(lon_diff > 0) or numpy.all(lon_diff < 0):
            return False
        else:
            raise ValueError('neither input array is strictly monotonic')
        # }}}

    def _most_monotonic(self, lat, lon):
        # {{{
        """
        Returns True if lat is "more monotonic" than lon in terms of the
        difference between the total number of degrees traversed and the net
        number of degrees traversed (or if both are equally as monotonic in
        this sense), and False otherwise.

        Authors
        -------
        Greg Streletz
        """
        lat_diff = numpy.diff(lat)
        lat_score = numpy.sum(numpy.fabs(lat_diff)) - abs(numpy.sum(lat_diff))

        lon_diff = numpy.diff(lon)
        lon_diff = numpy.where(lon_diff > 180, lon_diff - 360, lon_diff)
        lon_diff = numpy.where(lon_diff < -180, lon_diff + 360, lon_diff)
        lon_score = numpy.sum(numpy.fabs(lon_diff)) - abs(numpy.sum(lon_diff))

        if lat_score <= lon_score:
            return True
        else:
            return False
        # }}}

    def _most_steps_in_same_direction(self, lat, lon):
        # {{{
        """
        Returns True if lat is has more steps in the same direction (either
        steps that increase the latitude or steps that decrease the latitude)
        than lon (or the same number as lon), and False otherwise.

        Authors
        -------
        Greg Streletz
        """
        lat_changes = numpy.diff(lat)
        lat_changes = lat_changes[lat_changes != 0.0]  # ignore flat regions
        lat_changedirs = lat_changes / numpy.fabs(lat_changes)
        num_lat_positive = numpy.sum(numpy.where(lat_changedirs == 1, 1, 0))
        num_lat_negative = numpy.sum(numpy.where(lat_changedirs == -1, 1, 0))
        max_steps_lat = max(num_lat_positive, num_lat_negative)

        lon_changes = numpy.diff(lon)
        lon_changes = \
            numpy.where(lon_changes > 180, lon_changes - 360, lon_changes)
        lon_changes = \
            numpy.where(lon_changes < -180, lon_changes + 360, lon_changes)
        lon_changes = lon_changes[lon_changes != 0.0]  # ignore flat regions
        lon_changedirs = lon_changes / numpy.fabs(lon_changes)
        num_lon_positive = numpy.sum(numpy.where(lon_changedirs == 1, 1, 0))
        num_lon_negative = numpy.sum(numpy.where(lon_changedirs == -1, 1, 0))
        max_steps_lon = max(num_lon_positive, num_lon_negative)

        if max_steps_lat >= max_steps_lon:
            return True
        else:
            return False
        # }}}

    def _fewest_direction_changes(self, lat, lon):
        # {{{
        """
        Returns True if lat is has fewer changes in direction (from increasing
        in value to decreasing in value, or vice versa) than lon (or if both
        have to same number of changes), and False otherwise.

        Authors
        -------
        Greg Streletz
        """
        lat_changes = numpy.diff(lat)
        lat_changes = lat_changes[lat_changes != 0.0]  # ignore flat regions
        lat_changedirs = lat_changes / numpy.fabs(lat_changes)
        lat_changedir_changes = numpy.diff(lat_changedirs)
        num_lat_dirchanges = \
            len(lat_changedir_changes[lat_changedir_changes != 0.0])

        lon_changes = numpy.diff(lon)
        lon_changes = \
            numpy.where(lon_changes > 180, lon_changes - 360, lon_changes)
        lon_changes = \
            numpy.where(lon_changes < -180, lon_changes + 360, lon_changes)
        lon_changes = lon_changes[lon_changes != 0.0]  # ignore flat regions
        lon_changedirs = lon_changes / numpy.fabs(lon_changes)
        lon_changedir_changes = numpy.diff(lon_changedirs)
        num_lon_dirchanges = \
            len(lon_changedir_changes[lon_changedir_changes != 0.0])

        if num_lat_dirchanges <= num_lon_dirchanges:
            return True
        else:
            return False
        # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
