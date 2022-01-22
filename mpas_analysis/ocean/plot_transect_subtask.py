# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
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

from geometric_features import FeatureCollection

from mpas_tools.ocean.transects import get_outline_segments

from mpas_analysis.shared.plot import plot_vertical_section_comparison, \
    savefig, add_inset

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.climatology import compute_climatology, \
    get_remapped_mpas_climatology_file_name

from mpas_analysis.shared.constants import constants


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

    computeTransectsSubtask : ``ComputeTransectsSubtask``
        The subtask for computing the MPAS climatology that this subtask
        will plot

    plotObs : bool, optional
        Whether to plot against observations.

    controlConfig :  ``MpasAnalysisConfigParser``
        Configuration options for a control run (if any), ignored if
        ``plotObs == True``

    horizontalBounds : [float, float]
        The distance in km of the left- and right-hand edges of the plot,
        chosen automatically by default or if the list is empty

    outFileLabel : str
        The prefix on each plot and associated XML file

    fieldNameInTitle : str
        The name of the field being plotted, as used in the plot title

    mpasFieldName : str
        The name of the variable in the MPAS timeSeriesStatsMonthly output

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
                 computeTransectsSubtask, plotObs=True,
                 controlConfig=None, horizontalBounds=None):
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

        computeTransectsSubtask : ``ComputeTransectsSubtask``
            The subtask for computing the MPAS climatology that this subtask
            will plot

        plotObs : bool, optional
            Whether to plot against observations.

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any), ignored if
            ``plotObs == True``

        horizontalBounds : [float, float], optional
            The distance in km of the left- and right-hand edges of the plot,
            chosen automatically by default or if the list is empty
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.season = season
        self.transectName = transectName
        self.computeTransectsSubtask = computeTransectsSubtask
        self.plotObs = plotObs
        self.controlConfig = controlConfig
        if horizontalBounds is not None and len(horizontalBounds) == 2:
            self.horizontalBounds = horizontalBounds
        else:
            self.horizontalBounds = None
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
        self.run_after(computeTransectsSubtask)
        # }}}

    def set_plot_info(self, outFileLabel, fieldNameInTitle, mpasFieldName,
                      refFieldName, refTitleLabel, unitsLabel,
                      imageCaption, galleryGroup, groupSubtitle, groupLink,
                      galleryName, configSectionName, verticalBounds,
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

        verticalBounds : list
            the min and max for the vertical axis, or an emtpy list if the
            range automatically determined by matplotlib should be used

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
        if len(verticalBounds) == 0:
            self.verticalBounds = None
        else:
            self.verticalBounds = verticalBounds
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
        if self.outFileLabel != '':
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
        or a control run
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
            self.computeTransectsSubtask.get_remapped_file_name(
                season=season, comparisonGridName=transectName)

        remappedModelClimatology = xr.open_dataset(remappedFileName)

        # now the observations or control run
        if self.plotObs:
            verticalComparisonGridName = \
                self.computeTransectsSubtask.verticalComparisonGridName
            remappedFileName = \
                self.computeTransectsSubtask.obsDatasets.get_out_file_name(
                    transectName,
                    verticalComparisonGridName)
            remappedRefClimatology = xr.open_dataset(remappedFileName)

            # if Time is an axis, take the appropriate average to get the
            # climatology
            if 'Time' in remappedRefClimatology.dims:
                monthValues = constants.monthDictionary[season]
                remappedRefClimatology = compute_climatology(
                    remappedRefClimatology, monthValues, maskVaries=True)

        elif self.controlConfig is not None:
            climatologyName = self.computeTransectsSubtask.climatologyName
            remappedFileName = \
                get_remapped_mpas_climatology_file_name(
                    self.controlConfig, season=season,
                    componentName=self.componentName,
                    climatologyName=climatologyName,
                    comparisonGridName=transectName)
            remappedRefClimatology = xr.open_dataset(remappedFileName)
            controlStartYear = self.controlConfig.getint('climatology',
                                                         'startYear')
            controlEndYear = self.controlConfig.getint('climatology',
                                                       'endYear')
            if controlStartYear != self.startYear or \
                    controlEndYear != self.endYear:
                self.refTitleLabel = '{}\n(years {:04d}-{:04d})'.format(
                    self.refTitleLabel, controlStartYear, controlEndYear)

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

        remap = self.computeTransectsSubtask.remap

        if remap:
            x = 1e-3*remappedModelClimatology.x
            z = remappedModelClimatology.z

            # set lat and lon in case we want to plot versus these quantities
            lat = remappedModelClimatology.lat
            lon = remappedModelClimatology.lon

            if len(lat.dims) > 1:
                lat = lat[:, 0]

            if len(lon.dims) > 1:
                lon = lon[:, 0]

            # z is masked out with NaNs in some locations (where there is land)
            # but this makes pcolormesh unhappy so we'll zero out those
            # locations
            z = z.where(z.notnull(), 0.)

        else:
            x = 1e-3*remappedModelClimatology.dNode
            z = None
            lon = remappedModelClimatology.lonNode
            lat = remappedModelClimatology.latNode

            remappedModelClimatology['dNode'] = x

            # flatten the x, lon and lat arrays because this is what
            # vertical_section is expecting
            x = xr.DataArray(data=x.values.ravel(), dims=('nx',))
            lon = xr.DataArray(data=lon.values.ravel(), dims=('nx',))
            lat = xr.DataArray(data=lat.values.ravel(), dims=('nx',))

        # This will do strange things at the antemeridian but there's little
        # we can do about that.
        lon_pm180 = numpy.mod(lon + 180., 360.) - 180.

        if self.horizontalBounds is not None:
            mask = numpy.logical_and(
                x.values >= self.horizontalBounds[0],
                x.values <= self.horizontalBounds[1])
            inset_lon = lon_pm180.values[mask]
            inset_lat = lat.values[mask]
        else:
            inset_lon = lon_pm180.values
            inset_lat = lat.values

        fc = FeatureCollection()
        fc.add_feature(
            {"type": "Feature",
             "properties": {"name": self.transectName,
                            "author": 'Xylar Asay-Davis',
                            "object": 'transect',
                            "component": 'ocean',
                            "tags": ''},
             "geometry": {
                 "type": "LineString",
                 "coordinates": list(map(list, zip(inset_lon, inset_lat)))}})

        modelOutput = remappedModelClimatology[self.mpasFieldName]

        if remap:
            triangulation_args = None
            dOutlineModel = None
            zOutlineModel = None
        else:
            triangulation_args = self._get_ds_triangulation(
                remappedModelClimatology)
            dOutlineModel, zOutlineModel = \
                get_outline_segments(remappedModelClimatology)

        if self.controlConfig is None:
            dOutlineRef = None
            zOutlineRef = None
            dOutlineDiff = None
            zOutlineDiff = None
        else:
            dOutlineRef = dOutlineModel
            zOutlineRef = zOutlineModel
            dOutlineDiff = dOutlineModel
            zOutlineDiff = zOutlineModel

        if remappedRefClimatology is None:
            refOutput = None
            bias = None
        else:
            # todo: add a check to make sure the reference is on the same MPAS
            #   mesh
            refOutput = remappedRefClimatology[self.refFieldName]
            bias = modelOutput - refOutput

        filePrefix = self.filePrefix
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        title = '{}\n({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, season, self.startYear,
                self.endYear)

        xs = [x]
        xLabels = ['Distance (km)']
        yLabel = 'Depth (m)'

        # define the axis labels and the data to use for the upper
        # x axis or axes, if such additional axes have been requested

        upperXAxes = config.get('transects', 'upperXAxes')
        numUpperTicks = config.getint('transects', 'numUpperTicks')
        upperXAxisTickLabelPrecision = config.getint(
            'transects',
            'upperXAxisTickLabelPrecision')

        add_lat = False
        add_lon = False

        if upperXAxes in ['lat', 'both']:
            add_lat = True
        if upperXAxes in ['lon', 'both']:
            add_lon = True
        if upperXAxes == 'greatestExtent':
            if self._lat_greater_extent(lat, lon):
                add_lat = True
            else:
                add_lon = True
        elif upperXAxes == 'strictlyMonotonic':
            if self._strictly_monotonic(lat):
                add_lat = True
            elif self._strictly_monotonic(lon):
                add_lon = True
            else:
                raise ValueError('Neither lat nor lon is strictly monotonic.')
        elif upperXAxes == 'mostMonotonic':
            if self._lat_most_monotonic(lat, lon):
                add_lat = True
            else:
                add_lon = True
        elif upperXAxes == 'mostStepsInSameDirection':
            if self._lat_most_steps_in_same_direction(lat, lon):
                add_lat = True
            else:
                add_lon = True
        elif upperXAxes == 'fewestDirectionChanges':
            if self._lat_fewest_direction_changes(lat, lon):
                add_lat = True
            else:
                add_lon = True
        else:
            raise ValueError('invalid option for upperXAxes')

        if add_lat:
            xs.append(lat)
            xLabels.append('Latitude')
        if add_lon:
            xs.append(lon)
            xLabels.append('Longitude')

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
            labelContours = config.getboolean(
                'transects', 'labelContoursOnContourComparisonPlots')
        else:
            labelContours = config.getboolean('transects',
                                              'labelContoursOnHeatmaps')

        contourLabelPrecision = config.getint('transects',
                                              'contourLabelPrecision')

        if config.has_option(configSectionName, 'titleFontSize'):
            titleFontSize = config.getint(configSectionName, 'titleFontSize')
        else:
            titleFontSize = None

        if config.has_option(configSectionName, 'defaultFontSize'):
            defaultFontSize = config.getint(configSectionName,
                                            'defaultFontSize')
        else:
            defaultFontSize = None

        if config.has_option(configSectionName, 'axisFontSize'):
            axisFontSize = config.getint(configSectionName, 'axisFontSize')
        else:
            axisFontSize = None

        # construct a three-panel comparison plot for the transect, or a
        # single-panel contour comparison plot if compareAsContours is True

        fig, axes, suptitle = plot_vertical_section_comparison(
            config,
            modelOutput,
            refOutput,
            bias,
            configSectionName,
            xCoords=xs,
            zCoord=z,
            triangulation_args=triangulation_args,
            xOutlineModel=dOutlineModel,
            zOutlineModel=zOutlineModel,
            xOutlineRef=dOutlineRef,
            zOutlineRef=zOutlineRef,
            xOutlineDiff=dOutlineDiff,
            zOutlineDiff=zOutlineDiff,
            colorbarLabel=self.unitsLabel,
            xlabels=xLabels,
            ylabel=yLabel,
            title=title,
            modelTitle='{}'.format(mainRunName),
            refTitle=self.refTitleLabel,
            diffTitle=self.diffTitleLabel,
            numUpperTicks=numUpperTicks,
            upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
            invertYAxis=False,
            backgroundColor='#d9bf96',
            invalidColor='#d9bf96',
            xLim=self.horizontalBounds,
            yLim=self.verticalBounds,
            compareAsContours=compareAsContours,
            lineStyle=contourLineStyle,
            lineColor=contourLineColor,
            comparisonContourLineStyle=comparisonContourLineStyle,
            comparisonContourLineColor=comparisonContourLineColor,
            labelContours=labelContours,
            contourLabelPrecision=contourLabelPrecision,
            plotTitleFontSize=titleFontSize,
            defaultFontSize=defaultFontSize,
            axisFontSize=axisFontSize)

        # shift the super-title a little to the left to make room for the inset
        pos = suptitle.get_position()
        suptitle.set_position((pos[0] - 0.05, pos[1]))

        if not remap:
            # add open ocean or ice shelves
            d = remappedModelClimatology.dNode.values.ravel()
            ssh = remappedModelClimatology.ssh.values.ravel()
            if 'landIceFraction' in remappedModelClimatology:
                # plot ice in light blue
                color = '#e1eaf7'
            else:
                # plot open ocean in white
                color = 'white'
            mask = ssh < 0.
            for ax in axes:
                ax.fill_between(d, ssh, numpy.zeros(ssh.shape), where=mask,
                                interpolate=True, color=color,
                                edgecolor='black', linewidth=1.)

        # make a red start axis and green end axis to correspond to the dots
        # in the inset
        for ax in axes:
            ax.spines['left'].set_color('red')
            ax.spines['right'].set_color('green')
            ax.spines['left'].set_linewidth(4)
            ax.spines['right'].set_linewidth(4)

        add_inset(fig, fc, width=1.5, height=1.5, xbuffer=0.1, ybuffer=0.1)

        savefig(outFileName, config, tight=False)

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

    def _lat_greater_extent(self, lat, lon):
        # {{{
        """
        Returns True if lat has a greater extent (in degrees) than lon (or the
        same extent as lon), and False otherwise.

        Authors
        -------
        Greg Streletz, Xylar Asay-Davis
        """
        lat_extent = numpy.max(lat) - numpy.min(lat)

        lon_diff = numpy.diff(lon)
        lon_jump = numpy.where(
            numpy.logical_or(lon_diff > 180, lon_diff < -180), True, False)

        if numpy.all(numpy.logical_not(lon_jump)):
            # transect does not cross periodic boundary

            lon_extent = numpy.max(lon) - numpy.min(lon)

        else:
            # transect crosses periodic boundary at least once, so min and
            # max cannot be used to calculate the longitudinal extent

            # calculate the indices at which the transect crosses the
            # right-hand periodic boundary and enters the leftmost region of
            # the domain
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

            else:
                # if max >= min, the transect has covered the
                # entire longitudinal range (360 degrees)

                lon_extent = 360.0

        if lat_extent >= lon_extent:
            return True
        else:
            return False
        # }}}

    def _strictly_monotonic(self, coord):
        # {{{
        """Whether the coordinate is strictly monotonic"""
        # Authors
        # -------
        # Greg Streletz, Xylar Asay-Davis

        coord_diff = numpy.diff(coord.values)
        coord_diff = numpy.where(coord_diff > 180, coord_diff - 360, coord_diff)
        coord_diff = numpy.where(coord_diff < -180, coord_diff + 360,
                                 coord_diff)
        return numpy.all(coord_diff > 0) or numpy.all(coord_diff < 0)
        # }}}

    def _lat_most_monotonic(self, lat, lon):
        # {{{
        """
        Returns True if lat is "more monotonic" than lon in terms of the
        difference between the total number of degrees traversed and the net
        number of degrees traversed (or if both are equally as monotonic in
        this sense), and False otherwise.
        """
        # Authors
        # -------
        # Greg Streletz, Xylar Asay-Davis
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

    def _lat_most_steps_in_same_direction(self, lat, lon):
        # {{{
        """
        Returns True if lat is has more steps in the same direction (either
        steps that increase the latitude or steps that decrease the latitude)
        than lon (or the same number as lon), and False otherwise.
        """
        # Authors
        # -------
        # Greg Streletz, Xylar Asay-Davis

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

    def _lat_fewest_direction_changes(self, lat, lon):
        # {{{
        """
        Returns True if lat is has fewer changes in direction (from increasing
        in value to decreasing in value, or vice versa) than lon (or if both
        have to same number of changes), and False otherwise.

        Authors
        -------
        Greg Streletz, Xylar Asay-Davis
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

    def _get_ds_triangulation(self, dsTransectTriangles):
        """get matplotlib Triangulation from triangulation dataset"""

        nTransectTriangles = dsTransectTriangles.sizes['nTransectTriangles']
        dNode = dsTransectTriangles.dNode.isel(
            nSegments=dsTransectTriangles.segmentIndices,
            nHorizBounds=dsTransectTriangles.nodeHorizBoundsIndices)
        x = dNode.values.ravel()

        zTransectNode = dsTransectTriangles.zTransectNode
        y = zTransectNode.values.ravel()

        tris = numpy.arange(3 * nTransectTriangles).reshape(
            (nTransectTriangles, 3))
        triangulation_args = dict(x=x, y=y, triangles=tris)

        return triangulation_args

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
