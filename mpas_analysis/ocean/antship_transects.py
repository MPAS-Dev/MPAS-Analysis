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
from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.compute_transects_subtask import \
    ComputeTransectsSubtask, TransectsObservations

from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask

from mpas_analysis.shared.io.utility import build_obs_path

from collections import OrderedDict


class AntshipTransects(AnalysisTask):
    """
    Plot model output and compare it against Antarctic ship transects (https://doi.org/10.17882/99787)
    """
    # Authors
    # -------
    # Irena Vankova

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
        # Irena Vankova

        tags = ['transect', 'antarctic', 'publicObs']

        # call the constructor from the base class (AnalysisTask)
        super(AntshipTransects, self).__init__(
            config=config, taskName='antshipTransects',
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
            verticalComparisonGrid = config.getnumpy(
                sectionName, 'verticalComparisonGrid')

        verticalBounds = config.getexpression(sectionName, 'verticalBounds')

        horizontalBounds = config.getexpression(
            sectionName, 'horizontalBounds')

        #observationsDirectory = build_obs_path(
        #    config, 'ocean', 'woceSubdirectory')

        observationsDirectory = build_obs_path(
            config, 'ocean', 'antshipSubdirectory')

        origObsFileNames = \
            {'Belling-BelgicaTrough1-2007': 'Belling-BelgicaTrough1-2007.nc',
             'Belling-BelgicaTrough2-2007': 'Belling-BelgicaTrough2-2007.nc',
             'Belling-Cross1-2007': 'Belling-Cross1-2007.nc',
             'Belling-Cross2-2007': 'Belling-Cross2-2007.nc',
             'Belling-Cross3-2007': 'Belling-Cross3-2007.nc',
             'Belling-EltaninBay-Along-2007': 'Belling-EltaninBay-Along-2007.nc',
             'Belling-GeorgeVI-AlongW1-2007': 'Belling-GeorgeVI-AlongW1-2007.nc',
             'Belling-GeorgeVI-AlongW2-2007': 'Belling-GeorgeVI-AlongW2-2007.nc',
             'Amundsen-AcrossPIG1-2020': 'Amundsen-AcrossPIG1-2020.nc',
             'Amundsen-AcrossPIG2-2020': 'Amundsen-AcrossPIG2-2020.nc',
             'Amundsen-AcrossPIG3-2020': 'Amundsen-AcrossPIG3-2020.nc',
             'Amundsen-E-Dodson1-1994': 'Amundsen-E-Dodson1-1994.nc',
             'Amundsen-E-Dodson2-1994': 'Amundsen-E-Dodson2-1994.nc',
             'Amundsen-E-PIG-2020': 'Amundsen-E-PIG-2020.nc',
             'Amundsen-E-ThwaitesCrosson-1994': 'Amundsen-E-ThwaitesCrosson-1994.nc',
             'Amundsen-W-PIG-1994': 'Amundsen-W-PIG-1994.nc',
             'Amundsen-W-PIG-2009': 'Amundsen-W-PIG-2009.nc',
             'Amundsen-W-PIG-2014': 'Amundsen-W-PIG-2014.nc',
             'Amundsen-W-PIG-2019': 'Amundsen-W-PIG-2019.nc',
             'Ross-C1-2004': 'Ross-C1-2004.nc',
             'Ross-C2-2004': 'Ross-C2-2004.nc',
             'Ross-C3-2004': 'Ross-C3-2004.nc',
             'Ross-C4-2004': 'Ross-C4-2004.nc',
             'Ross-C5-2004': 'Ross-C5-2004.nc',
             'Ross-C6-2004': 'Ross-C6-2004.nc',
             'Ross-C7-2004': 'Ross-C7-2004.nc',
             'Ross-Shelf1-1994': 'Ross-Shelf1-1994.nc',
             'Amery-AlongFront-2002': 'Amery-AlongFront-2002.nc',
             'Amery-AlongFront-2006': 'Amery-AlongFront-2006.nc',
             'Amery-Lon73.5-1992': 'Amery-Lon73.5-1992.nc',
             'Amery-Lon73.5-2000': 'Amery-Lon73.5-2000.nc',
             'Amery-Lon75.5-1992': 'Amery-Lon75.5-1992.nc',
             'Amery-Lon75.5-2008': 'Amery-Lon75.5-2008.nc',
             'DML-RoiB-E-2008': 'DML-RoiB-E-2008.nc',
             'DML-RoiB-E-2019': 'DML-RoiB-E-2019.nc',
             'DML-RoiB-C-1996': 'DML-RoiB-C-1996.nc',
             'DML-RoiB-C-2006': 'DML-RoiB-C-2006.nc',
             'DML-RoiB-W-2003': 'DML-RoiB-W-2003.nc',
             'DML-Munin-Lon18-2003': 'DML-Munin-Lon18-2003.nc',
             'DML-Fimbul-E-Lon6.25-1991': 'DML-Fimbul-E-Lon6.25-1991.nc',
             'DML-Fimbul-C-Lon0-1992': 'DML-Fimbul-C-Lon0-1992.nc',
             'DML-Fimbul-C-Lon0-1996': 'DML-Fimbul-C-Lon0-1996.nc',
             'DML-Fimbul-C-Lon0-2008': 'DML-Fimbul-C-Lon0-2008.nc',
             'Weddell-RiiserLarsen-C1-1995': 'Weddell-RiiserLarsen-C1-1995.nc',
             'Weddell-RiiserLarsen-C2-1985': 'Weddell-RiiserLarsen-C2-1985.nc',
             'Weddell-RiiserLarsen-C2-2009': 'Weddell-RiiserLarsen-C2-2009.nc',
             'Weddell-RiiserLarsen-C2-2010': 'Weddell-RiiserLarsen-C2-2010.nc',
             'Weddell-RiiserLarsen-C3-1985': 'Weddell-RiiserLarsen-C3-1985.nc',
             'Weddell-RiiserLarsen-C4-1985': 'Weddell-RiiserLarsen-C4-1985.nc',
             'Weddell-RiiserLarsen-C5-1985': 'Weddell-RiiserLarsen-C5-1985.nc',
             'Weddell-RiiserLarsen-C6-1987': 'Weddell-RiiserLarsen-C6-1987.nc',
             'Weddell-StancombBrunt-1990': 'Weddell-StancombBrunt-1990.nc',
             'Weddell-Filchner-C1-2009': 'Weddell-Filchner-C1-2009.nc',
             'Weddell-Filchner-C2-2017': 'Weddell-Filchner-C2-2017.nc',
             'Weddell-Filchner-C3-1995': 'Weddell-Filchner-C3-1995.nc',
             'Weddell-Filchner-C4-1977': 'Weddell-Filchner-C4-1977.nc',
             'Weddell-Filchner-C4-1990': 'Weddell-Filchner-C4-1990.nc',
             'Weddell-Filchner-C5-1978': 'Weddell-Filchner-C5-1978.nc',
             'Weddell-FilchnerTrough-1995': 'Weddell-FilchnerTrough-1995.nc',
             'Weddell-FilchnerTrough-1999': 'Weddell-FilchnerTrough-1999.nc',
             'Weddell-Ronne-C1-1999': 'Weddell-Ronne-C1-1999.nc',
             'Weddell-FRISAlongFront-1995': 'Weddell-FRISAlongFront-1995.nc',
             'Weddell-FRISAlongFront-2018': 'Weddell-FRISAlongFront-2018.nc',
             'Weddell-Larsen-C1-1992': 'Weddell-Larsen-C1-1992.nc',
             'Weddell-Larsen-C2-1992': 'Weddell-Larsen-C2-1992.nc',
             'Weddell-Larsen-AlongBreak-1992': 'Weddell-Larsen-AlongBreak-1992.nc'}

        obsFileNames = {}
        for transectName in horizontalBounds:
            found = False
            for name in origObsFileNames:
                if transectName.startswith(name):
                    obsFileNames[transectName] = origObsFileNames[name]
                    found = True
                    break
            if not found:
                raise ValueError(f'Keys for horizontalBounds must start '
                                 f'with one of {list(origObsFileNames)}')

        for transectName in obsFileNames:
            fileName = '{}/{}'.format(observationsDirectory,
                                      obsFileNames[transectName])
            obsFileNames[transectName] = fileName

        fields = \
            {'temperature':
                {'mpas': 'timeMonthly_avg_activeTracers_temperature',
                 'obs': 'potentialTemperature',
                 'titleName': 'Potential Temperature',
                 'units': r'$\degree$C'},
             'salinity':
                {'mpas': 'timeMonthly_avg_activeTracers_salinity',
                 'obs': 'salinity',
                 'titleName': 'Salinity',
                 'units': r'PSU'},
             'potentialDensity':
                {'mpas': 'timeMonthly_avg_potentialDensity',
                 'obs': 'potentialDensity',
                 'titleName': 'Potential Density',
                 'units': r'kg m$^{-3}$'},
             'potentialDensityContour':
                 {'mpas': 'timeMonthly_avg_potentialDensity',
                  'obs': 'potentialDensity',
                  'titleName': 'Potential Density Contours',
                  'units': r'kg m$^{-3}$'}}

        transectCollectionName = 'AntShip_transects'
        if horizontalResolution not in ['obs', 'mpas']:
            transectCollectionName = \
                f'{transectCollectionName}_{horizontalResolution}km'

        transectsObservations = TransectsObservations(config, obsFileNames,
                                                      horizontalResolution,
                                                      transectCollectionName)

        computeTransectsSubtask = ComputeTransectsSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='AntShip',
            transectCollectionName=transectCollectionName,
            variableList=[field['mpas'] for field in fields.values()],
            seasons=seasons,
            obsDatasets=transectsObservations,
            verticalComparisonGridName=verticalComparisonGridName,
            verticalComparisonGrid=verticalComparisonGrid)

        plotObs = controlConfig is None
        if plotObs:

            refTitleLabel = 'Observations (ship CTD)'

            diffTitleLabel = 'Model - Observations'

        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = f'Control: {controlRunName}'

            diffTitleLabel = 'Main - Control'

        fieldNameDict = {'temperature': 'temperatureTransect',
                         'salinity': 'salinityTransect',
                         'potentialDensity': 'potentialDensityTransect',
                         'potentialDensityContour':
                             'potentialDensityContourTransect'}

        for fieldName in fields:
            for transectName in obsFileNames:
                for season in seasons:
                    outFileLabel = fieldNameDict[fieldName]
                    if plotObs:
                        refFieldName = fields[fieldName]['obs']
                    else:
                        refFieldName = fields[fieldName]['mpas']

                    fieldNameUpper = fieldName[0].upper() + fieldName[1:]
                    titleName = fields[fieldName]['titleName']
                    fieldNameInTitle = \
                        f'{titleName} from {transectName.replace("_", " ")}'

                    # make a new subtask for this season and comparison grid
                    subtask = PlotTransectSubtask(
                        self, season, transectName, fieldName,
                        computeTransectsSubtask, plotObs, controlConfig,
                        horizontalBounds[transectName])

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=fieldNameInTitle,
                        mpasFieldName=fields[fieldName]['mpas'],
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=fields[fieldName]['units'],
                        imageCaption=f'{fieldNameInTitle} {season}',
                        galleryGroup='Antarctic Ship Transects',
                        groupSubtitle=None,
                        groupLink='antship',
                        galleryName=titleName,
                        configSectionName=f'antship{fieldNameUpper}Transects',
                        verticalBounds=verticalBounds)

                    self.add_subtask(subtask)
