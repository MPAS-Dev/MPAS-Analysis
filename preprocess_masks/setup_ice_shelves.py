#!/usr/bin/env python
"""
This script creates region groups for ice shelves
"""
from geometric_features import GeometricFeatures, FeatureCollection

import matplotlib.pyplot as plt

plot = True

iceShelfNames = ['Abbot',
                 'Amery',
                 'Atka',
                 'Aviator',
                 'Bach',
                 'Baudouin',
                 'Borchgrevink',
                 'Brahms',
                 'Brunt_Stancomb',
                 'Campbell',
                 'Cheetham',
                 'Conger_Glenzer',
                 'Cook',
                 'Cosgrove',
                 'Crosson',
                 'Dennistoun',
                 'Dibble',
                 'Dotson',
                 'Drygalski',
                 'Edward_VIII',
                 'Ekstrom',
                 'Ferrigno',
                 'Filchner',
                 'Fimbul',
                 'Fitzgerald',
                 'Frost',
                 'GeikieInlet',
                 'George_VI',
                 'Getz',
                 'Gillet',
                 'Hamilton',
                 'Hannan',
                 'HarbordGlacier',
                 'Helen',
                 'Holmes',
                 'HolmesWest',
                 'Hull',
                 'Jelbart',
                 'Land',
                 'Larsen_B',
                 'Larsen_C',
                 'Larsen_D',
                 'Larsen_E',
                 'Larsen_F',
                 'Larsen_G',
                 'Lazarev',
                 'Lillie',
                 'Mariner',
                 'Matusevitch',
                 'Mendelssohn',
                 'Mertz',
                 'Moscow_University',
                 'Moubray',
                 'Mulebreen',
                 'Myers',
                 'Nansen',
                 'Nickerson',
                 'Ninnis',
                 'Nivl',
                 'Noll',
                 'Nordenskjold',
                 'Pine_Island',
                 'PourquoiPas',
                 'Prince_Harald',
                 'Publications',
                 'Quar',
                 'Rayner_Thyer',
                 'Rennick',
                 'Richter',
                 'Riiser-Larsen',
                 'Ronne',
                 'Western_Ross',
                 'Eastern_Ross',
                 'Shackleton',
                 'Shirase',
                 'Slava',
                 'SmithInlet',
                 'Stange',
                 'Sulzberger',
                 'Suvorov',
                 'Swinburne',
                 'Thwaites',
                 'Tinker',
                 'Totten',
                 'Tracy_Tremenchus',
                 'Tucker',
                 'Underwood',
                 'Utsikkar',
                 'Venable',
                 'Verdi',
                 'Vigrid',
                 'Vincennes',
                 'Voyeykov',
                 'West',
                 'Wilkins',
                 'Wilma_Robert_Downer',
                 'Withrow',
                 'Wordie',
                 'Wylde',
                 'Zubchatyy']

combinedIceShelves = {'Filchner-Ronne': ['Filchner', 'Ronne'],
                      'Ross': ['Western_Ross', 'Eastern_Ross'],
                      'Antarctica': ['AntarcticPenninsulaIMBIE',
                                     'WestAntarcticaIMBIE',
                                     'EastAntarcticaIMBIE'],
                      'Peninsula': ['AntarcticPenninsulaIMBIE'],
                      'West Antarctica': ['WestAntarcticaIMBIE'],
                      'East Antarctica': ['EastAntarcticaIMBIE']}

nIMBIEBasins = 27
for basinNumber in range(1, nIMBIEBasins+1):
    basinName = 'Antarctica_IMBIE{}'.format(basinNumber)
    combinedIceShelves['IMBIE{}'.format(basinNumber)] = [basinName]

# create a GeometricFeatures object that points to a local cache of geometric
# data and knows which branch of geometric_feature to use to download
# missing data
gf = GeometricFeatures()

# create a FeatureCollection containing all ice shelves and combined ice-shelf
# regions
fc = FeatureCollection()

# build analysis regions from combining ice shelves from regions with the
# appropriate tags
for shelfName in combinedIceShelves:
    subNames = combinedIceShelves[shelfName]
    print(shelfName)

    print(' * merging features')
    fcShelf = gf.read(componentName='iceshelves', objectType='region',
                      tags=subNames, allTags=False)

    print(' * combining features')
    fcShelf = fcShelf.combine(featureName=shelfName)

    # merge the feature for the basin into the collection of all basins
    fc.merge(fcShelf)

# build ice shelves from regions with the appropriate tags
for shelfName in iceShelfNames:
    print(shelfName)

    print(' * merging features')
    fcShelf = gf.read(componentName='iceshelves', objectType='region',
                      tags=[shelfName])

    print(' * combining features')
    fcShelf = fcShelf.combine(featureName=shelfName)

    # merge the feature for the basin into the collection of all basins
    fc.merge(fcShelf)

# save the feature collection to a geojson file
fc.to_geojson('iceShelves20200621.geojson')

if plot:
    fc.plot(projection='southpole')
    plt.show()
