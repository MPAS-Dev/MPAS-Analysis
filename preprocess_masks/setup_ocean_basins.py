#!/usr/bin/env python
"""
This script creates a geojson file with that includes the Global Ocean as well
as Atlantic, Pacific, Indian, Arctic, Southern Ocean, Equatorial (global
15S-15N), and Mediterranean basins
"""

import copy
import shapely.geometry
import shapely.ops

from geometric_features import GeometricFeatures, FeatureCollection


def build_ocean_basins(gf):
    """
    Builds features defining the major ocean basins

    Parameters
    ----------
    gf : ``GeometricFeatures``
        An object that knows how to download and read geometric featuers

    Returns
    -------
    fc : ``FeatureCollection``
        The new feature collection
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    fc = FeatureCollection()
    fc.set_group_name(groupName='OceanBasinRegionsGroup')

    # build ocean basins from regions with the appropriate tags
    for oceanName in ['Atlantic', 'Pacific', 'Indian', 'Arctic',
                      'Southern_Ocean', 'Mediterranean']:

        basinName = '{}_Basin'.format(oceanName)
        print(oceanName)

        print(' * merging features')
        fcBasin = gf.read(componentName='ocean', objectType='region',
                          tags=[basinName])

        print(' * combining features')
        fcBasin = fcBasin.combine(featureName=basinName)

        fc.merge(fcBasin)

    # add the global ocean, global ocean between 65S and 65S, and
    # equatorial region
    fc.merge(gf.read(componentName='ocean', objectType='region',
                     featureNames=['Global Ocean',
                                   'Global Ocean 65N to 65S',
                                   'Global Ocean 15S to 15N']))

    return fc


def remove_small_polygons(fc, minArea):
    """
    A helper function to remove small polygons from a feature collection

    Parameters
    ----------
    fc : ``FeatureCollection``
        The feature collection to remove polygons from

    minArea : float
        The minimum area (in square degrees) below which polygons should be
        removed

    Returns
    -------
    fcOut : ``FeatureCollection``
        The new feature collection with small polygons removed
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    fcOut = FeatureCollection()

    removedCount = 0
    for feature in fc.features:
        geom = feature['geometry']
        if geom['type'] not in ['Polygon', 'MultiPolygon']:
            # no area to check, so just add it
            add = True
        else:
            add = False
            featureShape = shapely.geometry.shape(geom)
            if featureShape.type == 'Polygon':
                if featureShape.area > minArea:
                    add = True
                else:
                    removedCount += 1
            else:
                # a MultiPolygon
                outPolygons = []
                for polygon in featureShape:
                    if polygon.area > minArea:
                        outPolygons.append(polygon)
                    else:
                        removedCount += 1
                if len(outPolygons) > 0:
                    outShape = shapely.ops.cascaded_union(outPolygons)
                    feature['geometry'] = shapely.geometry.mapping(outShape)
                    add = True
        if add:
            fcOut.add_feature(copy.deepcopy(feature))
        else:
            print("{} has been removed.".format(
                    feature['properties']['name']))

    print(' * Removed {} small polygons'.format(removedCount))

    return fcOut


def main():
    # create a GeometricFeatures object that knows which geometric data we want
    gf = GeometricFeatures()

    fcOceanBasins = build_ocean_basins(gf)
    fcOceanBasins.to_geojson('oceanBasins20200621.geojson')


if __name__ == '__main__':
    main()