#!/usr/bin/env python
"""
This script creates a geojson file with Antarctic regions similar to
Timmermann et al. 2013
"""
from geometric_features import GeometricFeatures

regions = ['Southern Ocean', 'Southern Ocean 60S', 'Eastern Weddell Sea Shelf',
           'Eastern Weddell Sea Deep', 'Western Weddell Sea Shelf',
           'Western Weddell Sea Deep', 'Weddell Sea Shelf', 'Weddell Sea Deep',
           'Bellingshausen Sea Shelf', 'Bellingshausen Sea Deep',
           'Amundsen Sea Shelf', 'Amundsen Sea Deep', 'Eastern Ross Sea Shelf',
           'Eastern Ross Sea Deep', 'Western Ross Sea Shelf',
           'Western Ross Sea Deep', 'East Antarctic Seas Shelf',
           'East Antarctic Seas Deep']

# create a GeometricFeatures object that knows which geometric data we want
gf = GeometricFeatures()

# create a FeatureCollection to which we will add all regions
fc = gf.read(componentName='ocean', objectType='region', featureNames=regions)

# save the feature collection to a geojson file
fc.to_geojson('antarcticRegions20200621.geojson')
