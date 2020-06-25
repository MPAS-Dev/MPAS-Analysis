#!/usr/bin/env python
"""
Create a geojson file with the standard transport transects
"""
from geometric_features import GeometricFeatures

# create a GeometricFeatures object that knows which geometric data we want
gf = GeometricFeatures()

# create a FeatureCollection to which we will add all regions
fc = gf.read(componentName='ocean', objectType='transect',
             tags=['standard_transport_sections'])

# save the feature collection to a geojson file
fc.to_geojson('transportTransects20200621.geojson')
