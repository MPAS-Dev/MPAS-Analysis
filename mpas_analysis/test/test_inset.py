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

import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from geometric_features import GeometricFeatures
from mpas_analysis.test import TestCase
from mpas_analysis.shared.plot.inset import add_inset


class TestPaths(TestCase):
    def test_add_inset_region(self):
        # Set up the figure and axes
        fig, ax = plt.subplots(figsize=(10, 5),
                               subplot_kw=dict(projection=ccrs.PlateCarree()))

        # Add coastlines to the map
        ax.coastlines()

        gf = GeometricFeatures()
        fc = gf.read(componentName='ocean', objectType='region',
                     featureNames=['North Atlantic Ocean'])

        add_inset(fig, fc)

        # Show the plot
        plt.savefig('inset_region.png')

    def test_add_inset_transect(self):
        # Set up the figure and axes
        fig, ax = plt.subplots(figsize=(10, 5),
                               subplot_kw=dict(projection=ccrs.PlateCarree()))

        # Add coastlines to the map
        ax.coastlines()

        gf = GeometricFeatures()
        fc = gf.read(componentName='ocean', objectType='transect',
                     featureNames=['Drake Passage'])

        add_inset(fig, fc)

        # Show the plot
        plt.savefig('inset_transect.png')

    def test_add_inset_point(self):
        # Set up the figure and axes
        fig, ax = plt.subplots(figsize=(10, 5),
                               subplot_kw=dict(projection=ccrs.PlateCarree()))

        # Add coastlines to the map
        ax.coastlines()

        gf = GeometricFeatures()
        fc = gf.read(componentName='ocean', objectType='point',
                     featureNames=['Equatorial_Pacific_W155.0_N10.0'])

        add_inset(fig, fc)

        # Show the plot
        plt.savefig('inset_point.png')
