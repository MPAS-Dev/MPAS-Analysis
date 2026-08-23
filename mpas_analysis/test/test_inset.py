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
import matplotlib.pyplot as plt  # noqa: E402
import cartopy.crs as ccrs  # noqa: E402
import pytest  # noqa: E402

from geometric_features import GeometricFeatures  # noqa: E402
from mpas_analysis.shared.plot.inset import add_inset  # noqa: E402


@pytest.mark.parametrize(
    'object_type,feature_name,filename',
    [('region', 'North Atlantic Ocean', 'inset_region.png'),
     ('transect', 'Drake Passage', 'inset_transect.png'),
     ('point', 'Equatorial_Pacific_W155.0_N10.0', 'inset_point.png')])
def test_add_inset(tmp_path, object_type, feature_name, filename):
    # Set up the figure and axes
    fig, ax = plt.subplots(figsize=(10, 5),
                           subplot_kw=dict(projection=ccrs.PlateCarree()))
    try:
        # Add coastlines to the map
        ax.coastlines()

        gf = GeometricFeatures()
        fc = gf.read(componentName='ocean', objectType=object_type,
                     featureNames=[feature_name])

        add_inset(fig, fc)

        out_filename = tmp_path / filename
        fig.savefig(out_filename)
        assert out_filename.exists()
    finally:
        plt.close(fig)
