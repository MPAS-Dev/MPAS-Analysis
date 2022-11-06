#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

from setuptools import setup, find_packages
import os
import re


install_requires = \
    ['bottleneck',
     'cartopy>=0.18.0',
     'cmocean',
     'dask',
     'gsw',
     'lxml',
     'matplotlib >=3.0.2',
     'netcdf4',
     'numpy',
     'pandas',
     'pillow',
     'progressbar2',
     'pyproj',
     'python-dateutil',
     'requests',
     'scipy',
     'setuptools',
     'shapely',
     'six',
     'xarray>=0.14.1']

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'mpas_analysis', '__init__.py')) as f:
    init_file = f.read()

version = re.search(r'{}\s*=\s*[(]([^)]*)[)]'.format('__version_info__'),
                    init_file).group(1).replace(', ', '.')

setup(name='mpas_analysis',
      version=version,
      description='Analysis for Model for Prediction Across Scales (MPAS) '
                  'simulations.',
      url='https://github.com/MPAS-Dev/MPAS-Analysis',
      author='MPAS-Analysis Developers',
      author_email='mpas-developers@googlegroups.com',
      license='BSD',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3.11',
          'Topic :: Scientific/Engineering',
      ],
      packages=find_packages(),
      package_data={'mpas_analysis': ['*.cfg',
                                      'obs/analysis_input_files',
                                      'obs/sose_10000km_input_files',
                                      'obs/observational_datasets.xml'],
                    'mpas_analysis.configuration': ['*.cfg'],
                    'mpas_analysis.shared.html': ['templates/*'],
                    'mpas_analysis.test': ['test*/*', 'test*/*/*'],
                    'mpas_analysis.shared.plot':
                        ['ScientificColourMaps7/*/*.xml',
                         'SciVisColorColormaps/*.xml']},
      install_requires=install_requires,
      entry_points={'console_scripts':
                    ['mpas_analysis = mpas_analysis.__main__:main',
                     'download_analysis_data = '
                     'mpas_analysis.download_data:download_analysis_data']})
