#!/usr/bin/env python
# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

from setuptools import setup, find_packages
import warnings

isrelease = True

version = '0.7.5'

if not isrelease:
    import subprocess
    try:
        pipe = subprocess.Popen(
            ["git", "describe", "--always", "--match", "v[0-9]*"],
            stdout=subprocess.PIPE)
        (version, stderr) = pipe.communicate()
    except:
        warnings.warn("WARNING: Couldn't get git revision, using generic "
                      "version string")

setup(name='mpas_analysis',
      version=version,
      description='Analysis for Model for Prediction Across Scales (MPAS) '
                  'simulations.',
      url='https://github.com/MPAS-Dev/MPAS-Analysis',
      author='MPAS-Analysis Developers',
      author_email='mpas-developers@googlegroups.com',
      license='BSD',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering',
      ],
      packages=find_packages(),
      package_data={'mpas_analysis': ['config.default'],
                    'mpas_analysis.shared.html': ['templates/*'],
                    'mpas_analysis.test': ['test*/*', 'test*/*/*'],
                    'mpas_analysis.obs': ['analysis_input_files', 'observational_datasets.xml']},
      install_requires=['numpy', 'scipy', 'matplotlib', 'netCDF4', 'xarray',
                        'dask', 'bottleneck', 'basemap', 'lxml', 'nco',
                        'pyproj', 'pillow', 'cmocean', 'progressbar2',
                        'requests'],
      scripts=['run_mpas_analysis', 'download_analysis_data.py'])
