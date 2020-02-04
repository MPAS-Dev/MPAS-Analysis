#!/usr/bin/env bash
#
# Copyright (c) 2017, UT-BATTELLE, LLC
# All rights reserved.
#
# This software is released under the BSD license detailed
# in the file licenses/A-PRIME_LICENSE
#
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

SOURCE_DIR="../.."

CURRENT="Copyright (c)"

ALWAYS_IGNORE=(-not -path "*.git*" \
               -not -path "${SOURCE_DIR}/docs/*"  \
               -not -path "${SOURCE_DIR}/licenses/*" \
               -not -path "${SOURCE_DIR}/build/*" \
               -not -path "${SOURCE_DIR}/dist/*" \
               -not -path "${SOURCE_DIR}/.idea/*" \
               -not -path "*.egg-info/*" \
               -not -path "${SOURCE_DIR}/conda/recipe/meta.yaml" \
               -not -path "${SOURCE_DIR}/setup.cfg")

FILE_IGNORE=(-not -iname "*.ocean" -not -iname "*_in" \
             -not -iname "*.md" -not -iname "*.json" -not -iname "*.txt" \
             -not -iname "*.png" -not -iname "*.jpg" -not -iname "*.svg" \
             -not -iname "config.*" -not -iname "*.nc" \
             -not -iname "*.pyc" \
             -not -iname "*.sl" -not -iname "*.ps1" -not -iname "*.yml"\
             -not -iname "*.pdf" -not -iname "*.xml" -not -iname "*.m" \
             -not -iname "*.ipynb")

PYTHON_IGNORE=(-not -iname "__init__.py")
