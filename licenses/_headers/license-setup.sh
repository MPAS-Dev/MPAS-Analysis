#!/usr/bin/env bash
#
# Copyright (c) 2017, UT-BATTELLE, LLC
# All rights reserved.
#
# This software is released under the BSD license detailed
# in the file licenses/A-PRIME_LICENSE
#
# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html

SOURCE_DIR="../.."

CURRENT="Copyright (c)"

ALWAYS_IGNORE=(-not -path "*.git*" -not -path "*docs/*" -not -path "build/*" \
               -not -path "dist/*" -not -path "*.egg-info/*" \
               -not -path "licenses/*" )

FILE_IGNORE=(-not -iname "*.md" -not -iname "*.json" -not -iname "*.txt" \
             -not -iname "*.png" -not -iname "*.jpg" -not -iname "*.svg" \
             -not -iname "config.*" -not -iname "*.nc"\
             -not -iname "*.pyc" \
             -not -iname "*.sl" -not -iname "*.ps1" -not -iname "*.yml"    )

PYTHON_IGNORE=(-not -iname "__init__.py")
