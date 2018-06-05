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
#
# Get the source dir and ignore variables
source license-setup.sh

OLD="Copyright (c) 2017, UT"
NEW="Copyright (c) 2017-2018, UT"

echo "--------------------------------------------------------------------------------"
echo "    THESE FILES HAVE AN OUTDATED LICENSE HEADER:"
echo "--------------------------------------------------------------------------------"
find ${SOURCE_DIR} -type f "${ALWAYS_IGNORE[@]}" \
    "${FILE_IGNORE[@]}" \
    "${PYTHON_IGNORE[@]}" \
    | xargs -r grep -l "$OLD" \
    | sort

# Also, bump any copyright statements in the documentation
find ${SOURCE_DIR}/docs -type f \
    -not -path "*_build*" \
    -not -path "*source*" \
    "${FILE_IGNORE[@]}" \
    | xargs -r grep -l "$OLD" \
    | sort


echo "--------------------------------------------------------------------------------"
echo "    UPDATING THE LICENSE HEADERS:"
echo "--------------------------------------------------------------------------------"
find ${SOURCE_DIR} -type f "${ALWAYS_IGNORE[@]}" \
    "${FILE_IGNORE[@]}" \
    | xargs -r grep -l "$OLD" \
    | sort \
    | while read SRC
do
    echo "Bumping $SRC"
    sed -i "s/$OLD/$NEW/g" ${SRC}
done

# Also, bump any copyright statements in the documentation
find ${SOURCE_DIR}/docs -type f \
    -not -path "*_build*" \
    -not -path "*source*" \
    "${FILE_IGNORE[@]}" \
    | xargs -r grep -l "$OLD" \
    | sort \
    | while read SRC
do
    echo "Bumping $SRC"
    sed -i "s/$OLD/$NEW/g" ${SRC}
done
