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
#
# Get the source dir and ignore variables
source license-setup.sh

echo "--------------------------------------------------------------------------------"
echo "    ATTEMPTING TO PREPEND A LICENSE HEADER ONTO THESE FILES:"
echo "--------------------------------------------------------------------------------"
find ${SOURCE_DIR} -type f "${ALWAYS_IGNORE[@]}" \
    "${FILE_IGNORE[@]}" \
    "${PYTHON_IGNORE[@]}" \
    | xargs -r grep -L "$CURRENT" \
    | sort

echo "--------------------------------------------------------------------------------"
echo "    BEGIN PREPENDING:"
echo "--------------------------------------------------------------------------------"


############################################################
# Prepend license header to python files with a shebang and an encoding header.
# Will ignore files with a current license header.
############################################################
GET=( -iname "*.py" )

find ${SOURCE_DIR} -type f \( "${GET[@]}" \) "${ALWAYS_IGNORE[@]}" "${PYTHON_IGNORE[@]}" \
    | xargs -r grep -l --max-count=1 "#!" \
    | xargs -r grep -l --max-count=1 "# -\*-" \
    | xargs -r grep -L "$CURRENT" \
    | while read SRC
do
    BN=`basename ${SRC}`
    echo HEADING ${SRC}
    cat ${SRC} | head -n 2 > /tmp/licHead
    cat header-py >> /tmp/licHead
    cat ${SRC} | tail -n +3 >> /tmp/licHead
    chmod --reference=${SRC} /tmp/licHead
    mv /tmp/licHead ${SRC}
done


#######################################################################
# Prepend license header to python files with encoding header.
# Will ignore files with a current license header.
#######################################################################
GET=( -iname "*.py" )

find ${SOURCE_DIR} -type f \( "${GET[@]}" \) "${ALWAYS_IGNORE[@]}" "${PYTHON_IGNORE[@]}" \
    | xargs -r grep -l --max-count=1 "# -\*-" \
    | xargs -r grep -L "$CURRENT" \
    | while read SRC
do
    BN=`basename ${SRC}`
    echo HEADING ${SRC}
    cat ${SRC} | head -n 1 > /tmp/licHead
    cat header-py >> /tmp/licHead
    cat ${SRC} | tail -n +2 >> /tmp/licHead
    chmod --reference=${SRC} /tmp/licHead
    mv /tmp/licHead ${SRC}
done


############################################################
# Prepend license header to python files without a shebang.
# Will ignore files with a current license header.
############################################################
GET=( -iname "*.py" )

find ${SOURCE_DIR} -type f \( "${GET[@]}" \) "${ALWAYS_IGNORE[@]}" "${PYTHON_IGNORE[@]}" \
    | xargs -r grep -L "#!" \
    | xargs -r grep -L "$CURRENT" \
    | while read SRC
do
    BN=`basename ${SRC}`
    echo HEADING ${SRC}
    cp header-py /tmp/licHead
    cat ${SRC} >> /tmp/licHead
    chmod --reference=${SRC} /tmp/licHead
    mv /tmp/licHead ${SRC}
done


#######################################################################
# Prepend license header to python, bash, and sh files with a shebang.
# Will ignore files with a current license header.
#######################################################################
GET=( -iname "*.py" -or -iname "*.sh" -or -iname "*.bash" -or -iname "*.csh" -or -iname "*.pbs")

find ${SOURCE_DIR} -type f \( "${GET[@]}" \) "${ALWAYS_IGNORE[@]}" "${PYTHON_IGNORE[@]}" \
    | xargs -r grep -l --max-count=1 "#!" \
    | xargs -r grep -L "$CURRENT" \
    | while read SRC
do
    BN=`basename ${SRC}`
    echo HEADING ${SRC}
    cat ${SRC} | head -1 > /tmp/licHead
    cat header-py >> /tmp/licHead
    cat ${SRC} | tail -n +2 >> /tmp/licHead
    chmod --reference=${SRC} /tmp/licHead
    mv /tmp/licHead ${SRC}
done


####################################################
# Prepend license header to html files.
# Will ignore files with a current license header.
####################################################
GET=( -iname "*.html")

find ${SOURCE_DIR} -type f \( "${GET[@]}" \) "${ALWAYS_IGNORE[@]}" \
    | xargs -r grep -L "$CURRENT" \
    | while read SRC
do
    BN=`basename ${SRC}`
    echo HEADING ${SRC}
    cp header-html /tmp/licHead
    cat ${SRC} >> /tmp/licHead
    chmod --reference=${SRC} /tmp/licHead
    mv /tmp/licHead ${SRC}
done

####################################################
# Prepend license header to css and js files.
# Will ignore files with a current license header.
####################################################
GET=( -iname "*.css" -or -iname "*.js" )

find ${SOURCE_DIR} -type f \( "${GET[@]}" \) "${ALWAYS_IGNORE[@]}" \
    | xargs -r grep -L "$CURRENT" \
    | while read SRC
do
    BN=`basename ${SRC}`
    echo HEADING ${SRC}
    cp header-css /tmp/licHead
    cat ${SRC} >> /tmp/licHead
    chmod --reference=${SRC} /tmp/licHead
    mv /tmp/licHead ${SRC}
done

echo "--------------------------------------------------------------------------------"
echo "    DONE PREPENDING!"


MISSED=`find ${SOURCE_DIR} -type f "${ALWAYS_IGNORE[@]}" \
    "${FILE_IGNORE[@]}" \
    "${PYTHON_IGNORE[@]}" \
    | xargs -r grep -L "$CURRENT"`

if [ "$MISSED" ]
then
    echo ""
    echo "    WARNING: There is no method to prepend a license header onto the following"
    echo "             files. Please either manually add the license head to these files,"
    echo "             add a method to prepend a license header onto these type of files,"
    echo "             or set them to be ignored in 'license-setup.sh'. "
    echo "--------------------------------------------------------------------------------"
    find ${SOURCE_DIR} -type f "${ALWAYS_IGNORE[@]}" \
        "${FILE_IGNORE[@]}" \
        "${PYTHON_IGNORE[@]}" \
        | xargs -r grep -L "$CURRENT" \
        | sort
    echo "--------------------------------------------------------------------------------"
else
    echo "--------------------------------------------------------------------------------"
fi
