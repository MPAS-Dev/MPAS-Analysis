#!/usr/bin/env bash

set -e

machine=cori

branch=$(git symbolic-ref --short HEAD)

rm -rf ${machine}_test_suite
rm -rf /global/cscratch1/sd/xylar/analysis_testing/${machine}/${branch}
rm -rf /global/cfs/cdirs/e3sm/www/xylar/analysis_testing/${machine}/${branch}
