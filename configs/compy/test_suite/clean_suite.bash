#!/usr/bin/env bash

set -e

machine=cori

branch=$(git symbolic-ref --short HEAD)

rm -rf ${machine}_test_suite
rm -rf /compyfs/asay932/analysis_testing/${machine}/${branch}
rm -rf /global/cfs/cdirs/e3sm/www/xylar/analysis_testing/${machine}/${branch}
