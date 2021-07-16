#!/usr/bin/env bash

set -e

machine=$1

branch=$(git symbolic-ref --short HEAD)

rm -rf ${machine}_test_suite
rm -rf /lcrc/group/e3sm/ac.xylar/analysis_testing/${machine}/${branch}
rm -rf /lcrc/group/e3sm/public_html/diagnostic_output/ac.xylar/analysis_testing/${machine}/${branch}
