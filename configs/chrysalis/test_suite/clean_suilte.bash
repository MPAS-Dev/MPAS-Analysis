#!/usr/bin/env bash

set -e

branch=$(git symbolic-ref --short HEAD)

rm -rf chrysalis_test_suite
rm -rf /lcrc/group/acme/ac.xylar/analysis_testing/${branch}
rm -rf /lcrc/group/acme/public_html/diagnostic_output/ac.xylar/analysis_testing/${branch}
