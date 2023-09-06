#!/bin/bash
echo "**** Invoking dist_source_me"
source ./util/dist_source_me

set -e

cd regression_tests
echo "**** Starting Regression Tests"
./scripts/run_tests.py -test all

echo "**** Tests finished"
