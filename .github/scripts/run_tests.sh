#!/bin/bash

echo -e "## Running tests\n"

echo -e "## Invoking dist_source_me\n"

echo '```'
source ./util/dist_source_me
echo '```'

set -e

cd regression_tests

echo -e "\n## Installing regression test requirements\n"

echo '```'
python -m pip install -r requirements.txt
echo '```'

echo -e "\n## Starting regression tests\n"

echo '```'
python -m pytest --import-mode=importlib -v
echo '```'

echo -e "\n## Tests finished!\n"
