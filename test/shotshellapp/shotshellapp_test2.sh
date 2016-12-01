#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

DATA_DIR=$APPDIR/../test/wrappers

# shotshellapp only works in serial for now:
mpirun -np 1 ${APPDIR}/shotshellapp 2 ${DATA_DIR}/shotshell.exo ${DATA_DIR}/shotshell-v.exo y

# Compare the values for the field
$APPDIR/apptest_cmp field_gold2.txt field2.txt 1e-9
