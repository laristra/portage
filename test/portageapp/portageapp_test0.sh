#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# 2d 1st order cell-centered remap of linear func
mpirun -np 4 $APPDIR/portageapp 0 4 8 y

# Compare the values for the field
$APPDIR/apptest_cmp field_gold0.txt field0.txt 1e-12
