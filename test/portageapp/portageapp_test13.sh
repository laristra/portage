#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# 2d 2nd order node-centered remap of quadratic func
mpirun -np 1 $APPDIR/portageapp 13 3 4 y

# Compare the values for the field
$APPDIR/apptest_cmp field_gold13.txt field13.txt 1e-12
