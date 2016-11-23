#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# test a simple 2d remap, save a field
mpirun -np 1 $APPDIR/portageapp 13 3 4 y

# Compare the values for the field
$APPDIR/apptest_cmp field_gold13.txt field13.txt 1e-12
