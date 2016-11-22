#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# test a simple 2d remap, save a field
mpirun -np 4 $APPDIR/portageapp 1 4 8 y

# Compare the values for the field
$APPDIR/apptest_cmp field_gold1.txt field.txt 1e-12
