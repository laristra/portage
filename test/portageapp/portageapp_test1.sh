#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# 2d 2nd order cell-centered remap of linear func
mpirun -np 4 $APPDIR/portageapp 1 4 8 y

# Compare the values for the field
$APPDIR/apptest_cmp field_gold1.txt field1.txt 1e-12
