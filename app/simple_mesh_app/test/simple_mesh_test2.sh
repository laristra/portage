#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# 3d 2nd order cell-centered remap of linear func
mpirun -np 1 $APPDIR/simple_mesh_app 2 4 5

# Compare the values for the field
$APPDIR/apptest_cmp field_gold2.txt field2.txt 1e-12
