#!/bin/bash


# Exit on error
set -e
# Echo each command
set -x

# 3d 2nd order node-centered remap of quad func
mpirun -np 1 $APPDIR/simple_mesh_app 8 4 5

# Compare the values for the field
$APPDIR/apptest_cmp field_gold8.txt field8.txt 1e-12
