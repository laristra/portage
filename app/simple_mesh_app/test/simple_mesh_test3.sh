#!/bin/bash


# Exit on error
set -e
# Echo each command
set -x

# 3d 2nd order cell-centered remap of quad func
mpirun -np 1 $APPDIR/simple_mesh_app 3 4 5

# Compare the values for the field
$APPDIR/apptest_cmp field_gold3.txt field3.txt 1e-12
