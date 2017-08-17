#!/bin/bash


# Exit on error
set -e
# Echo each command
set -x

# 3d 1st order cell-centered remap of quad func
mpirun -np 1 $APPDIR/simple_mesh_app 1 4 5

# Compare the values for the field
$APPDIR/apptest_cmp field_gold1.txt field1.txt 1e-12
