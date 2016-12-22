#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# 3d 2nd order cell-centered remap of quad func on non-conformal mesh
mpirun -np 1 $APPDIR/simple_mesh_app 4 4 5

# Compare the values for the field
$APPDIR/apptest_cmp field_gold4.txt field4.txt 1e-12
