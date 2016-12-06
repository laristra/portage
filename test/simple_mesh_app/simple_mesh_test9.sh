#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# 3d 2nd order node-centered remap of quad func on a non-conformal mesh
mpirun -np 1 $APPDIR/simple_mesh 9 4 5

# Compare the values for the field
$APPDIR/apptest_cmp field_gold9.txt field9.txt 1e-12
