#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# 3d 1st order node-centered remap of linear func
mpirun -np 1 $APPDIR/simple_mesh_app 5 4 5

# Compare the values for the field
$APPDIR/apptest_cmp field_gold5.txt field5.txt 1e-12
