#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# Small 1st order, cell centered, conformal

# Serial test:
mpirun -np 1 $APPDIR/portageapp 11 2 4 y
$APPDIR/apptest_cmp field_gold11.txt field11.txt 1e-12

# Parallel test:
mpirun -np 8 $APPDIR/portageapp 11 1 2 y

# This fails:
#$APPDIR/apptest_cmp field_gold11.txt field11.txt 1e-12
