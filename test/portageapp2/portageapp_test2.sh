#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# test a simple 2d remap
mpirun -np 4 $APPDIR/portageapp 0 4 16 y

# in a future version, one or more comparators can go here

