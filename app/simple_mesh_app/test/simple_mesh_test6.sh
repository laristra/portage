#!/bin/bash
: <<'END'
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
END


# Exit on error
set -e
# Echo each command
set -x

# 3d 1st order node-centered remap of quad func
mpirun -np 1 $APPDIR/simple_mesh_app 6 4 5

# Compare the values for the field
$APPDIR/apptest_cmp field_gold6.txt field6.txt 1e-12
