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

# 3d 2nd order cell-centered remap of quad func on non-conformal mesh
${RUN_COMMAND} $TESTAPPDIR/simple_mesh_app 4 4 5

# Compare the values for the field
$CMPAPPDIR/apptest_cmp field_gold4.txt field4.txt 1e-12
