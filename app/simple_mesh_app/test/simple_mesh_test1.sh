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

# 3d 1st order cell-centered remap of quad func
${RUN_COMMAND} $TESTAPPDIR/simple_mesh_app 1 4 5

# Compare the values for the field
$CMPAPPDIR/apptest_cmp field_gold1.txt field1.txt 1e-12
