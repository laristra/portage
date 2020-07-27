#!/bin/bash
: <<'END'
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
END

# Echo each command
set -x

FILENAME=simple_mm_test.gmv

rm -f ${FILENAME}*

${RUN_COMMAND} $TESTAPPDIR/simple_mesh_mm_app \
  --dim=2 \
  --nsourcecells=3 \
  --ntargetcells=2 \
  --material_fields=0,1 \
  --material_file=non_convex.bvf \
  --remap_order=1 \
  --results_file=${FILENAME}

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_simple_mm_test.gmv ${FILENAME} 1e-8

