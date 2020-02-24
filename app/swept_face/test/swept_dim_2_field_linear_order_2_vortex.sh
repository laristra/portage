#!/bin/bash
: <<'END'
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
END

set -e
set -x

epsilon=1.e-10

# run test
${RUN_COMMAND} $TESTAPPDIR/swept_face_demo \
  --dim=2 \
  --ncells=10 \
  --remap_order=2 \
  --field="x+2y" \
  --ntimesteps=4 \
  --scale_by=10 \
  --output_meshes=false \
  --result_file="data_dim_2_field_linear_order_2_vortex" \
  --keep_source=false \
  --simple=false

# compare values with related gold file
$CMPAPPDIR/apptest_cmp \
  "gold_dim_2_field_linear_order_2_vortex.txt" \
  "data_dim_2_field_linear_order_2_vortex.txt" \
  ${epsilon}