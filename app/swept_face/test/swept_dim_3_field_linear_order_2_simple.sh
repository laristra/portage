#!/bin/bash
: <<'END'
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
END

set -e
set -x

${RUN_COMMAND} $TESTAPPDIR/swept_face_demo \
  --dim=3 \
  --ncells=10 \
  --remap_order=2 \
  --field="x+2y+3z" \
  --ntimesteps=4 \
  --scale_by=10 \
  --output_meshes=false \
  --result_file="swept_dim_3_field_linear_order_2_simple" \
  --keep_source=false \
  --simple=true
