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

DATA_DIR=.

mpirun -np 1 ${TESTAPPDIR}/portageapp_rgmd_jali \
    --problem=ball \
    --dim=2 \
    --nsourcecells=4 \
    --ntargetcells=100 \
    --material_fields="x,x" \
    --output_meshes=y \
    --remap_order=2 \
    --limiter=barth_jespersen \
    --ball_center="1.15,0.625" \
    --ball_radii="0.64"\
    --field_filename="ball_limiter"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_portageapp_rgmd_jali_ball_limiter.gmv ball_limiter0_iteration_0.gmv 1e-9

