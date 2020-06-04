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
    --problem=tjunction \
    --dim=2 \
    --nsourcecells=10 \
    --ntargetcells=10 \
    --material_fields="1,x,y" \
    --intersect=n \
    --perturb_source=pseudorandom \
    --source_convex_cells=n \
    --output_meshes=y

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_portageapp_rgmd_jali_tjunction_2d_swept.gmv target_mm_0.gmv 1e-9

