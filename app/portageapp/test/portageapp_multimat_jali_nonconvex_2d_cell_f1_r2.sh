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

${TESTAPPDIR}/portageapp_multimat_jali \
    --dim=2 \
    --source_file=non_convex.exo \
    --ntargetcells=10 \
    --material_fields=2*x-3*y,2*x-3*y \
    --material_file=non_convex.bvf \
    --remap_order=2 \
    --all_convex=n \
    --results_file=multimat_jali_nonconvex_2d_cell_f1_r2.gmv

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_multimat_jali_nonconvex_2d_cell_f1_r2.gmv multimat_jali_nonconvex_2d_cell_f1_r2.gmv  1e-8

