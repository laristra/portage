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

# 3D, 2nd order accurate remap of quadratic cell-centered function

# SERIAL RUN

mpirun -np 1 $APPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--conformal=y \
--entity_kind=cell --field="x*x+y*y + z*z" \
--remap_order=2 \
--results_file="jali_rect_field_3d_cell_f2_r2.txt"

# Compare the values for the field
$APPDIR/apptest_cmp GOLD_jali_rect_3d_cell_f2_r2.txt jali_rect_field_3d_cell_f2_r2.txt 1e-12

# PARALLEL RUN (note reverse ranks is false here)

mpirun -np 4 $APPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--reverse_ranks=n --conformal=n \
--entity_kind=cell --field="x*x+y*y+z*z" \
--remap_order=2 \
--results_file="jali_rect_field_3d_cell_f2_r2_nc.txt"

# Compare the values for the field
$APPDIR/apptest_cmp GOLD_jali_rect_3d_cell_f2_r2_nc.txt.0 jali_rect_field_3d_cell_f2_r2_nc.txt.0 1e-12
$APPDIR/apptest_cmp GOLD_jali_rect_3d_cell_f2_r2_nc.txt.1 jali_rect_field_3d_cell_f2_r2_nc.txt.1 1e-12
$APPDIR/apptest_cmp GOLD_jali_rect_3d_cell_f2_r2_nc.txt.2 jali_rect_field_3d_cell_f2_r2_nc.txt.2 1e-12
$APPDIR/apptest_cmp GOLD_jali_rect_3d_cell_f2_r2_nc.txt.3 jali_rect_field_3d_cell_f2_r2_nc.txt.3 1e-12
