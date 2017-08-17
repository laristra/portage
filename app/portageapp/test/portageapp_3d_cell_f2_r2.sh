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
--entity_kind=cell --field_order=2 \
--remap_order=2 \
--output_results=y

# Compare the values for the field
$APPDIR/apptest_cmp GOLD-field_3d_cell_f2_r2.txt field_3d_cell_f2_r2.txt 1e-12

# PARALLEL RUN (note reverse ranks is false here)

mpirun -np 4 $APPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--reverse_ranks=n --conformal=n \
--entity_kind=cell --field_order=2 \
--remap_order=2 \
--output_results=y

# Compare the values for the field
$APPDIR/apptest_cmp GOLD-field_3d_cell_f2_r2_nc.txt.0 field_3d_cell_f2_r2_nc.txt.0 1e-12
$APPDIR/apptest_cmp GOLD-field_3d_cell_f2_r2_nc.txt.1 field_3d_cell_f2_r2_nc.txt.1 1e-12
$APPDIR/apptest_cmp GOLD-field_3d_cell_f2_r2_nc.txt.2 field_3d_cell_f2_r2_nc.txt.2 1e-12
$APPDIR/apptest_cmp GOLD-field_3d_cell_f2_r2_nc.txt.3 field_3d_cell_f2_r2_nc.txt.3 1e-12
