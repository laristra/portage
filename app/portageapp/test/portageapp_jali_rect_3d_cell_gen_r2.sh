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

# 3D, 2nd order accurate remap of general cell-centered function

# SERIAL RUN, conforming meshes

mpirun -np 1 $TESTAPPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--conformal=y \
--entity_kind=cell --field="sin(10*x)+cos(10*y)+tan(z)" \
--remap_order=2 \
--results_file="jali_rect_3d_cell_gen_r2.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_cell_gen_r2.txt jali_rect_3d_cell_gen_r2.txt 1e-12


# PARALLEL RUN, conforming meshes

mpirun -np 4 $TESTAPPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--conformal=y \
--entity_kind=cell --field="sin(10*x)+cos(10*y)+tan(z)" \
--remap_order=2 \
--results_file="jali_rect_3d_cell_gen_r2.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_cell_gen_r2.txt.0 jali_rect_3d_cell_gen_r2.txt.0 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_cell_gen_r2.txt.1 jali_rect_3d_cell_gen_r2.txt.1 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_cell_gen_r2.txt.2 jali_rect_3d_cell_gen_r2.txt.2 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_cell_gen_r2.txt.3 jali_rect_3d_cell_gen_r2.txt.3 1e-12


# PARALLEL RUN, non-conforming meshes

mpirun -np 4 $TESTAPPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--conformal=n \
--entity_kind=cell --field="sin(10*x)+cos(10*y)+tan(z)" \
--remap_order=2 \
--results_file="jali_rect_3d_cell_gen_r2_nc.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_cell_gen_r2_nc.txt.0 jali_rect_3d_cell_gen_r2_nc.txt.0 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_cell_gen_r2_nc.txt.1 jali_rect_3d_cell_gen_r2_nc.txt.1 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_cell_gen_r2_nc.txt.2 jali_rect_3d_cell_gen_r2_nc.txt.2 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_cell_gen_r2_nc.txt.3 jali_rect_3d_cell_gen_r2_nc.txt.3 1e-12
