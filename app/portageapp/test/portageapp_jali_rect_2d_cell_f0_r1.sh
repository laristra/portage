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

# 2D, 1st order accurate remap of constant cell-centered function

# SERIAL RUN (conformal)

mpirun -np 1 $TESTAPPDIR/portageapp_jali \
--dim=2 --nsourcecells=5 --ntargetcells=7 \
--conformal=y \
--entity_kind=cell --field="73.98" \
--remap_order=1 \
--results_file="jali_rect_2d_cell_f0_r1.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_2d_cell_f0_r1.txt jali_rect_2d_cell_f0_r1.txt 1e-12

# SERIAL RUN (non-conformal)

mpirun -np 1 $TESTAPPDIR/portageapp_jali \
--dim=2 --nsourcecells=5 --ntargetcells=7 \
--conformal=n \
--entity_kind=cell --field="73.98" \
--remap_order=1 \
--results_file="jali_rect_2d_cell_f0_r1_nc.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_2d_cell_f0_r1_nc.txt jali_rect_2d_cell_f0_r1_nc.txt 1e-12

# PARALLEL RUN (non-conformal)

mpirun -np 4 $TESTAPPDIR/portageapp_jali \
--dim=2 --nsourcecells=5 --ntargetcells=7 \
--conformal=n \
--entity_kind=cell --field="73.98" \
--remap_order=1 \
--results_file="jali_rect_2d_cell_f0_r1_nc.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_2d_cell_f0_r1_nc.txt.0 jali_rect_2d_cell_f0_r1_nc.txt.0 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_2d_cell_f0_r1_nc.txt.1 jali_rect_2d_cell_f0_r1_nc.txt.1 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_2d_cell_f0_r1_nc.txt.2 jali_rect_2d_cell_f0_r1_nc.txt.2 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_2d_cell_f0_r1_nc.txt.3 jali_rect_2d_cell_f0_r1_nc.txt.3 1e-12
