#!/bin/bash

# Exit on error
set -e
# Echo each command
set -x

# 3D, 1st order accurate remap of linear cell-centered function

# SERIAL RUN

$APPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--conformal=y \
--entity_kind=cell --field_order=1 \
--remap_order=1 \
--output_results=y

# Compare the values for the field
$APPDIR/apptest_cmp GOLD-field_3d_cell_f1_r1.txt field_3d_cell_f1_r1.txt 1e-12

# PARALLEL RUN

mpirun -np 4 $APPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--conformal=n \
--entity_kind=cell --field_order=1 \
--remap_order=1 \
--output_results=y

# Compare the values for the field
$APPDIR/apptest_cmp GOLD-field_3d_cell_f1_r1_nc.txt.0 field_3d_cell_f1_r1_nc.txt.0 1e-12
$APPDIR/apptest_cmp GOLD-field_3d_cell_f1_r1_nc.txt.1 field_3d_cell_f1_r1_nc.txt.1 1e-12
$APPDIR/apptest_cmp GOLD-field_3d_cell_f1_r1_nc.txt.2 field_3d_cell_f1_r1_nc.txt.2 1e-12
$APPDIR/apptest_cmp GOLD-field_3d_cell_f1_r1_nc.txt.3 field_3d_cell_f1_r1_nc.txt.3 1e-12
