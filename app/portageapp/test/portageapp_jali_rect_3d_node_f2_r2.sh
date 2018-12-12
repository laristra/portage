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

# 3D, 2nd order accurate remap of quadratic node-centered function

# SERIAL RUN, conforming meshes

mpirun -np 1 $TESTAPPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--conformal=y \
--entity_kind=node --field="x*x+y*y+z*z" \
--remap_order=2 \
--results_file="jali_rect_3d_node_f2_r2.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_node_f2_r2.txt jali_rect_3d_node_f2_r2.txt 1e-12


# PARALLEL RUN, conforming meshes (Note, reverse_ranks is false here)  

mpirun -np 4 $TESTAPPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--reverse_ranks=n --conformal=y \
--entity_kind=node --field="x*x+y*y+z*z" \
--remap_order=2 \
--results_file="jali_rect_3d_node_f2_r2.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_node_f2_r2.txt.0 jali_rect_3d_node_f2_r2.txt.0 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_node_f2_r2.txt.1 jali_rect_3d_node_f2_r2.txt.1 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_node_f2_r2.txt.2 jali_rect_3d_node_f2_r2.txt.2 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_node_f2_r2.txt.3 jali_rect_3d_node_f2_r2.txt.3 1e-12


# PARALLEL RUN non-conforming meshes (Note, reverse_ranks is false here)  

mpirun -np 4 $TESTAPPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--reverse_ranks=n --conformal=n \
--entity_kind=node --field="x*x+y*y+z*z" \
--remap_order=2 \
--results_file="jali_rect_3d_node_f2_r2_nc.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_node_f2_r2_nc.txt.0 jali_rect_3d_node_f2_r2_nc.txt.0 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_node_f2_r2_nc.txt.1 jali_rect_3d_node_f2_r2_nc.txt.1 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_node_f2_r2_nc.txt.2 jali_rect_3d_node_f2_r2_nc.txt.2 1e-12
$CMPAPPDIR/apptest_cmp GOLD_jali_rect_3d_node_f2_r2_nc.txt.3 jali_rect_3d_node_f2_r2_nc.txt.3 1e-12
