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

# PARALLEL

mpirun -np 4 ${TESTAPPDIR}/portageapp_jali \
    --dim=3 \
    --conformal=y \
    --nsourcecells=6 \
    --mesh_min=-0.5 \
    --mesh_max=0.5 \
    --target_file=${DATA_DIR}/3dvoro-voronoized-10k.exo \
    --entity_kind=cell \
    --field="atan(1e50*x)/atan(1e50)" \
    --remap_order=2 \
    --limiter=barth_jespersen \
    --results_file="jali_unstruc_3d_cell_f8_r8.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_unstruc_3d_cell_f8_r8.txt.0 jali_unstruc_3d_cell_f8_r8.txt.0 1e-9
$CMPAPPDIR/apptest_cmp GOLD_jali_unstruc_3d_cell_f8_r8.txt.1 jali_unstruc_3d_cell_f8_r8.txt.1 1e-9
$CMPAPPDIR/apptest_cmp GOLD_jali_unstruc_3d_cell_f8_r8.txt.2 jali_unstruc_3d_cell_f8_r8.txt.2 1e-9
$CMPAPPDIR/apptest_cmp GOLD_jali_unstruc_3d_cell_f8_r8.txt.3 jali_unstruc_3d_cell_f8_r8.txt.3 1e-9
