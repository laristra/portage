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

# SERIAL
mpirun -np 1 ${TESTAPPDIR}/portageapp_jali \
    --source_file=${DATA_DIR}/coarse_poly2D.exo \
    --target_file=${DATA_DIR}/fine_poly2D.exo \
    --entity_kind=cell \
    --field="x + y" \
    --remap_order=2 \
    --limiter=barth_jespersen \
    --results_file="jali_unstruc_2d_cell_f1_r8.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_unstruc_2d_cell_f1_r8.txt jali_unstruc_2d_cell_f1_r8.txt 1e-9

