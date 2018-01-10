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

mpirun -np 1 ${APPDIR}/portageapp_jali \
    --source_file=${DATA_DIR}/shotshell.exo \
    --target_file=${DATA_DIR}/shotshell-v.exo \
    --entity_kind=cell \
    --field="x+y" \
    --remap_order=2 \
    --results_file="jali_unstruc_2d_cell_f1_r2.txt"

# Compare the values for the field
$APPDIR/apptest_cmp GOLD_jali_unstruc_2d_cell_f1_r2.txt jali_unstruc_2d_cell_f1_r2.txt 1e-9


# PARALLEL

mpirun -np 4 ${APPDIR}/portageapp_jali \
    --source_file=${DATA_DIR}/shotshell.exo \
    --target_file=${DATA_DIR}/shotshell-v.exo \
    --entity_kind=cell \
    --field="x+y" \
    --remap_order=2 \
    --results_file="jali_unstruc_2d_cell_f1_r2.txt"

# Compare the values for the field
$APPDIR/apptest_cmp GOLD_jali_unstruc_2d_cell_f1_r2.txt.0 jali_unstruc_2d_cell_f1_r2.txt.0 1e-9
$APPDIR/apptest_cmp GOLD_jali_unstruc_2d_cell_f1_r2.txt.1 jali_unstruc_2d_cell_f1_r2.txt.1 1e-9
$APPDIR/apptest_cmp GOLD_jali_unstruc_2d_cell_f1_r2.txt.2 jali_unstruc_2d_cell_f1_r2.txt.2 1e-9
$APPDIR/apptest_cmp GOLD_jali_unstruc_2d_cell_f1_r2.txt.3 jali_unstruc_2d_cell_f1_r2.txt.3 1e-9

