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
    --source_file=${DATA_DIR}/sector5.exo \
    --target_file=${DATA_DIR}/sector6.exo \
    --entity_kind=node \
    --field="x+y" \
    --remap_order=2 \
    --results_file="jali_radial_3d_node_f1_r2.txt"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_jali_radial_3d_node_f1_r2.txt jali_radial_3d_node_f1_r2.txt 1e-9

