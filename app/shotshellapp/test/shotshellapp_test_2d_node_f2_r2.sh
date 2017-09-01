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

DATA_DIR=$APPDIR/../test/wonton

# SERIAL

mpirun -np 1 ${APPDIR}/shotshellapp_jali \
    --source=${DATA_DIR}/shotshell.exo \
    --target=${DATA_DIR}/shotshell-v.exo \
    --entity_kind=node \
    --field_order=2 \
    --remap_order=2 \
    --output_results=y

# Compare the values for the field
$APPDIR/apptest_cmp GOLD-field_2d_node_f2_r2.txt field_2d_node_f2_r2.txt 1e-9

# PARALLEL - DOES NOT WORK

# mpirun -np 4 ${APPDIR}/shotshellapp_jali \
#     --source=${DATA_DIR}/shotshell.exo \
#     --target=${DATA_DIR}/shotshell-v.exo \
#     --entity_kind=node \
#     --field_order=2 \
#     --remap_order=2 \
#     --output_results=y

# # Compare the values for the field
# $APPDIR/apptest_cmp GOLD-field_2d_node_f2_r2.txt.0 field_2d_node_f2_r2.txt.0 1e-9
# $APPDIR/apptest_cmp GOLD-field_2d_node_f2_r2.txt.1 field_2d_node_f2_r2.txt.1 1e-9
# $APPDIR/apptest_cmp GOLD-field_2d_node_f2_r2.txt.2 field_2d_node_f2_r2.txt.2 1e-9
# $APPDIR/apptest_cmp GOLD-field_2d_node_f2_r2.txt.3 field_2d_node_f2_r2.txt.3 1e-9
