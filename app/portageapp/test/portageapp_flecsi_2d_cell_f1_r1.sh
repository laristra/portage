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

# 2D, 1st order accurate remap of linear cell-centered function

# SERIAL RUN

mpirun -np 1 $APPDIR/portageapp_flecsi \
       6 6 1 y # nx, ny, order, dump_output

# Compare the values for the field
$APPDIR/apptest_cmp GOLD_flecsi_rect_2d_cell_f1_r1.txt flecsi_field_2d_cell_f1_r1.txt 1e-12

