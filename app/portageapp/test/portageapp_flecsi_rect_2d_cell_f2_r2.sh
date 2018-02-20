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

# 2D, 2nd order accurate remap of quadratic cell-centered function

# SERIAL RUN

mpirun -np 1 $TESTAPPDIR/portageapp_flecsi \
       6 6 2 y # nx, ny, order, dump_output

# Compare the values for the field
$CMPAPPDIR/apptest_cmp GOLD_flecsi_rect_2d_cell_f2_r2.txt flecsi_field_2d_cell_f2_r2.txt 1e-12
