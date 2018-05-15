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

# if we ever need to do this for real, then the paths below will need to be modified
# app/simple_mesh_app/simple_mesh_mm_app --nsourcecells=3 --ntargetcells=5 \
# --material_file=app/simple_mesh_app/test/triple_point_3x3_matdata.dat
