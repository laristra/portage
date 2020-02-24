#!/bin/bash
: <<'END'
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
END

# Exit on error
#set -e

# Echo each command
set -x

# This should probably be changed at some time
# In particular $4 is a field expression and can contain non alphanumeric
# characters that can be a problem for filenames, but at the moment it works
# and avoids file clobbers in a human readable way as opposed to hashes.
FILE="field_rank_$1_dim_$2_cells_$3_field_$4_order_$5"
TOLERANCE_2D=5.e-11
TOLERANCE_3D=5.e-9

rm -f ${FILE}*

# SERIAL RUN
mpirun -np 1 $TESTAPPDIR/swept_face_demo \
  --dim=$2 \
  --ncells=$3 \
  --field=$4 \
  --remap_order=$5 \
  --ntimesteps=4 \
  --scale_by=10 \
  --output_meshes=false \
  --result_file="${FILE}" \
  --keep_source=false \
  --simple=true

# PARALLEL RUN
mpirun -np $1 $TESTAPPDIR/swept_face_demo \
  --dim=$2 \
  --ncells=$3 \
  --field=$4 \
  --remap_order=$5 \
  --ntimesteps=4 \
  --scale_by=10 \
  --output_meshes=false \
  --result_file="${FILE}" \
  --keep_source=false \
  --simple=true

# COMPARE  
if [[ "$2" -eq "2" ]]; then
  TOLERANCE=${TOLERANCE_2D}
else 
  TOLERANCE=${TOLERANCE_3D}
fi

$CMPAPPDIR/distributed_cmp ${FILE} ${TOLERANCE}

