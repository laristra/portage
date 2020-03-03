#!/bin/bash
: <<'END'
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
END

# Echo each command
set -x

FILE="field_rank_$1_dim_$2_cells_$3_field_$6_order_$5"

rm -f ${FILE}*.txt*

# SERIAL RUN
mpirun -np 1 ${APPTEST} \
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
mpirun -np $1 --oversubscribe ${APPTEST} \
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
${COMPARE} "${FILE}_dist.txt" 5.e-2

