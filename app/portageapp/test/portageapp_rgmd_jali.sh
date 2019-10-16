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

FILENAME=field.txt
TOLERANCE=3.e-13

rm -f ${FILENAME}*

# SERIAL RUN
mpirun -np 1 $TESTAPPDIR/portageapp_rgmd_jali \
  --problem=$2 \
  --dim=$3 \
  --nsourcecells=$4 \
  --ntargetcells=$5 \
  --material_fields=$6 \
  --remap_order=$7 \
  --field_filename=${FILENAME}
 
# PARALLEL RUN
mpirun -np $1 $TESTAPPDIR/portageapp_rgmd_jali \
  --problem=$2 \
  --dim=$3 \
  --nsourcecells=$4 \
  --ntargetcells=$5 \
  --material_fields=$6 \
  --remap_order=$7 \
  --field_filename=${FILENAME}

# COMPARE  
$CMPAPPDIR/distributed_cmp ${FILENAME} ${TOLERANCE}

