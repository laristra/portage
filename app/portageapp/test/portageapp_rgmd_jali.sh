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

FILENAME=field_$1_$2_$3_$4_$5_$6_$7.txt

TOLERANCE_2D=5.e-11
TOLERANCE_3D=5.e-9

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
if [[ "$3" -eq "2" ]]
then
  TOLERANCE=${TOLERANCE_2D}
else 
  TOLERANCE=${TOLERANCE_3D}
fi

$CMPAPPDIR/distributed_cmp ${FILENAME} ${TOLERANCE}

