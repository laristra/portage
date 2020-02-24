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
FILENAME=field_$1_$2_$3_$4_$5_$6.txt

TOLERANCE_2D=5.e-11
TOLERANCE_3D=5.e-9

rm -f ${FILENAME}*

# SERIAL RUN
mpirun -np 1 $TESTAPPDIR/swept_face_demo \
  --dim=$2 \
  --ncells=$3 \
  --field=$4 \
  --remap_order=$5 \
  --field_filename=${FILENAME} \
  --simple=$6 \
  --ntimesteps=1
 
# PARALLEL RUN
mpirun -np $1 $TESTAPPDIR/swept_face_demo \
  --dim=$2 \
  --ncells=$3 \
  --field_expression=$4 \
  --remap_order=$5 \
  --field_filename=${FILENAME} \
  --simple=$6 \
  --ntimesteps=1

# COMPARE  
if [[ "$2" -eq "2" ]]
then
  TOLERANCE=${TOLERANCE_2D}
else 
  TOLERANCE=${TOLERANCE_3D}
fi

$CMPAPPDIR/distributed_cmp ${FILENAME} ${TOLERANCE}

