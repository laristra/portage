#!/bin/bash

# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
# https://github.com/laristra/portage/blob/master/LICENSE

set -e
set -x

suffix=${1}
epsilon=1.e-10

# run test
mpirun -np 4 ${ROOT_DIR}/part-remap "${ROOT_DIR}/test/input_${suffix}.json"

# compare values with related gold file
${COMPARE} \
  "gold_${suffix}_temperature.dat" \
  "data_${suffix}_temperature.dat" \
  ${epsilon}

if [[ ! "${suffix}" =~ "blocks" ]]; then
  ${COMPARE} \
    "gold_${suffix}_density.dat" \
    "data_${suffix}_density.dat" \
    ${epsilon}
fi