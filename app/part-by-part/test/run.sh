#!/bin/bash

# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
# https://github.com/laristra/portage/blob/master/LICENSE

set -e
set -x

test_case=${1}
epsilon=1.e-10

# run test
mpirun -np 4 ${ROOT_DIR}/part-remap "${ROOT_DIR}/test/input_${test_case}.json"

# compare values with gold file
${COMPARE} \
  "gold_${test_case}_temperature.dat" \
  "data_${test_case}_temperature.dat" \
  ${epsilon}

if [[ ! "${test_case}" =~ "blocks" ]]; then
  ${COMPARE} \
    "gold_${test_case}_density.dat" \
    "data_${test_case}_density.dat" \
    ${epsilon}
fi