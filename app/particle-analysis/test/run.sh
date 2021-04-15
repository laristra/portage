#!/bin/bash

# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
# https://github.com/laristra/portage/blob/master/LICENSE

set -e
set -x

epsilon=1.e-10

mpirun -np 1 ${ROOT_DIR}/particle-analysis "${ROOT_DIR}/test/gaussian-radial.json"

${COMPARE} "gold_source.dat" "source.dat" ${epsilon}
${COMPARE} "gold_exact.dat"  "exact.dat"  ${epsilon}
${COMPARE} "gold_remap.dat"  "remap.dat"  ${epsilon}