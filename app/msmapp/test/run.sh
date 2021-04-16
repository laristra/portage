#!/bin/bash

# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
# https://github.com/laristra/portage/blob/master/LICENSE

set -e
set -x

epsilon=1.e-10

mpirun -np 1 ${ROOT_DIR}/msmapp "${ROOT_DIR}/test/example_input"

${COMPARE} "gold_diagnostics.dat" "diagnostics.dat" ${epsilon}
