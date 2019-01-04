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

# 2d 0-th order remap of 0-th order function distribution 0 (grid)
${RUN_COMMAND} ${TESTAPPDIR}/swarmapp 6 23 19 0 12345678 2>&1 > output-6-0

# Compare the values for the field
python compare.py gold-outfield-6-0.csv outfield-6-0.csv 0
