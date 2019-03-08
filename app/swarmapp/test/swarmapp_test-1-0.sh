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

# serial id for file names
prob=1
id="$prob-0"

# 2d 0-th order remap of 0-th order function distribution 0 (grid)
${RUN_COMMAND} ${TESTAPPDIR}/swarmapp $prob 23 19 0 12345678 2>&1 > output-$id

# Compare the values for the field
../compare gold-outfield-$id.csv outfield-$id.csv
result1=$?

# Compare the values of the error
grep 'L2 NORM OF ERROR' gold-output-$id > line1-$id
grep 'L2 NORM OF ERROR' output-$id > line2-$id
diff line1-$id line2-$id
result2=$?

test $result1 == 0 -a $result2 == 0
