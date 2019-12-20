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

# 2nd order momentum CCH remapwil limiters of quadratic density and velocity
${RUN_COMMAND} $TESTAPPDIR/momentumapp 8 6 2 1 2 2

# Compare the values for the field
$CMPAPPDIR/apptest_cmp errors_gold1.txt errors1.txt 1e-10 1e-4
