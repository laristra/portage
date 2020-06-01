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

# 2nd order momentum SGH remap of quadratic density and velocity
${RUN_COMMAND} $TESTAPPDIR/momentumapp2D_mm 8 6 1 "1+x, 1+y" "1+x*x, 1+x*x" "2*y*y"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp errorsMM_gold0.txt errorsMM_0.txt 1e-10 1e-4
