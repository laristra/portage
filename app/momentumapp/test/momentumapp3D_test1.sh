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

# 2nd order momentum CCH remap of quadratic density and velocity
${RUN_COMMAND} $TESTAPPDIR/momentumapp3D 10 10 10 2 1  "1+x*y*z" "0.5-y" "-0.5+x + (z/3)" "z*(1-z)"

# Compare the values for the field
$CMPAPPDIR/apptest_cmp errors3D_gold1.txt errors3D_1.txt 1e-10 1e-4
