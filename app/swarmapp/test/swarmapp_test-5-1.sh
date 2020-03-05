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
prob=5
id="$prob-1"

# 2d remap
${RUN_COMMAND} ${TESTAPPDIR}/swarmapp $prob 23 19 0 12345678 1.5 1 2>&1 > output-$id

# Compare the values for the field
../swarmcompare gold-outfield-$id.csv outfield-$id.csv
result1=$?

# Compare the values of the error
grep 'L2 NORM OF ERROR' gold-output-$id > line1-$id
grep 'L2 NORM OF ERROR' output-$id > line2-$id
v1=$(cut -f 2 -d '=' < line1-$id | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
v2=$(cut -f 2 -d '=' < line2-$id | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
result2=$(echo scale=16\; v3=$v1-$v2\; if \(v3\<0\) v3=-v3\; if \(v3\<10^-12\) 1 else 0| bc)

test $result1 == 0 -a $result2 == 1
