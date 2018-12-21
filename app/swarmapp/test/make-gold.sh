#!/bin/bash
problems="0 1 2 3 4 5 6 7 8 9"
nsrc="23"
ntgt="19"
distribution="0 1"
for dist in $distribution; do 
  for prob in $problems; do 
    echo "running problem " $prob " in distribution " $dist
    swarmapp $prob $nsrc $ntgt $dist 2>&1 > gold-output-$prob-$dist
    mv -f outfield$prob.txt gold-outfield-$prob-$dist.txt
    rm -f outfield$prob.csv
    rm -f infield$prob.txt
  done
done
