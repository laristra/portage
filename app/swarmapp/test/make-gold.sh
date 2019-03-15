#!/bin/bash
problems="0 1 2 3 4 5 6 7 8 9"
nsrc="23"
ntgt="19"
distribution="0 1"
for dist in $distribution; do 
  for prob in $problems; do 
    echo "running problem " $prob " in distribution " $dist
    swarmapp $prob $nsrc $ntgt $dist 12345678 2>&1 > gold-output-$prob-$dist
    mv -f outfield-$prob-$dist.csv gold-outfield-$prob-$dist.csv
    rm -f infield$prob.{csv}
  done
done

problems="10 11 12 13 14 15 16 17 18 19"
nsrc="9"
ntgt="7"
distribution="0 1"
for dist in $distribution; do 
  for prob in $problems; do 
    echo "running problem " $prob " in distribution " $dist
    swarmapp $prob $nsrc $ntgt $dist 12345678 2>&1 > gold-output-$prob-$dist
    mv -f outfield-$prob-$dist.csv gold-outfield-$prob-$dist.csv
    rm -f infield$prob.{csv}
  done
done
