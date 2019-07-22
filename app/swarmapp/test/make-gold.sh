#!/bin/bash
problems="0 1 2 3 4 5 6 7 8 9"
nsrc="23"
ntgt="19"
center="0 1"
for cent in $center; do 
  for prob in $problems; do 
    echo "running problem " $prob " with center " $cent
    swarmapp $prob $nsrc $ntgt 0 12345678 1.5 $cent 2>&1 > gold-output-$prob-$cent
    mv -f outfield-$prob-$cent.csv gold-outfield-$prob-$cent.csv
    rm -f infield$prob.{csv}
  done
done

problems="10 11 12 13 14 15 16 17 18 19"
nsrc="9"
ntgt="7"
center="0 1"
for cent in $center; do 
  for prob in $problems; do 
    echo "running problem " $prob " with center " $cent
    swarmapp $prob $nsrc $ntgt 0 12345678 1.5 $cent 2>&1 > gold-output-$prob-$cent
    mv -f outfield-$prob-$cent.csv gold-outfield-$prob-$cent.csv
    rm -f infield$prob.{csv}
  done
done
