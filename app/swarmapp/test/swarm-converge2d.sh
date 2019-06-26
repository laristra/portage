#!/bin/bash

rm -f results.csv
touch results.csv
echo "{" >> results.csv

if [[ $1 == "-h" ]]; then
    echo "swarm-convger2d.sh distribution example hscale center seed"
    exit 0
fi

dist=$1
example=$2
hscale=$3
center=$4
seed=$5

for i in 8 16 32 64 128 256 512; do
    nsrc=$((i+i/8))
    echo doing $example $nsrc "x" $i $dist $hscale $center $seed
    swarmapp $example $nsrc $i $dist $seed $hscale $center 2>&1 > output-$i
    mv infield$example.csv infield$example-$i.csv
    mv outfield-$example-$dist.csv outfield-$example-$dist-$i.csv
    echo "{" $(grep 'smoothing length' output-$i | cut -d ' ' -f 4) "," >> results.csv
    echo     $(grep 'Linf NORM OF ERROR' output-$i | cut -d ' ' -f 6) "," >> results.csv
    echo     $(grep 'L1 NORM OF ERROR' output-$i | cut -d ' ' -f 6) "," >> results.csv
    echo     $(grep 'L2 NORM OF ERROR' output-$i | cut -d ' ' -f 6) "}," >> results.csv
done

echo "}" >> results.csv
