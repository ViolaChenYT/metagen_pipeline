#!/bin/bash

input=$1
output=$2
output="../$output"
meta="metadata.csv"
cd "gtdb"
while IFS="," read -r r0 r1 r2 r3 r4 r5 r6 r7 r8
do
  name="$r3"
  ani="$r7"
  if [[ $(echo "($ani > 0.88) && ($ani < 0.94)" | tr -d $'\r' | bc -l) -eq 1 ]]; then
    file="$(awk -F',' '$1 ~ name {print $2}' name="$name" $meta)"
    echo $name $file $ani
    zcat $file >> $output
  fi
done < <(tail -n +2 "../$input")
echo "Done"
