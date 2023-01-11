#!/bin/bash

input=$1
std_ref=$3
output=$2

keyword="$4"

echo $keyword

meta="gtdb/metadata.csv"
printf "" > $output
cnt=0
while (($cnt < 5)) && IFS="," read -r r0 r1 r2 r3 r4 r5 r6 r7
do
  name="$r3"
  ani=$r7
	file="$(awk -F',' '$1 ~ name {print $2}' name="$name" $meta)"
  if [[ ($(echo "($ani < 0.935)" | tr -d $'\r' | bc -l) -eq 1) &&  ($(echo "($ani > 0.88)" | tr -d $'\r' | bc -l) -eq 1)  && !("$name" == *"$keyword"*) ]]; then #&& !("$name" == *"$keyword"*)
    echo $name $file $ani >> "gtdb/$input.files.txt"
    zcat "gtdb/$file" >> $output
    cnt=$((cnt+1))
  fi
done < <(tail -n +2 "$input")
# awk '/^>/{f=!d[$1];d[$1]=1}f' $output > $2
cat $std_ref >> $output
echo "Done"
