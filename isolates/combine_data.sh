#!/bin/bash

output_prefix=$1
input=$2
echo "usage: ./combine_data.sh <output_prefix> <metadata file>"

output1="${output_prefix}_1.fastq.gz"
output2="${output_prefix}_2.fastq.gz"

printf "" > $output1
printf "" > $output2

n_lines=400000
count=0

while IFS=$'  ' read -r r0 r1 r2  #(($cnt < 10)) && (($total_cnt < 25)) && 
do ##  r0: species, r1: file name, r2: url
  if [ ! -f $r1 ]; then
    continue
  fi
  if [[ $r0 == "Escherichia_coli" ]]; then
    if [[ $r1 == *_1.fastq.gz ]]; then
      # cat $r1 >> $output1
      count=$((count+1))
    else 
      # cat $r1 >> $output2
      count=$((count+1))
    fi
  else
    if [[ $r1 == *_1.fastq.gz ]]; then
      zcat $r1 | head -n $n_lines | gzip >> $output1
      count=$((count+1))
    else 
      zcat $r1 | head -n $n_lines | gzip >> $output2
      count=$((count+1))
    fi
  fi  
done < <(tail -n +2 "$input")
echo $count
