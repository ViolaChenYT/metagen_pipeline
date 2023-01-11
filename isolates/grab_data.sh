#!/bin/bash

input=$1

while IFS=$'  ' read -r r0 r1 r2  #(($cnt < 10)) && (($total_cnt < 25)) && 
do
  if [ ! -f "$r1" ] ; then
    wget $r2
    echo $r1 $r2
  fi
done < <(tail -n +2 "$input")
