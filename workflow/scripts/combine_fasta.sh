#!/bin/bash

input=$1 # takes a csv?
std_ref=$3 # standard reference file including protocol. quite cool
# sample=$2
output=$2

keyword="$4"
close=$5
threshold_hi=5 #$6
threshold_lo=0.5

utq_t=0.01 # threshold for unique to query
fmo_t=0.1 # threshold for f_match_orig

abund_species_record="$output.txt"

meta="gtdb/gtdb_metadata.csv"
printf "" > $output
printf "" > $abund_species_record
cnt=0
total_cnt=0 
while IFS="," read -r r0 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 r14 #(($cnt < 10)) && (($total_cnt < 25)) && 
do
  contig_name="$r9"
  name=${contig_name:1} # remove the suffix of _1 (or others)
  
  unique_to_query=$r4 # f_unique_to_query
  match_orig=$r11
  # echo $r0 $r1 $r2 $r3
  # echo $name $unique_to_query $match_orig
	file="$(awk -F, '$1 ~ name {print $2}' name="$name" $meta)"
  ani="$(awk -F, '$4 ~ name {print $8}' name="$name" $close)"
  if [[ $ani == "" ]]; then
    # echo $name $ani
    continue
  fi 
  ani_threshold_lo=0.75
  ani_threshold_hi=0.9
  # turn this on if using UHGG
  # species="$(awk '$1 ~ name {print $19,$20}' name="$name" $meta)"
  # if [[ "$name" == *"Escherichia coli"* ]]; then
  #   echo "$name" $cov $ani
  # fi
  # echo $name $file
  
  if (( $(echo "$unique_to_query > $utq_t" | tr -d $'\r' | bc -l) ))  && (( $(echo "$match_origin >= $fmo_t" | tr -d $'\r' | bc -l) )); then
  # if [[ ($(echo "($cov > $threshold_hi)" | tr -d $'\r' | bc -l) -eq 1) \
  #       && ($(echo "($ani > $ani_threshold_lo)" | tr -d $'\r' | bc -l) -eq 1) ]] \
  #       || [[ ($(echo "($ani > $ani_threshold_hi)" | tr -d $'\r' | bc -l) -eq 1) \
  #           && ($(echo "($cov > $threshold_lo)" | tr -d $'\r' | bc -l) -eq 1) ]]; then
    # if [[ $name == *$keyword* ]]; then
    #   continue
    # fi
    echo $name
    # printf "$ani\n"
    echo $name >> $abund_species_record
    # try string replacement there ig
    zcat "gtdb/${file}" >> $output
    cnt=$((cnt+1))
  fi
  total_cnt=$((total_cnt+1))
done < <(tail -n +2 "$input")
# awk '/^>/{f=!d[$1];d[$1]=1}f' $output > $2
# cat $std_ref >> $output
# rm $output
echo "Done"
