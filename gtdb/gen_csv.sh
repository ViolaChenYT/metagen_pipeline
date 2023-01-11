#!/bin/bash
output="metadata.csv"
printf "name,genome_filename,protein_filename\n" > $output
for f in `find "./gtdb_genomes_reps/" -type f`
do
	line="$(zcat $f | head -n 1)"
	# s="$(grep -oP '(?<= ).*?(?=,)' <<< $line)"
	s="$(echo $line | sed 's/,/ /g')"
	printf "$s," >> $output
	printf "$f,\n" >> $output
done
