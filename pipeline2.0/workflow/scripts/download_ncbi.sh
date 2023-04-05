#!/bin/bash

if [[ $1 == GCF* ]] ;
then
	name1="assembly_nuccore_refseq"
	name2="assembly_nuccore_insdc"
else
	name1="assembly_nuccore_insdc"
	name2="assembly_nuccore_refseq"
fi

#echo $1 $2 $name1

esearch -db assembly -query $1 | elink -target nuccore -name $name1 | efetch -format fasta > $2

if [[ ! -s $2 ]] ;
then
	echo "First download fail for $1, switch to db name $name2"
	esearch -db assembly -query $1 | elink -target nuccore -name $name2 | efetch -format fasta > $2