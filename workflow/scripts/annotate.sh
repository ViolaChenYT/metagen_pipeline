#!/bin/bash

# echo "usage: <species> | <ID> | ..."
if [[ $1 == "EC" ]]; then
	sed -i 's/gi|218698419|ref|NC_011750.1|/Chromosome/g' $2
	snpeff Escherichia_coli_iai39 $2 > $3

elif [[ $1 == "KP" ]]; then
	sed -i 's/NC_016847.1/pKPHS5/g' $2
	sed -i 's/NC_016846.1/pKPHS2/g' $2
	sed -i 's/NC_016845.1/Chromosome/g' $2
	sed -i 's/NC_016840.1/pKPHS4/g' $2
	sed -i 's/NC_016839.1/pKPHS3/g' $2

	snpeff Klebsiella_pneumoniae_subsp_pneumoniae_hs11286 $2 > $3

else
	echo "species not recognized, currently only recognize KP and EC"
	echo "exiting..."
fi

