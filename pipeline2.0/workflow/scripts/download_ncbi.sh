# @Author: jsgounot
# @Date:   2023-03-10 13:49:47
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-10 14:23:48

# Download an assembly NCBI sequence
# Take either refseq or insc depending for the identifier
# Switch if one fail, can happen when an ID has been removed from the db
# example: https://www.ncbi.nlm.nih.gov/assembly/GCA_902498915.1

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
fi