#!/bin/python

# Note:		need indexed bam file (samtools index in.bam)
# Usage:	python bamfilter.py in.bam out.bam

import sys, pysam

inputfile = snakemake.input[0]
outputfile = snakemake.output[0]

samfile = pysam.AlignmentFile(inputfile, "rb")
outfile = pysam.AlignmentFile(outputfile, "wb", template=samfile)

def old_filt():
	for read in samfile.fetch():
		if read.infer_query_length() != None and read.is_proper_pair:
			len_ops, num_ops = read.get_cigar_stats()
			if sum(len_ops[0:2])/read.infer_read_length() > 0.7 and len_ops[10]/sum(len_ops[0:3]) < 0.05:
				outfile.write(read)			
	outfile.close()
	samfile.close()

def loose_filt():
	for read in samfile.fetch():
		if read.infer_query_length() != None and read.is_proper_pair:
			len_ops, num_ops = read.get_cigar_stats()
			if sum(len_ops[0:2])/read.infer_read_length() > 0.5 and len_ops[10]/sum(len_ops[0:3]) < 0.05:
				outfile.write(read)			
	outfile.close()
	samfile.close()

if __name__ == "__main__":
	if len(sys.argv) < 2:
		old_filt()
	else:
		print("lose filter")
		loose_filt()
		
