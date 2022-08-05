#!/bin/python

# Note:		need indexed bam file (samtools index in.bam)
# Usage:	python bamfilter.py in.bam out.bam

import sys
import pysam

if len(sys.argv) > 2:
    inputfile = sys.argv[1]
    outputfile = sys.argv[2]
else:
    inputfile = snakemake.input[0]
    outputfile = snakemake.output[0]

samfile = pysam.AlignmentFile(inputfile, "rb")
outfile = pysam.AlignmentFile(outputfile, "wb", template=samfile)

E_coli_len = 5132068


def old_filt():
    for read in samfile.fetch('NC_011750.1'):
        if read.infer_query_length() != None and read.is_proper_pair:
            len_ops, num_ops = read.get_cigar_stats()
            if (sum(len_ops[0:2])/read.infer_read_length() > 0.7) and \
                    (len_ops[10]/sum(len_ops[0:3]) < 0.05):
                outfile.write(read)
    outfile.close()
    samfile.close()


def loose_filt():
    for read in samfile.fetch("NC_011750.1"):
        if read.is_proper_pair:
            len_ops, num_ops = read.get_cigar_stats()
            if (len_ops[10]/sum(len_ops[0:3]) < 0.1):
                outfile.write(read)
    outfile.close()
    samfile.close()


def no_filt():
    for read in samfile.fetch("NC_011750.1"):
        outfile.write(read)
    outfile.close()
    samfile.close()


if __name__ == "__main__":
    if snakemake.params[0] == "-s":
        old_filt()
    elif snakemake.params[0] == "-l":
        loose_filt()
    else:
        no_filt()
