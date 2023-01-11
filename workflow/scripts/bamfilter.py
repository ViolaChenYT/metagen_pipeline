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

def get_contigs():
    ref = snakemake.params[1]  # the reference
    contig_names = []
    with open(ref, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == ">":
                contig = (line.split()[0])[1:]
                contig_names.append(contig.strip())
                # print(contig)
        f.close()
    contig_names.sort()
    return contig_names

def old_filt(contigs):
    for contig in contigs:
        print(contig, end=", ")
        cnt = 0
        for read in samfile.fetch(contig):
            if read.infer_query_length() != None and read.is_proper_pair:
                len_ops, num_ops = read.get_cigar_stats()
                if (sum(len_ops[0:2])/read.infer_read_length() > 0.9) and (len_ops[10]/sum(len_ops[0:3]) < 0.05):
                    outfile.write(read)
                    cnt += 1
        print(cnt)
#   print(contig_set)
#   print(contigs)
    outfile.close()
    samfile.close()


def loose_filt(contigs):
    for contig in contigs:
        for read in samfile.fetch(contig):
            if read.is_proper_pair:
                len_ops, num_ops = read.get_cigar_stats()
                if (sum(len_ops[0:2])/read.infer_read_length() > 0.7) and (len_ops[10]/sum(len_ops[0:3]) < 0.05):
                    if read.is_secondary or read.is_supplementary:
                        pass
                    else:
                        outfile.write(read)
    outfile.close()
    samfile.close()


def no_filt(contigs):
    # print(contigs)
    for contig in contigs:
        cnt = 0
        print(contig, end=", ")
        for read in samfile.fetch(contig):
            outfile.write(read)
            cnt += 1
            # else:
        print(cnt)
    outfile.close()
    samfile.close()


if __name__ == "__main__":
    contigs = get_contigs()
    if snakemake.params[0] == "-s":
        old_filt(contigs)
    elif snakemake.params[0] == "-l":
        loose_filt(contigs)
    else:
        no_filt(contigs)
