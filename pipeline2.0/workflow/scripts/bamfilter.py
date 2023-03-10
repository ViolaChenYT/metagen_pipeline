# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-03-07 14:35:28
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-07 17:52:32

import pysam

# First pass
bamfile = snakemake.input[0]
bamfile = pysam.AlignmentFile(bamfile, "rb")

threshold_coverage = snakemake.params.tcov
threshold_similarity = snakemake.params.tsim
found = dict()

for read in bamfile.fetch(until_eof=True):
    iql = read.infer_query_length()
    ipp = read.is_proper_pair

    if iql != None and ipp:
        len_ops, num_ops = read.get_cigar_stats()
        mid = sum(len_ops[0:3])

        # (M + I + D / read_length) and NM / (M + I + D)
        cov_test  = mid / iql >= threshold_coverage
        sim_value = len_ops[10] / mid 
        sim_test  = sim_value <= threshold_similarity
         
        if cov_test and sim_test:
            name = read.query_name
            contig = read.reference_name

            best_sim, names = found.get(name, (101, None))
            if sim_value < best_sim:
                found[name] = (sim_value, set((contig,)))
            elif sim_value == best_sim:
                names.add(contig)

bamfile.close()

# Second pass
bamfile = snakemake.input[0]
bamfile = pysam.AlignmentFile(bamfile, "rb")

outfile = snakemake.output[0]
outfile = pysam.AlignmentFile(outfile, "wb", template=bamfile)

reads_count = {
    'total': 0,
    'first': 0,
    'second': 0,
    'last': 0
}

for read in bamfile.fetch(until_eof=True):
    iql = read.infer_query_length()
    ipp = read.is_proper_pair
    reads_count['total'] += 1

    if iql != None and ipp:
        len_ops, num_ops = read.get_cigar_stats()
        mid = sum(len_ops[0:3])
        reads_count['first'] += 1

        # (M + I + D / read_length) and NM / (M + I + D)
        cov_test  = mid / iql >= threshold_coverage
        sim_value = len_ops[10] / mid 
        sim_test  = sim_value <= threshold_similarity    
        
        if cov_test and sim_test:
            name = read.query_name
            contig = read.reference_name
            reads_count['second'] += 1

            if contig in found[name][1]:
                outfile.write(read)
                reads_count['last'] += 1

bamfile.close()
outfile.close()

logfile = snakemake.log[0]
with open(logfile, 'w') as f:
    for element, value in reads_count.items():
        prc = value * 100 / reads_count['total']
        f.write(f'{element} {value} {prc:.2f}%\n')
