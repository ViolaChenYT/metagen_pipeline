#!/bin/bash

sourmash sketch dna ../refs/Escherichia_coli_iai39.fasta
wget https://osf.io/748ew/download # for k-mer size 31
mv download gtdb_all_genomes.stb.gz
sourmash search Escherichia_coli_iai39.fasta.sig gtdb_all_genomes.stb.gz
