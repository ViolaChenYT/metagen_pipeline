snpeff need to be installed and added to path but not sure how to do it with snakemake

# Pre-requisite
- sourmash
- gtdb database: download from https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/gtdb_genomes_reps_r207.tar.gz and unzip into a directory `gtdb/gtdb_genome_reps/` would allow the script to use it automatically
- run `sourmash sketch fromfile gtdb_metadata.csv -p dna -o gtdb_sketch` to generate the searchable sketch database for sourmash

# testing
``` 
snakemake -c1 --use-conda data/output/ann.vcf
snakemake -c1 --use-conda data/output/coverage.txt
snakemake -c1 --use-conda data/output/coverage.summary
snakemake -c1 --use-conda data/output/ec.bam.pileup.gz
```
running `./run clean` cleans the output files from tests

running `./run` automatically runs current Snakefile with `config.json`

## Before running
change in `config.json` the directory of samples
