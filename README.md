snpeff need to be installed and added to path but not sure how to do it with snakemake

# testing
``` 
snakemake -c1 --use-conda data/output/ann.vcf
snakemake -c1 --use-conda data/output/coverage.txt
snakemake -c1 --use-conda data/output/coverage.summary
snakemake -c1 --use-conda data/output/ec.bam.pileup.gz
```
running `/data/clean` cleans the output files from tests

## Before running
change in `config.json` the directory of samples
