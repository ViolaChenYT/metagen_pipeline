snpeff need to be installed and added to path but not sure how to do it with snakemake

# testing
```{
snakemake -c1 --use-conda data/output/ann.vcf
snakemake -c1 --use-conda data/output/coverage.txt
snakemake -c1 --use-conda data/output/coverage.summary
}
```