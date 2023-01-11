#!/usr/env python

import pysam
import sys
import pandas as pd

def parse_vcf(fname):
    yield from pysam.VariantFile(fname).fetch()

def extract_vcf_stat(fname):
    return pd.DataFrame((
        {
            'contig': record.contig,
            'position': record.pos,
            'ref': record.alleles[0],
            'alt': alt,
            'SB': record.info["SB"],
            'DP': record.info["DP"],
            'AF': str(record.info["AF"]).split(',')[idx],
        }
        for record in parse_vcf(fname)
        for idx, alt in enumerate(record.alts)
    )
    )

def filter_vcf(vcfin):
  vcf_in_df = extract_vcf_stat(vcfin)
  mean_cov = vcf_in_df['DP'].mean()
  std_cov = vcf_in_df['DP'].std()
  mean_sb = vcf_in_df['SB'].mean()
  std_sb = vcf_in_df['SB'].std()
  vcf_in = pysam.VariantFile(vcfin)
#   print(mean_cov,mean_sb)
  print(vcf_in.header,end="")
  for record in vcf_in.fetch():
    if ((record.info['DP'] < mean_cov + 2 * std_cov) and \
        (record.info['SB'] < mean_sb + 3 * std_sb)):
        print(record,end="")
        pass
    

if __name__ == "__main__":
  filter_vcf(sys.argv[1])