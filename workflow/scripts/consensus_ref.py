import sys

ref = sys.argv[1]
vcf = sys.argv[2]
outputfile = sys.argv[3]


def parse_vcf(fname):
    yield from pysam.VariantFile(fname).fetch()


def extract_vcf_stat(fname):
    ''' extract important fields from a VCF file:
        - contig name
        - position
        - reference and alternate alleles
        - strand bias, coverage, allele frequency
    '''
    return pd.DataFrame((
        {
            'contig': record.contig,
            'position': int(record.pos),
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
