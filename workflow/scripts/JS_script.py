import pysam
import math
from collections import defaultdict
import pandas as pd
import numpy as np
import sys
from os.path import exists
'''
Originally written by Jean-Sebastien
Significantly modified by Viola Chen
July 22 2022
'''


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


def parse_bam(fname, positions):
    bamfile = pysam.AlignmentFile(fname, "rb")
    bases = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    positions = {contig: {pos - 1 for pos in cpos}
                 for contig, cpos in positions.items()}

    for col in bamfile.pileup():

        if col.pos not in positions[col.reference_name]:
            continue

        col.set_min_base_quality(0)

        for read in col.pileups:
            if read.is_del == 1:
                base = 'D'
            else:
                base = read.alignment.query_sequence[read.query_position_or_next]
            bases[col.reference_name][col.pos][base] += 1

        if col.pos == 0:
            break

    bases = pd.DataFrame((
        {'contig': contig, 'position': pos + 1, ** bases}
        for contig, cpos in bases.items()
        for pos, bases in cpos.items()
    ))

    bnames = [column for column in bases.columns if column not in (
        'contig', 'position')]
    bases[bnames] = bases[bnames].fillna(0).astype(int)
    bases['DP_bam'] = bases[bnames].sum(axis=1)

    return bases


def assign_row_kind(row):
    if row['iso'] and row['meta']:
        return 'TP'
    elif row['iso']:
        return 'FN'
    elif row['meta']:
        return 'FP'
    else:
        raise Exception(row)


def compare_iso_meta():
    # isoname = "WEB1015"
    # iso = f"output/{isoname}/Escherichia_coli_iai39.nofilt.lofreq.vcf"
    iso = snakemake.input.iso
    iso = extract_vcf_stat(iso)

    # metaname = "MBH086"
    # meta = f"output/{metaname}/combined.nofilt.lofreq.vcf"
    meta = snakemake.input.meta
    meta = extract_vcf_stat(meta)

    df = iso.merge(meta, on=['contig', 'position', 'ref',
                   'alt'], how='outer', suffixes=('_iso', '_meta'))

    positions = df.groupby('contig')['position'].apply(set).to_dict()
    # meta = f"output/{metaname}/combined.nofilt.bam"
    # meta = snakemake.input.bam
    # baminfo = parse_bam(meta, positions)

    # df = df.merge(baminfo, on=['contig', 'position'], how='left')
    # df['AF_bam'] = df.apply(
    #     lambda row: row[row['alt']] / row['DP_bam'], axis=1)

    df['iso'] = ~ df['AF_iso'].isnull()
    df['meta'] = ~ df['AF_meta'].isnull()
    df['kind'] = df.apply(assign_row_kind, axis=1)

    print(df['kind'].value_counts())

    df['AF_iso'] = df['AF_iso'].astype(float)
    df['AF_meta'] = df['AF_meta'].astype(float)

    outfile = snakemake.output[0]
    df.to_csv(outfile, sep='\t')

    colnames = ['SB_iso', 'DP_iso', 'AF_iso', 'SB_meta',
                'DP_meta', 'AF_meta']
    # print(df)
    print(df.groupby('kind')[colnames].mean())


def parse_instrain(filename):
    df = pd.read_csv(filename, sep="\t", header=0)
    # print(df["scaffold"] == "NC_011750.1")
    df = df[df["scaffold"] == "NC_011750.1"]
    df["DP_meta"] = df["position_coverage"]
    df["contig"] = df["scaffold"].astype(str)
    df["position"] = (df["position"] + 1).astype(str)
    return df


def compare_instrain_iso():
    # iso = "output/WEB1015/Escherichia_coli_iai39.nofilt.lofreq.vcf"
    iso = snakemake.input.iso
    # iso = extract_vcf_stat(iso)
    # positions = iso.groupby('contig')['position'].apply(set).to_dict()
    iso = extract_lst()
    iso["AF"] = 0
    # instrain = "output/MBH082/Escherichia_coli_iai39_profile_IS/output/Escherichia_coli_iai39_profile_IS_SNVs.tsv"
    instrain = snakemake.input.meta
    instrain = parse_instrain(instrain)

    df = iso.merge(instrain, on=["contig", "position"],
                   how="outer", suffixes=('_iso', '_IS'))
    df['iso'] = ~ df['AF'].isnull()
    df["meta"] = ~ df['ref_base'].isnull()
    df['kind'] = df.apply(assign_row_kind, axis=1)
    outfile = snakemake.output[0]
    print("sample:", outfile, "\nt", df['kind'].value_counts())
    
    # outfile = "output/MBH082/single_instrain_iso_nofilt.tsv"
    df.to_csv(outfile, sep='\t')


def extract_lst():
    lst = snakemake.input.iso
    lst = pd.read_csv(lst, sep="\t", header=2, dtype=str)
    # df = lst[lst["orig_position"] == lst["new_position"]]
    lst["position"] = lst["orig_position"].astype(str)
    lst["contig"] = (lst['# seq_id']).str[:11]
    return lst


def compare_simulated():
    meta = snakemake.input.meta
    meta = extract_vcf_stat(meta)
    meta["position"] = meta["position"].astype(str)

    standard = extract_lst()
    df = meta.merge(standard, on=["contig", "position"],
                    how='outer', suffixes=('_std', '_meta'))
    
    print(df)
    df['iso'] = ~ df['# seq_id'].isnull()
    df['meta'] = ~ df['AF'].isnull()
    df['kind'] = df.apply(assign_row_kind, axis=1)
    # print(df)
    print(df['kind'].value_counts())

    outfile = snakemake.output[0]
    df.to_csv(outfile, sep='\t')
    # print(df)


def check_bam_parsing():
    isolate = "WEB1012"
    iso = f"output/{isolate}/Escherichia_coli_iai39.tight_filt.lofreq.vcf"
    df = extract_vcf_stat(iso)
    df['AF'] = df['AF'].astype(float)

    positions = df.groupby('contig')['position'].apply(set).to_dict()

    iso = f"output/{isolate}/Escherichia_coli_iai39.tight_filt.bam"
    baminfo = parse_bam(iso, positions)

    df = df.merge(baminfo, on=['contig', 'position'], how='left')
    df['AF_bam'] = df.apply(
        lambda row: row[row['alt']] / row['DP_bam'], axis=1)

    df['diff'] = df['AF'] - df['AF_bam']
    df = df[df['diff'] > 0.0001]

    print(df[df['diff'] > 0.0001])
    assert df[df['diff'] > 0.0001].empty
    print(df)


if __name__ == "__main__":
    print(snakemake.input.meta)
    if snakemake.params[0] == "meta":
        compare_iso_meta()
    elif snakemake.params[0] == "instrain":
        compare_instrain_iso()
    else:
        compare_simulated()
