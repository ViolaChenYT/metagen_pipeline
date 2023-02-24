import pysam
import math
from collections import defaultdict
import pandas as pd
import numpy as np
import sys
from os.path import exists
import subprocess
import matplotlib.pyplot as plt
import json
import seaborn as sns
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

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################

## NEW CODE Feb 22, 2023

def get_isolate_meta():
  metadata_file = "isolates/metadata.txt"
  df = pd.read_csv(metadata_file, sep=" ",header=0)
  def f(filename):
    return filename.split("_")[0]
  df["name"] = df["filename"].apply(f)
  df.index = df["name"]
  return df.to_dict()

def extract_lst(snplst,length=11):
  lst = pd.read_csv(snplst, sep="\t", header=2, dtype=str)
  # print(lst)
    # df = lst[lst["orig_position"] == lst["new_position"]]
  lst["position"] = lst["orig_position"].astype(str)
  lst["contig"] = (lst['# seq_id']).str[:length]
  return lst

def read_vcf(vcf, snp):
  ''' take in unfiltered vcf, fitler based on coverage, and compare to SNP-list to output variant-calling performance'''
  def parse_vcf(fname):
    yield from pysam.VariantFile(fname).fetch()
  my_snps = pd.DataFrame((
        {
            'contig': record.contig,
            'position': str(record.pos),
            'ref': record.alleles[0],
            'alt': alt,
            'SB': record.info["SB"],
            'DP': record.info["DP"],
            'AF': str(record.info["AF"]).split(',')[idx],
        }
      for record in parse_vcf(vcf)
    for idx, alt in enumerate(record.alts) 
  ))
  contiglen = len(my_snps["contig"].iloc[1])
  true_snps = extract_lst(snp, length=contiglen)
  mean_cov =  my_snps["DP"].mean()
  my_snps = my_snps[my_snps["DP"] < 3 * mean_cov]
  # print((true_snps["contig"]), (my_snps["contig"]))
  df = my_snps.merge(true_snps, on=["contig","position"], how="outer",suffixes=("_true","_lofreq"))
  df['true'] = ~ df['# seq_id'].isnull()
  df['lofreq'] = ~ df['AF'].isnull()
  def assign_row_kind(row):
    if row['true'] and row['lofreq']: return 'TP'
    elif row['true']:  return 'FN'
    elif row['lofreq']: return 'FP'
    else: raise Exception(row)
  df['kind'] = df.apply(assign_row_kind, axis=1)
  # print(df)
  d = df['kind'].value_counts()
  if 'FP' not in d:
    d['FP'] = 0
  return d['FN'],d['FP'],d['TP']

def check_read_n_snp(result_dir, methods, fastq, snp=None,identifier='read'):
  bash_command = f"zcat {fastq} | grep '_180_' | wc -l"
  # print(bash_command)
  n_reads = int(subprocess.check_output(bash_command, shell=True, text=True))
  print("total", n_reads)
  dct = get_isolate_meta()["species"]
  # print(dct.keys())
  all_dict = dict()
  total = np.zeros((len(methods),),dtype=int)
  fp = np.zeros((len(methods),),dtype=int)
  snp_tp = np.zeros((len(methods),),dtype=int)
  snp_fp = np.zeros((len(methods),),dtype=int)
  snp_fn = np.zeros((len(methods),),dtype=int)
  for i in range(len(methods)):
    method = methods[i]
    bam_file = f"{result_dir}/{method}.filt.sorted.bam"
    vcf_file = f"{result_dir}/{method}.filt.vcf"
    bam = pysam.AlignmentFile(bam_file, "rb")
    species_count = dict()
    for read in bam.fetch():
      if read.is_supplementary or read.is_secondary:
      # print(".",end="",flush=True)
        continue
      else:
        total[i] += 1
        if not (read.query_name.startswith(identifier)):
          fp[i] += 1
          id = read.query_name.split(".")[0]
          if read.query_name[0] != 'S': species = read.query_name[:4]
          else: species = dct[id]
          if species not in species_count:
            species_count[species] = 1
          else: species_count[species] += 1
    if snp != None:
      # print("SNPs on")
      snp_fn[i],snp_fp[i],snp_tp[i] = read_vcf(vcf_file, snp)
    all_dict[method] =species_count
  with open(f'{result_dir}/FP_reads.json', 'w') as f:
    json.dump(all_dict, f)
  f.close()
  # print(total)
  fn = fp - total + n_reads * 2
  df = pd.DataFrame({"read_FN":fn, "read_FP":fp, "read_TP":total-fp})
  df['name'] = methods
  df["read_sensitivity"] = df["read_TP"] / (df['read_TP']+df['read_FN']) # tp / (tp + fn)
  df["read_precision"] = df['read_TP'] / (df['read_TP'] + df['read_FP']) # tp / (tp + fp)
  # print(df)
  snp_df = pd.DataFrame({"SNP_FN":snp_fn, "SNP_FP":snp_fp, "SNP_TP":snp_tp})
  snp_df['name'] = methods
  snp_df["SNP_sensitivity"] = snp_df["SNP_TP"] / (snp_df['SNP_TP']+snp_df['SNP_FN']) # tp / (tp + fn)
  snp_df["SNP_precision"] = snp_df['SNP_TP'] / (snp_df['SNP_TP'] + snp_df['SNP_FP']) 
  ans = pd.merge(df, snp_df, on="name")
  ans.index = ans['name']
  return ans

def get_contigs(species_name):
  '''
  @param: species_name is the full name of the species starting with contig(?) name
  @logic: first get the reference file name from metadata
  @returns: a list containing all the contigs from the reference genome
  '''
  metadata = "gtdb/gtdb_metadata.csv"
  meta_df = pd.read_csv(metadata)
  ref_filename = meta_df.loc[meta_df['name'] == f">{species_name}", 'genome_filename']
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

# plt.cla() clears the current axism, should use this
def handle_abund_species(sample,identifier="read"):
  abund_species_file = f"output/{sample}/abundant_species.csv"
  abund_df = pd.read_csv(abund_species_file)
  close_species_file = f"output/{sample}/close_species.csv"
  close_df = pd.read_csv(close_species_file)
  df = pd.merge(abund_df, close_df, on="name")
  df["name"] = df["name"].str.replace(">","")
  df = df[['name', "average_abund","median_abund","ani"]]
  selected_species = f"output/{sample}/abundant_species.fasta.txt"
  with open(selected_species) as f:
    text = f.read()
    lines = text.split('\n')
  df = df[df["name"].isin(lines)]
  df.index = df['name']
  for name in df['name']:
    specific_cont = get_contigs(name)
    bam_file = f"{result_dir}/abundant_species.filt.sorted.bam"
    # vcf_file = f"{result_dir}/{method}.filt.vcf"
    bam = pysam.AlignmentFile(bam_file, "rb")
    # species_count = dict()
    df.at[name, 'reads_gotten'] = 0
    for cont in specific_cont:
      for read in bam.fetch(cont):
        if read.query_name.startswith(identifier):
          df.at[name, 'reads_gotten'] += 1
  df.to_csv(f"output/{sample}/selected_species.csv")
  f.close()


if __name__ == "__main__":
  sample = snakemake.params.samp
  methods = snakemake.params.methods
  ref = snakemake.params.ref
  snp_lst = snakemake.params.snp_lst
  result_dir = f"output/{sample}/"
  fastq = f"data/{sample[:2]}1_reads_R1.fq.gz"
  if sample[0] == "1": identifier = "read"
  else: identifier = "ferg"

  # handle_abund_species(sample,identifier)

  df = check_read_n_snp(result_dir, methods, fastq, snp_lst,identifier)
  df = df.round(2)
  # plotting figures
  width = 0.25  # the width of the bars
  multiplier = 0
  x = np.arange(len(methods))
  read_df = df[["read_sensitivity","read_precision"]]
  read_df = read_df.stack().to_frame().reset_index()
  read_df = read_df.set_axis(['name', 'type', 'value'], axis=1)
  ax = sns.barplot(x = "name", y = "value", hue="type", data = read_df, width=0.5)
  for container in ax.containers:
    ax.bar_label(container)
  ax.legend(loc='lower right')
  ax.set_ylim(0, 1)
  plt.savefig(f"{result_dir}/read_mapping_perf.jpg")
  plt.cla()
  fig, ax = plt.subplots(constrained_layout=True)
  for what in ['SNP_sensitivity','SNP_precision']:
      attribute=what
      measurement=df[what]
      offset = width * multiplier
      rects = ax.bar(x + offset, measurement, width, label=attribute)
      ax.bar_label(rects, padding=3)
      multiplier += 1
  ax.set_xticks(x + width, methods)
  ax.legend(loc='lower right')
  ax.set_ylim(0, 1)
  plt.savefig(f"{result_dir}/SNP_perf.jpg")
  df.to_csv(f"{result_dir}/performance.csv")

