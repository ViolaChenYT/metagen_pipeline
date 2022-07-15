from math import nan
from turtle import position
import pysam
from pysam import VariantFile
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from os.path import exists
import sys
import seaborn as sns

def read_vcf(filename):
    ''' open a vcf file and return the reads
    does not return header information
    '''
    vcf = pysam.VariantFile(filename)
    records = vcf.fetch()
    return records

def compare_snps(isolate, instrain, inhouse):
  '''compare snps from isolate-reference vs
  -  instrain output SNVs
  -  inhouse called variants
  @param: isolate: path to isolate to reference vcf file
  @param: instrain: path of instrain outputed SNVs tsv file
  @param: inhouse: path to lofreq output vcf file
  '''
  all_snps = dict()
  # isolate_snvs = read_vcf(isolate)
  # for iso_snp in isolate_snvs:
  #   all_snps[iso_snp.pos] = f"{iso_snp.ref}->{iso_snp.alts[0]}"
  # n = len(all_snps)
  hit, miss = 0,0
  inhouse_snvs = read_vcf(inhouse)
  lf_snps = set()
  for snp in inhouse_snvs:
    lf_snps.add(snp.pos)
    if snp.pos in all_snps:
      hit += 1
    else: 
      miss += 1
  isolate_snvs = read_vcf(isolate)
  for snp in isolate_snvs:
    if not snp.pos in lf_snps:
      print(snp,end="")


def analyze(filename, label=""):
  vcf_file = read_vcf(filename)
  covs = []
  sbs = []
  afs = []
  for var in vcf_file:
    sbs.append(var.info["SB"])
    covs.append(var.info["DP"])
    afs.append(var.info["AF"])
  sbs = np.array(sbs)
  covs = np.array(covs)
  afs = np.array(afs)
  print(label)
  print("Strand bias:", np.mean(sbs), np.std(sbs), end="\t")
  print("Coverage:", np.mean(covs), np.std(covs))
  print("AF:", np.mean(afs), np.std(afs))
  return (covs, afs)

def check_bam(bam, varlist):
  var = read_vcf(varlist)
  varlist = set()
  for v in var:
    varlist.add(v.pos)
  bases = ['A','T','C','G']
  alleles = dict((pos,{'A':0,'T':0,'C':0,'G':0}) for pos in varlist)
  bamfile = pysam.AlignmentFile(bam,"rb")
  for col in bamfile.pileup():
    # print(col)
    if col.pos not in varlist:
      continue
    for read in col.pileups:
      base = read.alignment.query_sequence[read.query_position_or_next]
      alleles[col.pos][base] += 1
  cov = np.array([sum(d.values()) for d in alleles.values()])
  print(np.mean(cov),np.std(cov))
  def get_AF(d):
    total = sum(d.values())
    if total == 0:
      return nan
    allele = 0
    max_cnt = 0
    for b in bases:
      if d[b] > 0:
        allele += 1
      if d[b] > max_cnt:
        max_cnt = d[b]
    # if allele > 2:
    #   print(f"multiallelic")
    return max_cnt / total
  afs = np.array([get_AF(d) for d in alleles.values()])
  print(np.nanmean(afs), np.nanstd(afs))
  return cov, afs[~np.isnan(afs)]

def gen_data():
  fn_sample_cov, fn_sample_af = check_bam(f"/mnt/volume1/cpe/{sample}/ec_filt.bam",false_neg)
  print(len(fn_sample_af))
  tp_sample_cov, tp_sample_af = check_bam(f"/mnt/volume1/cpe/{sample}/ec_filt.bam",true_pos)
  print(len(tp_sample_af))

  fn_iso_cov, fn_iso_af = check_bam(f"/mnt/volume1/isolate/{isolate}/ec_filt.bam",false_neg)
  print(len(fn_iso_af))
  tp_iso_cov, tp_iso_af = check_bam(f"/mnt/volume1/isolate/{isolate}/ec_filt.bam",true_pos)
  print(len(tp_iso_af))

  fn_cov = [fn_sample_cov,fn_iso_cov]
  fn_af = [fn_sample_af,fn_iso_af]
  tp_cov = [tp_sample_cov,tp_iso_cov]
  tp_af = [tp_sample_af,tp_iso_af]
  np.save("fn_cov",fn_cov)
  np.save("fn_af",fn_af)
  np.save("tp_cov",tp_cov)
  np.save("tp_af",tp_af)

def plot_snp_cov(npy,vcf, figname,binsize=10000):
  [sample,iso] = np.load(npy)
  
  bins = np.empty((5200000//binsize,))
  idx = 0
  snps = read_vcf(vcf)
  print(len(sample))
  for var in snps:
    if idx >= len(sample):
      break
    bins[var.pos//binsize] += 1 #sample[idx]
    idx += 1
  plt.scatter(binsize * np.arange(len(bins)), bins)
  plt.ylabel("average no. of SNPs per bin")
  plt.savefig(figname)

if __name__ == "__main__":
  isolate = "WEB1012"
  iso = f"/mnt/volume1/pipeline/output/{isolate}/Escherichia_coli_iai39.filt.lofreq.vcf"
  sample = "MBH082"
  instrain = f"/mnt/volume1/instrain_out/{sample}/output/{sample}_SNVs.tsv"
  inhouse = f"/mnt/volume1/pipeline/output/{sample}/Escherichia_coli_iai39.filt.lofreq.vcf"

  plot_snp_cov("../scripts/fn_cov.npy","false_neg.vcf","false_neg_cov.png")
  # compare_snps(iso, instrain,inhouse)

  # analyze(iso)
  # analyze_false()

  # cov, afs = check_bam(f"/mnt/volume1/cpe/{sample}/ec_filt.bam",false_pos)
  # gen_data()
  # kind = "cov"
  # sample = "tp"
  # def expand(s):
  #   if s == "fn":
  #     return "False negatives"
  #   else: return "True positives"
  # cov = np.load(f"{sample}_{kind}.npy",allow_pickle=True)
  # n = len(cov[0])
  # print(n)
  # samp,iso = cov[0].reshape((n,1)), cov[1].reshape((n,1))
  # fn_cov = pd.DataFrame(np.hstack((np.log10(samp,where=samp>0),np.log10(iso,where=iso>0))))
  # fn_cov.columns = ['Sample','Isolate']
  # fn_cov.plot(kind="box")
  # # cov[0],cov[1] = np.log(cov[0]), np.log(cov[1])
  # # plt.boxplot(cov, manage_ticks=True)
  # # plt.xticks([1,2],['Sample','Isolate'])
  # plt.ylabel("Log Coverage")
  # plt.title(expand(sample))
  # plt.savefig(f"{sample}_{kind}.png")

  # iso = "/mnt/volume1/result/WEB1015/filt.bam.vcf"
  # sample = "MBH086"
  # instrain = f"/mnt/volume1/instrain_out/{sample}/output/{sample}_SNVs.tsv"
  # inhouse = f"/mnt/volume1/result/{sample}/filt.bam.vcf"
  # compare_snps(iso, instrain,inhouse)

  # iso = "/mnt/volume1/result/WEB1352/filt.bam.vcf"
  # sample = "pos_control"
  # instrain = f"/mnt/volume1/instrain_out/{sample}/output/{sample}_SNVs.tsv"
  # inhouse = f"/mnt/volume1/result/{sample}/filt.bam.vcf"
  # compare_snps(iso, instrain,inhouse)