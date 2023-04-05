import pandas as pd
import numpy as np
import sys

species = ["combined", "Escherichia_coli_iai39", "metagenome"]
filters = ["nofilt", "loose_filt", "filt"]  #
iso_filters = ["filt", "nofilt"]

header = ["FN", "FP", "TP", "Sensitivity", "Precision"]


def analyze_tsv(filename):
    df = pd.read_csv(filename, sep="\t", header=0)
    tp_df = df[df["kind"] == "TP"]
    fn_df = df[df["kind"] == "FN"]
    fp_df = df[df["kind"] == "FP"]
    fn, fp, tp = len(fn_df), len(fp_df), len(tp_df)
    # print(tp["DP_bam"].mean())
    sensitivity = tp / (tp+fn)
    precision = 1 - fp/(fp+tp)
    return (f"{fn}\t{fp}\t{tp}\t{sensitivity}\t{precision}")


def analyze_instrain(filename):
    df = pd.read_csv(filename, sep="\t", header=0)


if __name__ == "__main__":
    if len(snakemake.params) < 1:
        raise Exception("please enter sample name, eg. MBH082")
    sample = snakemake.params[0]
    outfile = snakemake.output[0]
    with open(outfile, "w") as out:
        out.write(sample + "\t" + '\t'.join(header)+"\n")
        for sp in species:
            for filtr in filters:
                for isofilt in iso_filters:
                    filename = f"output/{sample}/compare_{sp}.{filtr}.iso_{isofilt}.lofreq.tsv"
                    out.write(f"{sp}_{filtr}_iso{isofilt}"+"\t")
                    # print("yeeet")
                    stats = analyze_tsv(filename)
                    out.write(stats+"\n")
    out.close()
    for sp in species:
        filename = f"output/{sample}/{sp}_profile_IS/output/{species}_profile_IS_SNVs.tsv"
        print(f"instrain\t{sp}")
        analyze_tsv(filename)
