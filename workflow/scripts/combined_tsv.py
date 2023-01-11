import pandas as pd
import numpy as np
import sys

species = snakemake.params.species
filters = ["filt"]  #
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
        sample = sys.argv[1]
        outfile = sys.argv[2]
        # raise Exception("please enter sample name, eg. MBH082")
    else:
        sample = snakemake.params.sample
        outfile = snakemake.output[0]
    with open(outfile, "w") as out:
        out.write(sample + "\t" + '\t'.join(header)+"\n")
    
        for tsv in snakemake.input:
            name = tsv.split("/")[1] + ("").join(tsv.split("/")[2].split(".")[:2])[8:]
            out.write(f"{name}"+"\t")
            # print("yeeet")
            stats = analyze_tsv(tsv)
            out.write(stats+"\n")
            instrain = False
            if instrain:
                for sp in species:
                    for isofilt in iso_filters:
                        filename = f"output/{sample}/{sp}_instrain_iso_{isofilt}.lofreq.tsv"
                        out.write(f"instrain_{sp}_iso{isofilt}\t")
                        stats = analyze_tsv(filename)
                        out.write(stats+"\n")
    out.close()
