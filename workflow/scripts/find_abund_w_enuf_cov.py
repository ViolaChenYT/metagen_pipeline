import pysam
import gzip
import pandas as pd
import numpy as np


def test():
    inputfile = "/mnt/volume1/pipeline/output/spiked_kpneu/abundant_species.csv"
    data = pd.read_csv(inputfile)
    data.index = (data["name"])
    data = data[["average_abund", "f_unique_weighted"]]
    data = data.sort_values("f_unique_weighted", ascending=False)
    abund = data[data["average_abund"] > 5]
    abund = abund[abund["f_unique_weighted"] > 0.01]

    print("\n".join(abund.index.values))


if __name__ == "__main__":
    test()
