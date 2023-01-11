#!/usr/bin/python

import pandas as pd
import numpy as np
import sys
import pysam

bam_file = sys.argv[1]


def check_bam():
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for read in samfile.fetch():
        print(read.query_name)
        print(".", end="", flush=True)
        # print(read.get_cigar_stats())
        # print(read.get_reference_name())


if __name__ == "__main__":
    check_bam()
