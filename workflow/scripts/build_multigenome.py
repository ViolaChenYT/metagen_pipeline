import pandas as pd
import gzip
import sys

if __name__ == "__main__":
  df = pd.read_csv(sys.argv[2],header=0)
  query = sys.argv[1]
  print(df)
  with open(snakemake.output[0], "w") as out:
      for filename in df["filename"]:
        with gzip.open(filename, "rt") as handle:
          for record in SeqIO.parse(handle, "fasta"):
            out.write(record)
  out.close()