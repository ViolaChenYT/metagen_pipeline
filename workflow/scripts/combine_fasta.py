import pandas as pd
import numpy as np
import sys
import subprocess
from difflib import SequenceMatcher # for testing and checking

def filter(df, ani_t=0.95,utq_t=0.01,fmo_t=0.1, qani_t=0.7):
    mask1 = df['f_unique_to_query'] > utq_t
    mask2 = df['f_match_orig'] >= fmo_t
    
    return df[mask1 & mask2]

meta = "gtdb/gtdb_metadata.csv"

# input=$1 # takes a csv?
# std_ref=$3 # standard reference file including protocol. quite cool
# # sample=$2
# output=$2

# keyword="$4"
# close=$5
# threshold_hi=5 #$6
# threshold_lo=0.5
gather = sys.argv[1]
search = sys.argv[5]
outfile = sys.argv[2]
keyword = sys.argv[4]

gather = pd.read_csv(gather)
search = pd.read_csv(search)

gather = gather.merge(search, on='md5', how='left', suffixes=('_gather', '_search'))   
# print(gather.columns)
gather = filter(gather)
if 'name' in gather.columns:
    gather['name'] = gather['name'][:-1]
else:
    gather['name'] = gather['name_gather']

gather['gtdb_gid'] = gather['name_gather'].apply(lambda name: name.split(' ')[0])

def get_filename(name):
    pass

gather.to_csv(outfile, index=False)