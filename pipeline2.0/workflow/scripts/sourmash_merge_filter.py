import pandas as pd
import numpy as np

def filter_df(df, rani=0.95, utq=0.01, fmo=0.1, qani=0.7):
    maxani = df[df['ani'] >= rani]['ani'].max()
    if not np.isnan(maxani): df = df[df['ani'] < maxani]
        
    mask1 = df['f_unique_to_query'] > utq
    mask2 = df['f_match_orig'] >= fmo
    mask3 = df['ani'] >= qani
        
    return df[mask1 & mask2 & mask3]

gather = snakemake.input['gather'][0]
search = snakemake.input['search'][0]

gather = pd.read_csv(gather)
search = pd.read_csv(search)

gather = gather.merge(search, on='md5', how='left', suffixes=('_gather', '_search'))   
gather = filter_df(gather)

gather['gtdb_gid'] = gather['name_gather'].apply(lambda name: name.split(' ')[0])

outfile = snakemake.output[0]
gather.to_csv(outfile, index=False)