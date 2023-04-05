from collections import defaultdict

import pysam
import pandas as pd
import numpy as np

def parse_vcf(fname):
	vf = pysam.VariantFile(fname)
	for record in vf.fetch():

		info = dict(record.info)
		info.pop('DP4')

		record = {
			'contig': record.contig,
			'position': record.pos,
			'ref': record.ref,
			'alts': '/'.join(record.alts),
			** info
		}

		yield record

def vcf_to_df(fname):
	return pd.DataFrame(
		parse_vcf(fname)
		)

def read_fasta(fname):
	d = defaultdict(list)
	current = None

	handler = gzip.open(fname, 'rt') if fname.endswith('gz') else open(fname)
	with handler as f:
		for line in f:
			line = line.strip()
			if line.startswith('>'):
				current = line[1:].split()[0]
			elif current is None:
				raise Exception('No header found at the beggining of the file')
			else:
				d[current].append(line)

	return {
		contig: ''.join(seqs) for contig, seqs in d.items()	
	}

def add_target_ref(row, fdata):
	contig, position = row['contig'], row['position']
	if not contig in fdata: return 'CNF'
	return fdata[contig][position - 1]

def set_vcalling_cat(row, cname):
    # TP: Variant correctly called
    # FN: Variant exist but not called
    # FP: Variant called but does not exist
    # FP1: The variant is not called at all
    # FP2: The variant is different from what's been called
    # cname if which row we consider for the true reference base    
    # FP1 only happen for cname == alts_ref
    # for target_ref, since we do not look at variant calling, only FP2 might occur
    # in this version I do not split multi alternatives and consider all alts as a block
    
  ref = row[cname]
  var = row['alts_decoys2ref']

  if isinstance(var, float) and np.isnan(var):
    return 'FN'
  if ref != var:
    if isinstance(ref, float) and np.isnan(ref):
      return 'FP1'
    else:
      return 'FP2'
  return 'TP'

# ------------------------------------------------------------------------------------------
# We compare variants calling between our result
# and what we have with our target reads wo the decoy (the simple isolate mapping)

vcf_wi_decoys = snakemake.input['vcf_wi_decoys']
vcf_wo_decoys = snakemake.input['vcf_wo_decoys']

vcf_wi_decoys = vcf_to_df(vcf_wi_decoys)
vcf_wo_decoys = vcf_to_df(vcf_wo_decoys)

df = vcf_wo_decoys.merge(vcf_wi_decoys, on=['contig', 'position'], suffixes=('_target2ref', '_decoys2ref'), how='outer')

# in case the target is a mutated (or not) version the reference sequence
# we can also assess the ground truth by checking the base at each position
# it is also possible that the target != reference (another strain)
# results obtained this way assess the true ground truth but might be impacted
# by other factor such as variant caller performance
target_ref = snakemake.input['target_ref']
fdata = read_fasta(target_ref)

df['ref_target'] = df.apply(lambda row: add_target_ref(row, fdata), axis=1)

df['cat_target2ref'] = df.apply(lambda row: set_vcalling_cat(row, 'alts_target2ref'), axis=1)
df['cat_rtarget'] = df.apply(lambda row: set_vcalling_cat(row, 'ref_target'), axis=1)

# because in our VCF file, we have data for the mapping reference containing decoy
# we need to filter out the contigs which belong to the decoys
mapping_ref = snakemake.input['mapping_ref']
fdata = read_fasta(mapping_ref)
df['ismapref'] = df['contig'].isin(set(fdata))

outfile = snakemake.output[0]
df.to_csv(outfile, sep='\t', index=False)