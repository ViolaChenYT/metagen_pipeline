import gzip, os
from collections import defaultdict

import pandas as pd
import pysam

ref = snakemake.input['ref']
bamfile = snakemake.input['bam']
read = snakemake.input['r1']
refname = snakemake.params['refname']
clist = snakemake.input['clist']
outfile = snakemake.output[0]

def read_fasta(fname):
	d = defaultdict(list)
	current = None
	for line in gzip.open(fname, 'rt'):
		line = line.strip()
		if line.startswith('>'):
			current = line[1:]
		elif current is None:
			raise Exception('No header found at the beggining of the file')
		else:
			d[current].append(line)

	return {
		contig: ''.join(seqs) for contig, seqs in d.items()	
	}

refcontigs = {contig.split()[0] for contig in read_fasta(ref).keys()}

# ----------------------


bamfile = pysam.AlignmentFile(bamfile, "rb")

counts = defaultdict(lambda: defaultdict(int))
read_refs_set = set()

for aln in bamfile.fetch(until_eof=True):
    name = aln.query_name
    rname = name.split('::')[1]
    contig = aln.reference_name
    counts[rname][contig] += 1

    if rname == refname:
    	read_refs_set.add(name)

bamfile.close()
read_refs_set = read_refs_set

df = pd.DataFrame([
	{'read_ori': name , 'bam_contig': contig, 'count': value}
	for name, sub in counts.items() for contig, value in sub.items()
	])

# ----------------------

def read_fasta(fname):
	d = defaultdict(list)
	current = None
	for line in gzip.open(fname, 'rt'):
		line = line.strip()
		if line.startswith('>'):
			current = line[1:]
		elif current is None:
			raise Exception('No header found at the beggining of the file')
		else:
			d[current].append(line)

	return {
		contig: ''.join(seqs) for contig, seqs in d.items()	
	}

df['read_ref'] = df['read_ori'] == refname
df['ref_contig'] = df['bam_contig'].isin(refcontigs)

def setkind(row):
	ref_contig, read_ref = row['ref_contig'], row['read_ref']
	if ref_contig and read_ref: return 'TP'
	elif ref_contig and not read_ref: return 'FP'
	elif read_ref and not ref_contig: return 'FN'
	else: return 'TN'

df['category'] = df.apply(setkind, axis=1)

# ----------------------

cdic = {}
with open(clist) as f:
	for line in f:
		line = line.strip()
		gid = os.path.basename(line.split(':')[0])
		cid = ':'.join(line.split(':')[1:]).split()[0][1:]
		cdic[cid] = gid

for contig in refcontigs:
	cdic[contig] = refname

cdic['None'] = 'PUnmapped'
df['contig_gid'] = df['bam_contig'].fillna('None').apply(lambda contig: cdic[contig])


# ----------------------

# The difference between Unmapped and PUnmapped
# Unmapped: Reads count completly lost after filtering
# PUnmapped: Reads defined as unmapped during mapping process (before filtering)

num_reads = set()
with open(read) as f:
	for line in f:
		if line.startswith('@'):
			line = line.strip()[1:-2]
			num_reads.add(line)

count_missing = len(num_reads - read_refs_set)

df = pd.concat([df, pd.DataFrame([{
	'read_ori': refname,
	'bam_contig': 'None',
	'count': count_missing,
	'read_ref': True,
	'ref_contig': False,
	'category': 'FN',
	'contig_gid': 'Unmapped'
	}])])

# ----------------------

df.to_csv(outfile, index=False, sep='\t')