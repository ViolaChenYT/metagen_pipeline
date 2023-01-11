from collections import defaultdict
import re
import string
import pysam
import numpy as np
import pandas as pd
from argparse import ArgumentParser


def get_refs(ref_file):
    refs = dict()
    with open(ref_file, "r") as f:
        lines = f.readlines()
        current_line = ""
        contig = ""
        for i in range(len(lines)):

            line = lines[i]
            if line[0] == ">":
                contig = (line.split(" ")[0])[1:]
                if i > 0:
                    refs[contig] = current_line
                    current_line = ""
            else:
                current_line += line
        if current_line != "":
            refs[contig] = current_line
    f.close()
    return refs


def find_contig(contig, ref):
    for seq in ref:
        print(ref[seq].find(contig))
        return (ref[seq].find(contig))

def parse_vcf(fname):
    yield from pysam.VariantFile(fname).fetch()
    
def extract_vcf_stat(fname):
    ''' extract important fields from a VCF file:
        - contig name
        - position
        - reference and alternate alleles
        - strand bias, coverage, allele frequency
    '''
    return pd.DataFrame((
        {
            'contig': record.contig,
            'position': int(record.pos),
            'ref': record.alleles[0],
            'alt': alt,
            'SB': record.info["SB"],
            'DP': record.info["DP"],
            'AF': str(record.info["AF"]).split(',')[idx],
        }
        for record in parse_vcf(fname)
        for idx, alt in enumerate(record.alts)
    )
    )

def extract_lst(lst):
    ''' takes in the true SNP locations .lst file
        extracts contig name and original position of SNP
    '''
    lst = pd.read_csv(lst, sep="\t", header=2, dtype=str)
    # df = lst[lst["orig_position"] == lst["new_position"]]
    lst["og_position"] = lst["orig_position"].astype(int)
    lst["new_position"] = lst["new_position"].astype(int)
    lst["contig"] = lst['# seq_id'].str.slice(0, 11)
    return lst

def main(args):
    mummer = args.mummer
    mummer_dict = dict()
    with open(mummer, "r") as f:
      lines = f.readlines()
      contig = ""
      # build a dictionary using the mummer output
      for line in lines:
        if line[0] == ">":
          contig = line[1:].strip()
          mummer_dict[contig] = set()
        else:
            line = line.strip()
            og_pos, contig_pos, length = line.split()
            mummer_dict[contig].add((int(og_pos), int(contig_pos), int(length)))
          # somehow save the next few lines into dict until the next contig
    f.close()
    # print(mummer_dict)
    def get_orig_pos(contig, pos):
        for (ref, start, length) in mummer_dict[contig]:
            if pos >= start and pos < start + length:
                return pos - start + ref
    # read in the list of true SNPS
    true_snps = extract_lst(args.reference) # the offset should still be regarding the original position
    snps = set(true_snps["og_position"])
    # read in the list of called SNPs from VCF
    variant_file = args.query
    vcf = extract_vcf_stat(variant_file)
    vcf["remapped_pos"] = np.vectorize(get_orig_pos)(vcf["contig"], vcf["position"])
    vcf["TP"] = vcf["remapped_pos"].isin(snps)
    vcf["isnull"] = vcf["remapped_pos"].isnull()
    print(vcf["TP"].value_counts())
    vcf.to_csv(args.output, sep='\t')
    
if __name__ == "__main__":
    parser = ArgumentParser(
        description="from mummer output check accuracy of SNP calling with assembled genome")
    parser.add_argument(
        "-m", "--mummer", help="mummer output file with reference sequence as ref and assembled genome (de-repped) as query", required=True)
    parser.add_argument(
        "-r", "--reference", help="SNP list from the correct reference file to map to, *.lst", required=True)
    parser.add_argument(
        "-q", "--query", help="VCF file from assembled genome after binning and dereplication", required=True)
    parser.add_argument(
        "-o", "--output", help="output tsv file containing all SNP locations and attributes", required=True)
    args = parser.parse_args()
    main(args)


_DEBUG = 0


# "Fast Text Searching with Errors" by Sun Wu and Udi Manber
# TR 91-11, Dept of Computer Science, University of Arizona, June 1991.
# http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.20.8854

def WM_approx_Ham1_search(pattern, text):
    """Generate (Hamming_dist, start_offset)
    for matches with distance 0 or 1"""
    m = len(pattern)
    S_table = defaultdict(int)
    for i, c in enumerate(pattern):
        S_table[c] |= 1 << i
    R0 = 0
    R1 = 0
    mask = 1 << (m - 1)
    for j, c in enumerate(text):
        S = S_table[c]
        shR0 = (R0 << 1) | 1
        R0 = shR0 & S
        R1 = ((R1 << 1) | 1) & S | shR0
        if _DEBUG:
            print("j= %2d msk=%s S=%s R0=%s R1=%s"
                  % tuple([j] + map(bitstr, [mask, S, R0, R1])))
        if R0 & mask:  # exact match
            yield 0, j - m + 1
        elif R1 & mask:  # match with one substitution
            yield 1, j - m + 1


if _DEBUG:

    def bitstr(num, mlen=8):
        wstr = ""
        for i in xrange(mlen):
            if num & 1:
                wstr = "1" + wstr
            else:
                wstr = "0" + wstr
            num >>= 1
        return wstr


def Ham_dist(s1, s2):
    """Calculate Hamming distance between 2 sequences."""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def long_check(pattern, text):
    """Naively and understandably generate (Hamming_dist, start_offset)
    for matches with distance 0 or 1"""
    m = len(pattern)
    for i in xrange(len(text) - m + 1):
        d = Ham_dist(pattern, text[i:i+m])
        if d < 2:
            yield d, i


def Paul_McGuire_regex(pattern, text):
    searchSeqREStr = (
        '('
        + pattern
        + ')|('
        + ')|('.join(
            pattern[:i]
            + "[ACTGN]".replace(c, '')
            + pattern[i+1:]
            for i, c in enumerate(pattern)
        )
        + ')'
    )
    searchSeqRE = re.compile(searchSeqREStr)
    for match in searchSeqRE.finditer(text):
        locn = match.start()
        dist = int(bool(match.lastindex - 1))
        yield dist, locn


# if __name__ == "__main__":

#     genome1 = "TTTACGTAAACTAAACTGTAA"
#     #         01234567890123456789012345
#     #                   1         2

#     tests = [
#         (genome1, "ACGT ATGT ACTA ATCG TTTT ATTA TTTA"),
#         ("T" * 10, "TTTT"),
#         ("ACGTCGTAAAA", "TCGT"), # partial match can shadow an exact match
#         ]

#     nfailed = 0
#     for genome, patterns in tests:
#         print "genome:", genome
#         for pattern in patterns.split():
#             print pattern
#             a1 = list(WM_approx_Ham1_search(pattern, genome))
#             a2 = list(long_check(pattern, genome))
#             a3 = list(Paul_McGuire_regex(pattern, genome))
#             print a1
#             print a2
#             print a3
#             print a1 == a2, a2 == a3
#             nfailed += (a1 != a2 or a2 != a3)
#     print "***", nfailed


