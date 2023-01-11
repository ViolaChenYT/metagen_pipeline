************************
* ANIcalculator README *
************************
ANIcalculator (Average Nucleotide Identity Calculator)

Authors:
  Neha Varghese
  Amrita Pati
  Kostas Mavrommatis

Contact/questions:
  njvarghese@lbl.gov

Copyright (c) 2012-2014 : 
  DOE Joint Genomic Institute, Walnut Creek ,CA
  Lawrence Berkeley National Lab , Berkeley, CA

License: 
  LAWRENCE BERKELEY NATIONAL LAB END USER LICENSE AGREEMENT FOR NON-COMMERCIAL RESEARCH USE

Requirements:
  PERL version 5 or above

Compatibility:
  Linux

About:
  This tool will calculate the bidirectional average nucleotide identity (gANI) and Alignment Fraction (AF) between two genomes.
  Required input is the full path to the fna file (nucleotide sequence of genes in fasta format) of each query genome.
  Either the rRNA and tRNA genes can be excluded, or provided in a list with the -ignoreList option. 
  This is necessary as the presence of tRNA and/or rRNA genes in the fna will artificially inflate the ANI.

To run:
 ./ANIcalculator -help

Note:
  If you decide to move the 'ANIcalculator' binary to another location, 
  please ensure that the 'nsimscan' binary is placed in the same location as well.

