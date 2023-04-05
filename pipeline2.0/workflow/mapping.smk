import utils

'''
@ Summary of JS's doc:
- utils for compress/decompress
- alignment with bwa
- samtools fixmate and sort / index
- bam filter ** impt ** 
- samtools depth
- lofreq 
- lofreq faidx
'''

rule compress_ref:
    input:
        lambda wc: config['samples'][wc.sample]['ref']
    output:
        temp('mapping/{sample}/refs/decoys/ref.comp.fa.gz')
    shell:
        'cat {input} | gzip > {output}'

rule make_reference:
    input:
        decoys = 'sourmash/res_gtdb/samples/{sample}.merged.fa.gz',
        ref    = lambda wc: utils.get_compressed_ref(config, wc.sample)
    output:
        temp('mapping/{sample}/refs/decoys/ref.fa.gz')
    shell:
        'cat {input.decoys} {input.ref} > {output}'

rule bwa_index:
    input:
        rules.make_reference.output
    output:
        'mapping/{sample}/refs/decoys/ref.fa.gz.amb'
    conda:
        config['softparams']['conda']['bwasam']
    log:
        'mapping/{sample}/refs/decoys/ref.fa.gz.log'
    shell:
        'bwa index {input} 2> {log}'

rule bwa:
    input:
        r1 = lambda wc: config['samples'][wc.sample]['r1'],
        r2 = lambda wc: config['samples'][wc.sample]['r2'],
        ref = rules.make_reference.output,
        idx = rules.bwa_index.output
    output:
        temp('mapping/{sample}/aln.decoys.bam')
    conda:
        config['softparams']['conda']['bwasam']
    threads:
        16
    log:
        'mapping/{sample}/aln.log'
    shell:
        "(bwa mem -t {threads} {input.ref} {input.r1} {input.r2}"
        " | samtools view -b | samtools sort -n --threads {threads} > {output}) 2> {log}"

rule samtools_fixmate:
    input:
        'mapping/{sample}/aln.{bamprefix}.bam'
    output:
        temp('mapping/{sample}/aln.{bamprefix}.fixmate.bam')
    conda:
        config['softparams']['conda']['bwasam']
    shell:
        'samtools fixmate {input} {output}'

rule samtools_sort_coordinates:
    input:
        rules.samtools_fixmate.output
    output:
        temp('mapping/{sample}/aln.{bamprefix}.fixmate.sorted.bam')
    conda:
        config['softparams']['conda']['bwasam']
    shell:
        "samtools sort -o {output} {input}"

rule samtools_index:
    input:
        rules.samtools_sort_coordinates.output
    output:
        temp('mapping/{sample}/aln.{bamprefix}.fixmate.sorted.bam.bai')
    conda:
        config['softparams']['conda']['bwasam']
    shell:
        'samtools index {input}'

rule bam_filtering:
    input:
        bam = rules.samtools_sort_coordinates.output,
        bai = rules.samtools_index.output
    output:
        'mapping/{sample}/aln.{bamprefix}.fixmate.sorted.filtered.bam'
    params:
        tcov = 0.9,
        tsim = 0.05
    conda:
        config['softparams']['conda']['bwasam']
    log:
        'mapping/{sample}/bamfiltering.{bamprefix}.log'
    script:
        'scripts/bamfilter.py'

rule depth:
    input:
        rules.bam_filtering.output
    output:
        temp('mapping/{sample}/aln.{bamprefix}.depth.txt.gz')
    conda:
        config['softparams']['conda']['bwasam']
    shell:
        'samtools depth -aa {input} | gzip > {output}'

# ----------------------------------------------------------------------
# lofreq

rule gz2bgz:
    input: 
        'mapping/{sample}/refs/{refori}/ref.fa.gz'
    output:
        'mapping/{sample}/refs/{refori}/ref.fa.bgz'
    conda:
        config['softparams']['conda']['bwasam']
    shell:
        'gunzip -c {input} | bgzip  > {output}'

rule lofreq_faidx:
    input:
        rules.gz2bgz.output
    output:
        'mapping/{sample}/refs/{refori}/ref.fa.bgz.fai'
    conda:
        config['softparams']['conda']['lofreq']
    shell:
        'lofreq faidx {input}'

rule lofreq:
    input:
        ref = lambda wc: utils.get_ref_path_vcalling(wc),
        fai = lambda wc: utils.get_ref_path_vcalling(wc, '.fai'),
        bam = rules.bam_filtering.output
    output:
        rules.bam_filtering.output[0][:-4] + '.vcf'
    conda:
        config['softparams']['conda']['lofreq']
    log:
        rules.bam_filtering.output[0][:-4] + '.lofreq.log'
    shell:
        'lofreq call -f {input.ref} -o {output} --verbose {input.bam} 2> {log}'