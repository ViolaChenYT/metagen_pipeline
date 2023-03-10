import utils

# ----------------------------------------------------------------
# Intermediate files for wo_decoy

rule link_ref:
    input:
        lambda wc: utils.get_compressed_ref(config, wc.sample)
    output:
        'mapping/{sample}/refs/wodecoys/ref.fa.gz'
    shell:
        'ln -s $(pwd)/{input} {output}'

rule bwa_index_wo_decoys:
    input:
        tmp = lambda wc: utils.get_compressed_ref(config, wc.sample),
        ref = rules.link_ref.output
    output:
        'mapping/{sample}/refs/wodecoys/ref.fa.gz.amb'
    conda:
        config['softparams']['conda']['bwasam']
    log:
        'mapping/{sample}/refs/wodecoys/ref.fa.gz.log'
    shell:
        'bwa index {input.ref} 2> {log}'

rule bwa_wo_decoys:
    input:
        r1 = lambda wc: config['samples'][wc.sample]['r1'],
        r2 = lambda wc: config['samples'][wc.sample]['r2'],
        tmp = lambda wc: utils.get_compressed_ref(config, wc.sample),
        ref = rules.link_ref.output,
        idx = rules.bwa_index_wo_decoys.output
    output:
        temp('mapping/{sample}/aln.wodecoys.bam')
    conda:
        config['softparams']['conda']['bwasam']
    threads:
        16
    log:
        'mapping/{sample}/aln.wodecoys.log'
    shell:
        "(bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -b | samtools sort -n --threads {threads} > {output}) 2> {log}"

# ----------------------------------------------------------------
# Mapping assessment

rule assess_mapping:
    input:
        ref = lambda wc: utils.get_compressed_ref(config, wc.sample),
        bam = 'mapping/{sample}/aln.{suffix}.bam',
        r1  = lambda wc: utils.get_target_r1(config, wc.sample),
        clist = 'sourmash/res_gtdb/samples/{sample}.merged.clist.txt'
    params:
        refname = lambda wc: config['samples'][wc.sample]['target']
    output:
        'assessement/{sample}/mapping_assessement.{suffix}.tsv'
    conda:
        config['softparams']['conda']['bwasam']
    script:
        'scripts/mapping_assessement.py'

# ----------------------------------------------------------------
# Intermediate files for VCF assessment

rule bwa_target_reads:
    input:
        r1 = lambda wc: utils.get_target_r1(config, wc.sample),
        r2 = lambda wc: utils.get_target_r2(config, wc.sample),
        tmp = lambda wc: utils.get_compressed_ref(config, wc.sample),
        ref = rules.link_ref.output,
        idx = rules.bwa_index_wo_decoys.output
    output:
        temp('mapping/{sample}/aln.wodecoys_target.bam')
    conda:
        config['softparams']['conda']['bwasam']
    threads:
        16
    log:
        'mapping/{sample}/aln.wodecoys_target.log'
    shell:
        "(bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -b | samtools sort -n --threads {threads} > {output}) 2> {log}"


# ----------------------------------------------------------------
# VCF assessment

rule assess_vcf:
    input:
        vcf_wi_decoys = 'mapping/{sample}/aln.{kind}.fixmate.sorted.filtered.vcf',
        vcf_wo_decoys = 'mapping/{sample}/aln.wodecoys_target.fixmate.sorted.filtered.vcf',
        target_ref    = lambda wc: utils.get_target_mutated_ref(wc, config),
        mapping_ref   = lambda wc: config['samples'][wc.sample]['ref']
    output:
        'assessement/{sample}/vcalling_assessement.{kind}.tsv'
    conda:
        config['softparams']['conda']['bwasam']
    script:
        'scripts/vcalling_assessement.py'