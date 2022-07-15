from pathlib import Path

rule all:
  input:
    expand("output/{sample}/{species}.{filt}.{mapper}.annotated.vcf", \
        sample=config, species=("Escherichia_coli_iai39",), mapper=('lofreq',), \
        filt=('filt',)),
    expand("output/{sample}/{species}.{filt}.coverage.txt.gz", \
        sample=config, species=("Escherichia_coli_iai39",), mapper=('lofreq',), \
        filt=('filt',))

# rule link_ref:
#   input:
#     lambda wc: Path(config[wc.sample]['ref']).resolve()
#   output:
#     'refs/{sample}_Escherichia_coli_iai39.fasta'
#   shell:
#     'ln -s {input} {output}'
# 
#   doesn't make much sense if every single sample has its own reference  

REFERENCE = "refs/Escherichia_coli_iai39.fasta"

rule bowtie_build:
    input:
        "refs/{species}.fasta"
    output:
        "refs/{species}.fasta.1.bt2",
        "refs/{species}.fasta.2.bt2",
        "refs/{species}.fasta.3.bt2",
        "refs/{species}.fasta.4.bt2",
        "refs/{species}.fasta.rev.1.bt2",
        "refs/{species}.fasta.rev.2.bt2"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "bowtie2-build {input} {input}"

rule alignment:
    input:
        r1 = lambda wc: config[wc.sample]['r1'],
        r2 = lambda wc: config[wc.sample]['r2'],
        ref = REFERENCE,
        idx = rules.bowtie_build.output
    output:
        'output/{sample}/{species}.bam'
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "bowtie2 -p 6 -x {input.ref} -1 {input.r1} -2 {input.r2} | "
        "samtools view -@ 16 -b -u - | samtools sort -o {output}"

rule alignment_mine:
    input:
        r1 = "{homedir}/{sample}/merge_{sample}_R1.fastq.gz" ,
        r2 = "{homedir}/{sample}/merge_{sample}_R2.fastq.gz" ,
        ref = REFERENCE,
        idx = rules.bowtie_build.output
    output:
        "{homedir}/{sample}/{species}.bam"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "bowtie2 -p 6 -x {input.ref} -1 {input.r1} -2 {input.r2} | "
        "samtools view -@ 16 -b -u - | samtools sort -o {output}"

rule samtools_index:
    input:
        rules.alignment.output
    output:
        'output/{sample}/{species}.bam.bai'
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        'samtools index {input}'

rule bamfilter:
    input:
        rules.alignment.output,
        rules.samtools_index.output
    output:
        "output/{sample}/{species}.filt.bam"
    conda:
        'workflow/envs/mapping.yaml'
    script:
        'workflow/scripts/bamfilter.py'

# rule lose_filter:
#     input:
#         rules.alignment.output,
#         rules.samtools_index.output
#     output:
#         "output/{sample}/{species}.{filt}.bam"
#     conda:
#         'workflow/envs/mapping.yaml'
#     script:
#         'workflow/scripts/bamfilter.py -l'

rule lofreq:
    input:
        ref = REFERENCE,
        filt_bam = rules.bamfilter.output
    output:
        "output/{sample}/{species}.{filt}.lofreq.vcf"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "lofreq call -f {input.ref} -o {output} --verbose {input.filt_bam}"

rule pileup:
    input:
        ref = REFERENCE,
        bam = rules.alignment.output
    output:
        gz = "output/{sample}/{sample}.pileup.gz"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "samtools mpileup -E -f {input.ref} {input.bam} | gzip > {output.gz}"

rule coverage_txt:
    input:
        rules.bamfilter.output
    output:
        "output/{sample}/{species}.{filt}.coverage.txt.gz"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "genomeCoverageBed -d -ibam {input} | gzip > {output}"

rule coverage_avg:
    input:
        bam = rules.alignment.output,
        bai = rules.samtools_index.output
    output:
        "output/{sample}/{species}.{filt}.coverage.summary.txt"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "samtools depth -aa {input.bam} | awk '{{sum+=$3}} END {{print sum/NR}}' > {output}"

rule annotate_gene:
    input:
        "output/{sample}/{species}.{filt}.{mapper}.vcf"
    output:
        "output/{sample}/{species}.{filt}.{mapper}.annotated.vcf"
    params:
        snpeff_species = lambda wc: config[wc.sample]['snpeff_species']
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "(snpeff {params.snpeff_species} {input} > {output}) && mv snpEff* output/{wildcards.sample}/"
