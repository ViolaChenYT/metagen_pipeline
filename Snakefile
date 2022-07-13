from pathlib import Path

rule all:
    input:
        expand("mapping/{sample}.filt.{mapper}.annotated.vcf", sample=config, mapper=('lofreq',))

rule link_ref:
    input:
        lambda wc: Path(config[wc.sample]['ref']).resolve()
    output:
        'refs/{sample}.fasta'
    shell:
        'ln -s {input} {output}'

rule bowtie_build:
    input:
        "refs/{sample}.fasta"
    output:
        "refs/{sample}.fasta.1.bt2",
        "refs/{sample}.fasta.2.bt2",
        "refs/{sample}.fasta.3.bt2",
        "refs/{sample}.fasta.4.bt2",
        "refs/{sample}.fasta.rev.1.bt2",
        "refs/{sample}.fasta.rev.2.bt2"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "bowtie2-build {input} {input}"

rule alignment:
    input:
        r1 = lambda wc: config[wc.sample]['r1'],
        r2 = lambda wc: config[wc.sample]['r2'],
        ref = rules.link_ref.output,
        idx = rules.bowtie_build.output
    output:
        'mapping/{sample}.bam'
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "bowtie2 -p 6 -x {input.ref} -1 {input.r1} -2 {input.r2} | "
        "samtools view -@ 16 -b -u - | samtools sort -o {output}"

rule samtools_index:
    input:
        rules.alignment.output
    output:
        'mapping/{sample}.bam.bai'
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        'samtools index {input}'

rule bamfilter:
    input:
        rules.alignment.output,
        rules.samtools_index.output
    output:
        "mapping/{sample}.filt.bam"
    conda:
        'workflow/envs/mapping.yaml'
    script:
        'workflow/scripts/bamfilter.py'

rule lofreq:
    input:
        ref = rules.link_ref.output,
        filt_bam = rules.bamfilter.output
    output:
        "mapping/{sample}.filt.lofreq.vcf"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "lofreq call -f {input.ref} -o {output} --verbose {input.filt_bam}"

rule pileup:
    input:
        ref = rules.link_ref.output,
        bam = rules.alignment.output
    output:
        gz = "mapping/{sample}.pileup.gz"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "samtools mpileup -E -f {input.ref} {input.bam} | gzip > {output.gz}"

rule coverage_txt:
    input:
        rules.bamfilter.output
    output:
        "mapping/{sample}.filt.coverage.txt"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "genomeCoverageBed -d -ibam {input} > {output}"

rule coverage_avg:
    input:
        bam = rules.alignment.output,
        bai = rules.samtools_index.output
    output:
        "mapping/{sample}.filt.coverage.summary.txt"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "samtools depth -aa {input.bam} | awk '{{sum+=$3}} END {{print sum/NR}}' > {output}"

rule annotate_gene:
    input:
        "mapping/{sample}.filt.{mapper}.vcf"
    output:
        "mapping/{sample}.filt.{mapper}.annotated.vcf"
    params:
        snpeff_species = lambda wc: config[wc.sample]['snpeff_species']
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "java -Xmx8g -jar workflow/soft/snpEff/snpEff.jar {params.snpeff_species} {input} > {output}"