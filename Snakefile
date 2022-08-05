from pathlib import Path
import os

rule all:
  input:
    expand("output/{sample}/{species}.{filt}.bam", \
        sample=config, species=("metagenome"), \
        filt=("nofilt","filt","loose_filt")),
    # expand("output/{sample}/{species}.{filt}.{mapper}.vcf", \
    #     sample=config, species=("metagenome"), mapper=('lofreq',), \
    #     filt=("nofilt","filt","loose_filt")),
    # expand("output/{sample}/compare_{species}.{filt}.iso_{isofilt}.{mapper}.tsv", \
    #     sample=config, species=("Escherichia_coli_iai39","metagenome"), mapper=('lofreq',), \
    #     filt=( "nofilt","filt","loose_filt",), isofilt=("filt", "nofilt",)),
    # expand("output/{sample}/compare_result.tsv", sample=config)
    # "refs/metagenome.fasta"
    # expand("output/{sample}/{species}_profile_IS", \
        # sample=config, species=("combined","Escherichia_coli_iai39","metagenome",))
# rule link_ref:
#   input:
#     lambda wc: Path(config[wc.sample]['ref']).resolve()
#   output:
#     'refs/{sample}_Escherichia_coli_iai39.fasta'
#   shell:
#     'ln -s {input} {output}'
# 
#   doesn't make much sense if every single sample has its own reference  


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

rule alignment_isolate:
    input:
        r1 = "/mnt/volume1/isolate/{sample}/merge_{sample}_R1.fastq.gz" ,
        r2 = "/mnt/volume1/isolate/{sample}/merge_{sample}_R2.fastq.gz" ,
        ref = "refs/{species}.fasta",
        idx = rules.bowtie_build.output
    output:
        "output/{sample}/{species}.bam"
    conda:
        'workflow/envs/mapping.yaml'
    threads:
        8
    shell:
        "bowtie2 -a -p {threads} -x {input.ref} -1 {input.r1} -2 {input.r2} | "
        "samtools sort -o {output} -@ {threads}"

rule alignment:
    input:
        r1 = lambda wc: config[wc.sample]['r1'],
        r2 = lambda wc: config[wc.sample]['r2'],
        ref = "refs/{species}.fasta",
        idx = rules.bowtie_build.output
    output:
        'output/{sample}/{species}.bam'
    conda:
        'workflow/envs/mapping.yaml'
    threads:
        8
    shell:
        "bowtie2 -k 5 -p {threads} -x {input.ref} -1 {input.r1} -2 {input.r2} | "
        "samtools sort -o {output} -@ {threads}"



rule samtools_index:
    input:
        "output/{sample}/{species}.bam"
    output:
        'output/{sample}/{species}.bam.bai'
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        'samtools index {input}'

rule samtools_index_filtered:
    input:
        "output/{sample}/{species}.{filt}.bam"
    output:
        'output/{sample}/{species}.{filt}.bam.bai'
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        'samtools index {input}'

rule bamfilter:
    input:
        'output/{sample}/{species}.bam',
        'output/{sample}/{species}.bam.bai'
    output:
        "output/{sample}/{species}.filt.bam"
    conda:
        'workflow/envs/mapping.yaml'
    params:
        "-s"
    script:
        'workflow/scripts/bamfilter.py'

rule loose_filter:
    input:
        'output/{sample}/{species}.bam',
        rules.samtools_index.output
    output:
        "output/{sample}/{species}.loose_filt.bam"
    conda:
        'workflow/envs/mapping.yaml' 
    params:
        "-l"
    script:
        'workflow/scripts/bamfilter.py'

rule no_filter:
    input:
        'output/{sample}/{species}.bam',
        rules.samtools_index.output
    output:
        "output/{sample}/{species}.nofilt.bam"
    conda:
        'workflow/envs/mapping.yaml'
    params:
        "-n"
    script:
        'workflow/scripts/bamfilter.py'

rule lofreq:
    input:
        ref = "refs/{species}.fasta",
        filt_bam = "output/{sample}/{species}.{filt}.bam"
    output:
        "output/{sample}/{species}.{filt}.lofreq.vcf"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "lofreq call -f {input.ref} -o {output} --verbose {input.filt_bam}"

rule pileup:
    input:
        ref = "refs/{species}.fasta",
        bam = 'output/{sample}/{species}.bam'
    output:
        gz = "output/{sample}/{sample}.pileup.gz"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "samtools mpileup -E -f {input.ref} {input.bam} | gzip > {output.gz}"

rule coverage_txt:
    input:
        "output/{sample}/{species}.{filt}.bam"
    output:
        "output/{sample}/{species}.{filt}.coverage.txt.gz"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "genomeCoverageBed -d -ibam {input} | gzip > {output}"

rule coverage_avg:
    input:
        bam = 'output/{sample}/{species}.bam',
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
        "sed -i 's/NC_011750.1/Chromosome/g' {input} && "
        "(snpeff {params.snpeff_species} {input} > {output}) && "
        "mv ./snpEff* output/{wildcards.sample}/"

rule compare_snps:
    input:
        meta = "output/{sample}/{species}.{filt}.{mapper}.vcf",
        iso = lambda wc: os.path.join("output/", config[wc.sample]['isolate'] + "/Escherichia_coli_iai39." + wc.isofilt + "." + wc.mapper + ".vcf"),
        bam = "output/{sample}/{species}.{filt}.bam",
        bai = rules.samtools_index_filtered.output
    output:
        "output/{sample}/compare_{species}.{filt}.iso_{isofilt}.{mapper}.tsv"
    conda:
        'workflow/envs/mapping.yaml'
    params:
        "meta"
    script:
        "workflow/scripts/JS_script.py"

rule compare_instrain:
    input:
        meta = "output/{sample}/{species}_profile_IS/output/{species}_profile_IS_SNVs.tsv",
        iso = lambda wc: os.path.join("output/", config[wc.sample]['isolate'] + "/Escherichia_coli_iai39." + wc.isofilt + "." + wc.mapper + ".vcf")
    output:
        "output/{sample}/{species}_instrain_iso_{isofilt}.{mapper}.tsv"
    conda:
        'workflow/envs/mapping.yaml'
    params:
        "instrain"
    script:
        "workflow/scripts/JS_script.py"

rule prodigal:
    input:
        ref = "refs/{species}.fasta"
    output:
        "refs/{species}.genes.fna"
    shell:
        "prodigal -i {input.ref} -d {output}"

rule instrain:
    input:
        genes = rules.prodigal.output,
        ref = "refs/{species}.fasta",
        bam = 'output/{sample}/{species}.nofilt.bam'
    output:
        directory("output/{sample}/{species}_profile_IS")
    threads:
        8
    shell:
        "inStrain profile {input.bam} {input.ref} -g {input.genes} -o {output} -p 6"

rule compare_metrics:
    input:
    output:
        "output/{sample}/compare_result.tsv"
    params:
        "{sample}"
    conda:
        'workflow/envs/mapping.yaml'
    script:
        "combined_tsv.py"

# need to add rules for downloading database, unzip and create signatures

rule get_close_species:
    input:
        ref = "gtdb/Escherichia_coli_iai39.fasta.sig",
        db = "gtdb/databases.zip"
    output:
        "gtdb/Escherichia_coli_iai39.close_species.csv"
    conda:
        "workflow/envs/mapping.yaml"
    shell:
        "sourmash search {input.ref} {input.db} -o {output} --threshold=0"

rule build_metagenome:
    input:
        rules.get_close_species.output
    output:
        "refs/metagenome_tmp.fasta"
    shell:
        "./workflow/scripts/combine_fasta.sh {input} {output}"

rule dereplication:
    input:
        rules.build_metagenome.output
    output:
        "refs/metagenome.fasta"
    shell:
        "awk '/^>/{f=!d[$1];d[$1]=1}f' {input} > {output}"
