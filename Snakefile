import os
from snakemake.io import glob_wildcards
method_list = ("target_species","abundant_species","database", "assembly")
filter_list = ("filt","unfilt")

# thresholds = (1, 2, 3, 4, 5, 6, 7, 8, 9) # threshold for coverage
# thresholds = (0.75, 0.77, 0.8, 0.82, 0.85) # threshold for ANI
# thresholds = (0.0, 0.2, 0.4,0.6,0.8, 0.9)

# og_ref = "refs/" + species.replace(" ", "_") + ".fasta"
wildcard_constraints:
   species = '\w+',
   filt = '\w+'

rule all:
  input:
    # expand("output/{sample}/abundant_species.{threshold}.filt.vcf",sample=config,threshold=thresholds),
    # expand("output/{sample}/assembly.{thresh}.filt.vcf",sample=config,thresh=thresholds),
    # expand("output/{sample}/{species}.{filt}.bam",species=method_list,sample=config,filt=filter_list),
    # expand("output/{sample}/{species}.{filt}.vcf",species=method_list,sample=config, filt=filter_list),
    expand("output/{sample}/performance.csv",species=method_list,sample=config)
    # expand("output/{sample}/{species}.{filt}.indel.vcf",species=method_list,sample=config, filt=filter_list),
    # expand("refs/{sample}.assembly",sample=config),
    # expand("output/{sample}/{species}.{filt}.bam",species=method_list,sample=config,filt=filter_list),
    # expand("output/{sample}/{species}.{threshold}.{filt}.bam.bai",species=method_list,sample=config,filt=filter_list,threshold=thresholds),
    # # expand("output/{sample}/dRep_assembly/",sample=config)
    # expand("output/{sample}/{species}.{threshold}.{filt}.vcf", \
    #     sample=config, species=method_list, mapper=('lofreq',), \
    #     filt=filter_list,threshold=thresholds), 
    # # expand("output/{sample}/{species}_profile_IS", sample=config, species=method_list),
    # expand("output/{sample}/compare_{species}.{filt}.{mapper}.tsv", \
    #     sample=config, species=method_list, mapper=('lofreq',), \
    #     filt=filter_list, isofilt=filter_list),
    # # expand("output/{sample}/{species}_instrain.lofreq.tsv",sample=config,species=method_list, filt=filter_list)
    # expand("output/{sample}/compare_result.tsv", sample=config),

    # expand("output/{sample}/compare_{species}.{filt}.iso_{isofilt}.{mapper}.filt.tsv", sample=config, species=method_list, mapper=("lofreq",),filt=filter_list,isofilt=filter_list,),
    # expand("output/{sample}/compare_result_filtered.tsv", sample=config)

#################################################################################################3
## Idea:
#   learn from dRep and have a script that checks the dependencies are properly installed
#   and list down their locations
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
        "refs/{sample}.{species}.fasta"
    output:
        "refs/{sample}.{species}.fasta.1.bt2",
        "refs/{sample}.{species}.fasta.2.bt2",
        "refs/{sample}.{species}.fasta.3.bt2",
        "refs/{sample}.{species}.fasta.4.bt2",
        "refs/{sample}.{species}.fasta.rev.1.bt2",
        "refs/{sample}.{species}.fasta.rev.2.bt2"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "bowtie2-build {input} {input}"

rule alignment_bowtie:
    input:
        r1 = lambda wc: config[wc.sample]['r1'],
        r2 = lambda wc: config[wc.sample]['r2'],
        ref = "refs/{sample}.{species}.fasta",
        idx = rules.bowtie_build.output
    output:
        tmp = temp("output/{sample}/{species}_tmp.bam"),
        sec = temp("output/{sample}/{species}_secondary.bam"),
        bam = 'output/{sample}/{species}_bowtie.bam'
    conda:
        'workflow/envs/mapping.yaml'
    threads:
        8
    shell:
        "bowtie2 -a -p {threads} -x {input.ref} -1 {input.r1} -2 {input.r2} | samtools sort -o {output.tmp} -@ {threads} && "
        "samtools view -f 100 -bh {output.tmp} > {output.sec} && samtools merge {output.bam} {output.tmp} {output.sec}"

rule assembly:
    input:
        r1 = lambda wc: config[wc.sample]['r1'],
        r2 = lambda wc: config[wc.sample]['r2']
    output:
        fa = "output/{sample}/assembly/final.contigs.fa"
    benchmark:
        "benchmarks/{sample}/megahit_resource_usage.txt"
    threads:
        16
    shell:
        "megahit -1 {input.r1} -2 {input.r2} -f -o output/{wildcards.sample}/assembly -t {threads}"

rule mummmer_assembly_to_ref:
    input:
        ref = lambda wc: config[wc.sample]['ref'],
        assembly = rules.assembly.output.fa
    output:
        "output/{sample}/assembly/mummer.mums"
    shell:
        "mummer {input.ref} {input.assembly} > {output}"

rule get_nonref_contig:
    input:
        assembly = "output/{sample}/assembly/final.contigs.fa",
        mummer =  rules.mummmer_assembly_to_ref.output
    output:
        temp("output/{sample}/assembly/passed_contigs.txt")
    shell:
        "python refs/filter_contigs.py {input.assembly} {input.mummer} > {output}"

rule filter_assembly:
    input:
        nonrefcontig = rules.get_nonref_contig.output,
        assembly =  "output/{sample}/assembly/final.contigs.fa",
        ref = lambda wc: config[wc.sample]['ref'],
        abund = "output/{sample}/abundant_species.fasta"
    output:
        "output/{sample}/assembly.fasta"
    shell:
        "grep -w -A 1 -f {input.nonrefcontig} {input.assembly} > {output} && cat {input.abund} >> {output}"

ruleorder:  align_large > alignment

# rule align_assembly:
#     input:
#     output:

# rule binning:
#     input:
#         ref = "refs/{sample}.assembly/final.contigs.fa",
#         bam = rules.align_assembly.output.bam
#     output:
#         directory("output/{sample}/bins/")
#     params:
#         depth = "output/{sample}/assembly_alignment_depth.txt",
#         outdir = "output/{sample}/bins/bin"
#     shell:
#         "jgi_summarize_bam_contig_depths --outputDepth {params.depth} {input.bam} && "
#         "metabat2 -m 1500 -i {input.ref} -a {params.depth} --outFile {params.outdir}"


rule drep:
    input:
        "output/{sample}/bins/"
    output:
        directory("output/{sample}/dRep_assembly/")
    shell:
        "dRep dereplicate {output} -g {input}/*.fa"
    
# rule database_n_assembly: # not actively in use
#     input: 
#         ref = "refs/{sample}.assembly/final.contigs.fa",
#         database = "metaSNV/progenomes2_speciesReps_genomes.fna"
#     output:
#         ref = "output/{sample}/database_n_assembly.fna"
#     shell:
#         "cp {input.database} {output.ref} && cat {input.ref} >> {output.ref}"


rule index_db:
    input:
        db = "gtdb/uhgg_reps.fasta.gz" # db = "metaSNV/progenomes2_speciesReps_genomes.fna"
    output:
        "gtdb/database.mmi"
    shell:
        "minimap2 -x sr -t 10 -d {output} {input.db}"

# rule index_db_assembly: # not actively in use
#     input:
#         "output/{sample}/database_n_assembly.fna"
#     output:
#         "output/{sample}/database_n_assembly.mmi"
#     shell:
#         "minimap2 -x sr -t 10 -d {output} {input}"

rule align_large:
    input:
        r1 = lambda wc: config[wc.sample]['r1'],
        r2 = lambda wc: config[wc.sample]['r2'],
        ref = "gtdb/database.mmi"
    output:
        bam = temp('output/{sample}/database.bam')
    conda:
        'workflow/envs/mapping.yaml' # include minimap2
    benchmark:
        "benchmarks/{sample}/align_large_resource_usage.txt"
    params:
        lambda wc: wc.sample
    threads:
        16
    shell:
        "minimap2 -ax sr -t {threads} --split-prefix=foo {input.ref} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output.bam}"

rule alignment: 
    input:
        r1 = lambda wc: config[wc.sample]['r1'],
        r2 = lambda wc: config[wc.sample]['r2'],
        ref = "output/{sample}/{species}.fasta"
    output:
        bam = temp('output/{sample}/{species}.bam')
    conda:
        'workflow/envs/mapping.yaml' # include minimap2
    threads:
        12
    shell:
        "minimap2 -t {threads} -ax sr --secondary=yes {input.ref} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output.bam}"
# rule alignment_isolate:
#     input:
#         r1 = "/mnt/volume1/isolate/{iso_label}/merge_{iso_label}_R1.fastq.gz" ,
#         r2 = "/mnt/volume1/isolate/{iso_label}/merge_{iso_label}_R2.fastq.gz" ,
#         ref = "refs/{species}.fasta",
#         idx = rules.bowtie_build.output
#     output:
#         temp("output/{iso_label}/{species}.bam")
#     conda:
#         'workflow/envs/mapping.yaml'
#     threads:
#         8
#     shell:
#         "bowtie2 -a -p {threads} -x {input.ref} -1 {input.r1} -2 {input.r2} | "
#         "samtools sort -o {output} -@ {threads}"

rule samtools_index: 
    input:
        "output/{sample}/{placeholder}.bam"
    output:
        'output/{sample}/{placeholder}.bam.bai'
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        'samtools index {input}'
# rule samtools_index_filtered:
#     input:
#         "output/{sample}/{species}.{filt}.bam"
#     output:
#         'output/{sample}/{species}.{filt}.bam.bai'
#     conda:
#         'workflow/envs/mapping.yaml'
#     shell:
#         'samtools index {input}'

# rule samtools_index_filtered_sorted:
#     input:
#         "output/{sample}/{species}.{filt}.sorted.bam"
#     output:
#         'output/{sample}/{species}.{filt}.sorted.bam.bai'
#     conda:
#         'workflow/envs/mapping.yaml'
#     shell:
#         'samtools index {input}'

rule bamfilter:
    input:
        "output/{sample}/{species}.bam",
        "output/{sample}/{species}.bam.bai"
    output:
        "output/{sample}/{species}.{filt}.bam"
    conda:
        'workflow/envs/mapping.yaml'
    params:
        "{filt}",
        lambda wc: config[wc.sample]['ref']
    script:
        'workflow/scripts/bamfilter.py'

rule sort_bam:
    input: 
        bam = "output/{sample}/{species}.{filt}.bam"
    output:
        "output/{sample}/{species}.{filt}.sorted.bam"
    shell:
        "samtools sort {input.bam} -o {output}"

rule lofreq:
    input:
        ref = lambda wc: config[wc.sample]['ref'],
        filt_bam = "output/{sample}/{species}.{filt}.sorted.bam",
        idx = "output/{sample}/{species}.{filt}.sorted.bam.bai"
    output:
        "output/{sample}/{species}.{filt}.vcf"
    conda:
        'workflow/envs/mapping.yaml'
    threads:
        8
    shell: 
        "lofreq faidx {input.ref} && lofreq call-parallel --pp-threads {threads} -f {input.ref} -o {output} --verbose {input.filt_bam}"


rule get_lofreq_indel_qual:
    input:
        ref = lambda wc: config[wc.sample]['ref'],
        filt_bam = "output/{sample}/{species}.{filt}.bam",
        idx = "output/{sample}/{species}.{filt}.bam.bai"
    output:
        "output/{sample}/{species}.{filt}.indel.bam"
    conda:
        'workflow/envs/mapping.yaml'
    shell:
        "lofreq faidx {input.ref} && lofreq indelqual --dindel -f {input.ref} -o {output} {input.filt_bam}"


rule lofreq_indel:
    input:
        ref = lambda wc: config[wc.sample]['ref'],
        bam = rules.get_lofreq_indel_qual.output
    output:
        "output/{sample}/{species}.{filt}.indel.vcf"
    conda:
        'workflow/envs/mapping.yaml'
    threads:
        8
    shell: 
        "lofreq call --call-indels -f {input.ref} -o {output} --verbose {input.bam}"
        # "lofreq faidx {input.ref} && lofreq call-parallel --pp-threads {threads} -f {input.ref} -o {output} --verbose {input.filt_bam}"

rule pileup:
    input:
        ref = "refs/{species}.fasta",
        bam = rules.alignment.output.bam
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


rule process_result:
    input:
        expand("output/{sample}/{method}.filt.vcf", method=method_list, sample=config),
        expand("output/{sample}/{method}.filt.sorted.bam", method=method_list,sample=config)
    output:
        "output/{sample}/performance.csv",
        "output/{sample}/read_mapping_perf.jpg",
        "output/{sample}/SNP_perf.jpg",
        "output/{sample}/FP_reads.json"
    conda:
        'workflow/envs/plotting.yaml'
    params:
        samp = "{sample}",
        methods = list(method_list),
        ref = lambda wc: config[wc.sample]['ref'],
        snp_lst = lambda wc: config[wc.sample]['isolate']
    script:
        "workflow/scripts/process_vcf.py"

# rule compare_instrain:
#     input:
#         source = "output/{sample}/{species}_profile_IS",
#         meta = "output/{sample}/{species}_profile_IS/output/{species}_profile_IS_SNVs.tsv",
#         iso = lambda wc: config[wc.sample]["isolate"]
#         # iso = lambda wc: os.path.join("output/", config[wc.sample]['isolate'] + "/Escherichia_coli_iai39." + wc.isofilt + "." + wc.mapper + ".vcf")
#     output:
#         "output/{sample}/{species}_instrain.{mapper}.tsv"
#     conda:
#         'workflow/envs/mapping.yaml'
#     params:
#         "instrain"
#     script:
#         "workflow/scripts/JS_script.py"

rule prodigal:
    input:
        ref =lambda wc: config[wc.species]["ref"]
    output:
        "refs/{sample}.{species}.genes.fna"
    shell:
        "prodigal -i {input.ref} -d {output}"


rule instrain:
    input:
        genes = rules.prodigal.output,
        ref = lambda wc: config[wc.species]["ref"] ,
        bam = 'output/{sample}/{species}.nofilt.bam' # cuz inStrain does its own filtering, right?
    output:
        directory("output/{sample}/{species}_profile_IS")
    conda:
        "workflow/envs/instrain.yaml"
    benchmark:
        "benchmarks/{sample}/instrain_resource_usage_{species}.txt"
    threads:
        8
    shell:
        "inStrain profile --min_cov 10 {input.bam} {input.ref} -g {input.genes} -o {output} -p 8"

# rule compare_metrics:
#     input:
#         # result_files = glob_wildcards("output/{sample}/")
#         # expand(rules.compare_instrain.output, sample=config,species=method_list,isofilt=filter_list,mapper="lofreq"),
#         expand("output/{sample}/compare_{species}.{filt}.lofreq.tsv", species=method_list,filt=filter_list,sample = "{sample}")
#     output:
#         "output/{sample}/compare_result.tsv"
#     params:
#         sample = "{sample}",
#         species = list(method_list),
#         filt = False
#     conda:
#         'workflow/envs/mapping.yaml'
#     script:
#         "workflow/scripts/combined_tsv.py"

# need to add rules for downloading database, unzip and create signatures

rule get_abund_species:
    input:
        r1 = lambda wc: config[wc.sample]['r1'],
        db = "gtdb/gtdb_sketch"
    output:
        "output/{sample}/abundant_species.csv"
    params:
        sig = "output/{sample}/{sample}.sig"
    conda:
        "workflow/envs/sim.yaml"
    shell:
        "sourmash sketch dna -p abund {input.r1} -o {params.sig} &&"
        "sourmash gather {params.sig} {input.db} -o {output} --threshold=0.0"

rule get_close_species:
    input:
        ref = lambda wc: config[wc.sample]['ref'],
        db = "gtdb/gtdb_sketch"
    output:
        "output/{sample}/close_species.csv"
    conda:
        "workflow/envs/sim.yaml"
    params:
        sig = "output/{sample}/close_species.sig"
    shell:
        "sourmash sketch dna {input.ref} -o {params.sig} && "
        "sourmash search {params.sig} {input.db} -o {output} --threshold=0.0 --ignore-abundance"

rule build_metagenome_abund:
    input:
        csv = rules.get_abund_species.output,
        ref = lambda wc: config[wc.sample]['ref'],
        close = rules.get_close_species.output
    output:
        "output/{sample}/abundant_species.fasta"
    params:
        species = lambda wc: config[wc.sample]["species"],
        sample = lambda wc: wc.sample
        # t = lambda wc: wc.threshold
    shell:
        "./workflow/scripts/combine_fasta.sh {input.csv} {output} {input.ref} '{params.species}' {input.close}"

rule build_metagenome_close:
    input:
        csv = rules.get_close_species.output,
        ref = lambda wc: config[wc.sample]['ref']
    output:
        "output/{sample}/close_species.fasta"
    params:
        lambda wc: config[wc.sample]["species"]
    shell:
        "./workflow/scripts/combine_fasta_close.sh {input.csv} {output} {input.ref} '{params.species}' "

rule build_target_reference:
    input:
        ref = lambda wc: config[wc.sample]['ref']
    output:
        "output/{sample}/target_species.fasta"
    shell:
        "cat {input.ref} > {output}"

rule build_consensus_reference: # don't think it was ever in use
    input:
        ref = lambda wc: config[wc.sample]['ref'],
        vcf = "output/{sample}/target_species.filt.vcf"
    output:
        "output/{sample}/consensus.fasta"
    shell:
        "python workflow/scripts/consensus_ref.py {input.ref} {input.vcf} {output}"
