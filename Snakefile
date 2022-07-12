rule all:
  input:
    expand("/mnt/volume1/result/{sample}/ec.bam", sample=config)

rule bowtie_build:
  input:
    "data/{species}_ref.fna"
  output:
    "data/{species}_ref.fna.1.bt2",
    "data/{species}_ref.fna.2.bt2",
    "data/{species}_ref.fna.3.bt2",
    "data/{species}_ref.fna.4.bt2",
    "data/{species}_ref.fna.rev.1.bt2",
    "data/{species}_ref.fna.rev.2.bt2"
  shell:
    "bowtie2-build {input} {input}"

rule isolate_alignment:
  input:
    r1 = lambda wc: config[wc.sample]['r1'] ,
    r2 = lambda wc: config[wc.sample]['r2'] ,
    ref = "data/E_coli_ref.fna"
  output:
    bam_file = "{homedir}/isolate/{sample}/ec.bam"
  conda:
    'workflow/envs/mapping.yaml'
  shell:
    "bowtie2 -p 6 -x {input.ref} -1 {input.r1} -2 {input.r2} | "
    "samtools view -@ 16 -b -u - | "
    "samtools sort -o {output.bam_file}"


rule alignment:
  input:
    r1 = "{sample}/merge_{sample}_R1.fastq.gz",
    r2 = "{sample}/merge_{sample}_R2.fastq.gz",
    ref = "data/E_coli_ref.fna",
    bt1 = "data/E_coli_ref.fna.1.bt2",
    bt2 = "data/E_coli_ref.fna.2.bt2",
    bt3 = "data/E_coli_ref.fna.3.bt2",
    bt4 = "data/E_coli_ref.fna.4.bt2",
    rev1 = "data/E_coli_ref.fna.rev.1.bt2",
    rev2 = "data/E_coli_ref.fna.rev.2.bt2"
  output:
    bam_file = "{sample}/ec.bam",
  shell:
    "bowtie2 -p 6 -x {input.ref} -1 {input.r1} -2 {input.r2} | "
    "samtools view -@ 16 -b -u - | samtools sort -o {output.bam_file}"


rule samtools_index:
  input:
    '{sample}/ec.bam'
  output:
    '{sample}/ec.bam.bai'
  shell:
    'samtools index {input}'

rule filter:
  input:
    "{sample}/ec.bam",
    "{sample}/ec.bam.bai"
  output:
    "{sample}/ec_filt.bam"
  script:
    'workflow/scripts/bamfilter.py'

rule lofreq:
  input:
    ref = "data/E_coli_ref.fna",
    filt_bam = "{homedir}/cpe/{sample}/ec_filt.bam"
  output:
    "{homedir}/result/{sample}/filt.bam.vcf"
  shell:
    "lofreq call -f {input.ref} -o {output} --verbose {input.filt_bam}"

rule lofreq_test:
  input:
    ref = "data/E_coli_ref.fna",
    filt_bam = "{dir}/ec_filt.bam"
  output:
    "{dir}/output/filt.bam.vcf"
  shell:
    "lofreq call -f {input.ref} -o {output} --verbose {input.filt_bam}"

rule lofreq_isolate:
  input:
    ref = "data/E_coli_ref.fna",
    filt_bam = "{homedir}/isolate/{sample}/ec_filt.bam"
  output:
    "{homedir}/result/{sample}/filt.bam.vcf"
  shell:
    "lofreq call -f {input.ref} -o {output} --verbose {input.filt_bam}"

rule pileup:
  input:
    ref = "data/E_coli_ref.fna",
    bam = "{homedir}/result/{sample}/ec.bam"
  output:
    gz = "{homedir}/result/{sample}/ec.bam.pileup.gz"
  shell:
    "samtools mpileup -E -f {input.ref} {input.bam} | gzip > {output.gz}"

rule pileup_test:
  input:
    ref = "data/E_coli_ref.fna",
    bam = "{homedir}/ec.bam"
  output:
    gz = "{homedir}/output/ec.bam.pileup.gz"
  shell:
    "samtools mpileup -E -f {input.ref} {input.bam} | gzip > {output.gz}"

rule coverage_txt:
  input:
    "{homedir}/cpe/{sample}/ec_filt.bam"
  output:
    "{homedir}/result/{sample}/coverage.txt"
  shell:
    "genomeCoverageBed -d -ibam {input} > {output}"

rule coverage_avg:
  input:
    "{homedir}/cpe/{sample}/ec_filt.bam"
  output:
    "{homedir}/result/{sample}/coverage.summary"
  shell:
    "samtools depth -aa {input} | awk '{sum+=$3} END {print sum/NR}' > {output}"

rule coverage_txt_test:
  input:
    "{homedir}/ec_filt.bam"
  output:
    "{homedir}/output/coverage.txt"
  shell:
    "genomeCoverageBed -d -ibam {input} > {output}"

rule coverage_avg_test:
  input:
    "{homedir}/ec_filt.bam"
  output:
    "{homedir}/output/coverage.summary"
  shell:
    "samtools depth -aa {input} | awk '{sum+=$3} END {print sum/NR}' > {output}"

rule annotate_gene:
  input:
    "{homedir}/{dir}/{sample}/filt.bam.vcf"
  output:
     "{homedir}/{dir}/{sample}/ann.vcf"
  shell:
    "workflow/scripts/annotate.sh EC {input} {output}"

rule annotate_gene_test:
  input:
    "{homedir}/{dir}/filt.bam.vcf"
  output:
     "{homedir}/{dir}/ann.vcf"
  shell:
    "workflow/scripts/annotate.sh EC {input} {output}"


