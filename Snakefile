
config = {
  'WEB1352': {
    'r1': '/mnt/volume1/isolate/WEB1352/WEB1352-GCACGCGT_S7_L004_R1_001_HS006-PE-R00207.fastq.gz',
    'r2': '/mnt/volume1/isolate/WEB1352/WEB1352-GCACGCGT_S7_L004_R2_001_HS006-PE-R00207.fastq.gz'
  }
}

rule all:
  input:
    expand("/mnt/volume1/result/{sample}/ec.bam", sample=config)

rule isolate_alignment:
  input:
    r1 = lambda wc: config[wc.sample]['r1'] ,
    r2 = lambda wc: config[wc.sample]['r2'] ,
    ref = "mnt/volume1/ref/E_coli_ref.fna"
  output:
    bam_file = "/mnt/volume1/isolate/{sample}/ec.bam"
  shell:
    "bowtie2 -p 6 -x {input.ref} -1 {input.r1} -2 {input.r2} | samtools view -@ 16 -b -u - | samtools sort -o {output.bam_file}"
    "samtools index {output.bam_file}"

rule sample_alignment:
  input:
    r1 = "/mnt/volume1/cpe/{sample}/merge_{sample}_R1.fastq.gz",
    r2 = "/mnt/volume1/cpe/{sample}/merge_{sample}_R2.fastq.gz",
    ref = "/mnt/volume1/ref/E_coli_ref.fna"
  output:
    bam_file = "/mnt/volume1/cpe/{sample}/ec.bam"
  shell:
    "bowtie2 -p 6 -x {input.ref} -1 {input.r1} -2 {input.r2} | "
    "samtools view -@ 16 -b -u - | samtools sort -o {output.bam_file} && "
    "samtools index {output.bam_file}"

rule filter:
  input:
    "/mnt/volume1/cpe/{sample}/ec.bam"
  output:
    "/mnt/volume1/cpe/{sample}/ec_filt.bam"
  shell:
    "python /mnt/volume1/scripts/bamfilter.py {input} {output}"

rule lofreq:
  input:
    ref = "/mnt/volume1/ref/E_coli_ref.fna"
    filt_bam = "/mnt/volume1/cpe/{sample}/ec_filt.bam"
  output:
    "/mnt/volume1/result/{sample}/filt.bam.vcf"
  shell:
    "lofreq call -f {input.ref} -o {output} --verbose {output}"

rule pileup:
  input:
    ref = "/mnt/volume1/ref/E_coli_ref.fna"
    bam = "/mnt/volume1/result/{sample}/ec.bam"
  output:
    pileup = temp("/mnt/volume1/result/{sample}/ec.bam.pileup")
    gz = "/mnt/volume1/result/{sample}/ec.bam.pileup.gz"
  shell:
    "samtools mpileup -E -f {input.ref} {input.bam} > {output.pileup} &&"
    "gzip {output.pileup} && rm {output.pileup}"

rule coverage:
  input:
    "/mnt/volume1/cpe/{sample}/ec_filt.bam"
  output:
    txt = "/mnt/volume1/result/{sample}/coverage.txt"
    summary = "/mnt/volume1/result/{sample}/coverage.summary"
  shell:
    "genomeCoverageBed -d -ibam {input} > {output.txt} && "
    "samtools depth -aa {input} | "
    "awk '{sum+=$3} END {print sum/NR}' > {output.summary} "

rule annotate_gene:
  input:
    "/mnt/volume1/result/{sample}/filt.bam.vcf"
  output:
     "/mnt/volume1/result/{sample}/ann.vcf"
  shell:
    "/mnt/volume1/scripts/annotate.sh EC {input} {output}"


