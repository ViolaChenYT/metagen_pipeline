## Why is it called "simulate"


rule ncbi_ass: # download NCBI database??? reference?
    ''' download the specific NCBI reference for particular genome ID '''
    input:
        script = workflow.source_path('scripts/download_ncbi.sh')
    output:
        'simulation/download/{genomeid}.fa'
    conda:
        config['softparams']['conda']['entrez']
    resources:
        ncbi_load=1
    params:
        config['ncbi_api']
    shell:
        'bash {input.script} {wildcards.genomeid} {output} {params} && [[ -s {output} ]]'

rule mutate:
    ''' create mutated version of the sequence''' 
    input:
        rules.ncbi_ass.output
    output:
        'simulation/mutated/{sample}/{genomeid}.fa'
    conda:
        config['softparams']['conda']['biopython']
    params:
        mfreq = lambda wc: config['samples'][wc.sample]['genomes'][wc.genomeid][0] / 100
    script:
        'scripts/simulation_mutate.py' # has to be changed, i am not implementing this mutation

rule illumina:
    ''' art == illumina 454 and solid read simulator
    - why specifically art?
    - i can totally do this myself as well right????'''
    input:
        rules.mutate.output
    output:
        r1 = temp('simulation/illumina/{sample}/{genomeid}_R1.fq'),
        r2 = temp('simulation/illumina/{sample}/{genomeid}_R2.fq')
    params:
        instrument = 'HS25',
        seed = 1,
        length = 150,
        meansize = 300,
        std = 10,
        coverage = lambda wc: config['samples'][wc.sample]['genomes'][wc.genomeid][1],
        outname = 'simulation/illumina/{sample}/{genomeid}_R'
    conda:
        config['softparams']['conda']['art']
    log:
        'simulation/illumina/{sample}/{genomeid}.log'
    shell:
        'art_illumina -ss {params.instrument} -i {input} -p -l {params.length} -f {params.coverage} -m {params.meansize}\
         -s {params.std} -o {params.outname} -d "::{wildcards.genomeid}::" -na -rs {params.seed} > {log}'

rule merge_reads:
    ''' put reads together and gzip them
    '''
    input:
        lambda wc: [f'simulation/illumina/{wc.sample}/{genomeid}_{wc.rtype}.fq'
            for genomeid in config['samples'][wc.sample]['genomes']]
    output:
        'simulation/illumina/{sample}/merged_{rtype}.fq.gz'
    shell:
        'cat {input} | gzip > {output}'

rule fastANI_files_list:
    input:
        lambda wc: [f'simulation/download/{genomeid}.fa'
            for genomeid in config['samples'][wc.sample]['genomes']]
    output:
        temp('simulation/fastANI/{sample}.genomes.txt')
    shell:
        'ls -d {input} > {output}'

rule fastANI:
    ''' seem to be batch running fastANI????
    '''
    input:
        rules.fastANI_files_list.output
    output:
        'simulation/fastANI/{sample}.fastani.tsv'
    conda:
        config['softparams']['conda']['fastANI']
    log:
        'simulation/fastANI/{sample}.fastani.log'
    shell:
        'fastANI --ql {input} --rl {input} -o {output} 2> {log}'