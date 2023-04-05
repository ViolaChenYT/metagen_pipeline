from utils import get_db_config
import pandas as pd

SBT_LINK = 'https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.sbt.zip'
ZIP_LINK = 'https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.zip'

rule download_sbt:
    output:
        protected('sourmash/db/gtdb-rs207.genomic-reps.dna.k31.sbt.zip')
    params:
        link = SBT_LINK
    shell:
        'curl {params.link} > {output}'

rule download_zip:
    output:
        protected('sourmash/db/gtdb-rs207.genomic-reps.dna.k31.zip')
    params:
        link = ZIP_LINK
    shell:
        'curl {params.link} > {output}'

rule sourmash_gather_signature:
    input:
        r1 = lambda wc: config['samples'][wc.sample]['r1']
    output:
        temp("sourmash/samples/{sample}/{sample}.sig")
    conda:
        config['softparams']['conda']['sourmash']
    shell:
        "sourmash sketch dna -p abund {input.r1} -o {output} --force"

rule sourmash_gather:
    input:
        sig = rules.sourmash_gather_signature.output,
        db = get_db_config(config, 'sourmash_sbt', 'sourmash/db/gtdb-rs207.genomic-reps.dna.k31.sbt.zip')
    output:
        "sourmash/samples/{sample}/gather.csv"
    conda:
        config['softparams']['conda']['sourmash']
    shell:
        "sourmash gather {input.sig} {input.db} -o {output} --threshold=0.0 --quiet"

rule sourmash_search_signature:
    input:
        lambda wc: config['samples'][wc.sample]['ref']
    output:
        temp("sourmash/samples/{sample}/close_species.sig")
    conda:
        config['softparams']['conda']['sourmash']
    shell:
        'sourmash sketch dna {input} -o {output} --force'

rule sourmash_search:
    input:
        sig = rules.sourmash_search_signature.output,
        db = get_db_config(config, 'sourmash_zip', 'sourmash/db/gtdb-rs207.genomic-reps.dna.k31.zip')
    output:
        "sourmash/samples/{sample}/refsearch.csv"
    conda:
        config['softparams']['conda']['sourmash']
    shell:
        "sourmash search {input.sig} {input.db} -o {output} --threshold=0.0 --ignore-abundance --quiet"

checkpoint sourmash_merge_filter:
    input:
        gather = rules.sourmash_gather.output,
        search = rules.sourmash_search.output
    output:
        "sourmash/samples/{sample}/sourmash_res_filtered.csv",
    script:
        'scripts/sourmash_merge_filter.py'

rule empty_file:
    output:
        temp('sourmash/res_gtdb/empty.txt.gz')
    shell:
        'touch sourmash/res_gtdb/empty.txt && gzip sourmash/res_gtdb/empty.txt'

rule download_genome:
    input:
        script = workflow.source_path('scripts/download_ncbi.sh')
    output:
        'sourmash/res_gtdb/download/{gid}.fa'
    conda:
        config['softparams']['conda']['entrez']
    resources:
        ncbi_load=1
    shell:
        'bash {input.script} {wildcards.gid} {output} && [[ -s {output} ]]'

rule compress_genome:
    input:
        rules.download_genome.output
    output:
        'sourmash/res_gtdb/download/{gid}.fa.gz'
    shell:
        'gzip {input}'

def get_gtdb_genomes(wc):
    fname = checkpoints.sourmash_merge_filter.get(** wc).output[0]
    gtdb_gids = list(pd.read_csv(fname)['gtdb_gid'])
    if not gtdb_gids:
        return 'sourmash/res_gtdb/empty.txt.gz'
    return [
        f'sourmash/res_gtdb/download/{gid}.fa.gz'
        for gid in gtdb_gids
    ]

rule merge_download_genomes:
    input:
        get_gtdb_genomes
    output:
        temp('sourmash/res_gtdb/samples/{sample}.merged.fa.gz')
    shell:
        'cat {input} > {output}'

rule contig_list:
    input:
        get_gtdb_genomes
    output:
        'sourmash/res_gtdb/samples/{sample}.merged.clist.txt'
    shell:
        '(zgrep -H ">" {input} > {output}) || true'
        # I don't like it but it's the only way to make it work
        # with empty file