import os
from snakemake.io import glob_wildcards
from snakemake.utils import min_version
min_version("6.0")

### Default parameters
## --------------------------------------------------------------------------
method_list = ("target_species", "abundant_species", "database", "assembly")
filter_list = ("filt", "unfilt")


wildcard_constraints:
    species="\w+",
    filt="\w+",


### Load conda configurations
## --------------------------------------------------------------------------

config['softparams'] = {}
fname = workflow.source_path('conda_env.json')
with open(fname) as f:
    config['softparams']['conda'] = json.load(f)

### Load samples
## --------------------------------------------------------------------------

samples = config["samples"]
for sample in samples:
    config["samples"][sample]["r1"] = f"simulation/illumina/{sample}/merged_R1.fq.gz"
    config["samples"][sample]["r2"] = f"simulation/illumina/{sample}/merged_R2.fq.gz"

    ref = config["samples"][sample]["ref"]
    config["samples"][sample]["ref"] = (
        ref if os.path.isfile(ref) else f"simulation/download/{ref}.fa"
    )
# --------------------------------------------------------------------------------------------------------------------------------------------


rule all:
    input:
        # expand("simulation/fastANI/{sample}.fastani.tsv", sample=samples),
        expand(
            "assessment/{sample}/mapping_assessment.{kind}.tsv",
            sample=samples,
            kind=(
                "wodecoys.fixmate.sorted.filtered",
                "decoys",
                "decoys.fixmate.sorted.filtered",
            ),
        ),
        # expand(
        #     "assessment/{sample}/vcalling_assessment.{kind}.tsv",
        #     sample=samples,
        #     kind=("wodecoys", "decoys"),
        # )

###############################################################################
module simulate:
    snakefile:
        "workflow/simulate.smk"
    config:
        config


use rule * from simulate as sim_*


module sourmash:
    snakefile:
        "workflow/sourmash.smk"
    config:
        config


use rule * from sourmash as sourmash_*


module mapping:
    snakefile:
        "workflow/mapping.smk"
    config:
        config


use rule * from mapping as mapping_*


module assessment:
    snakefile:
        "workflow/assessment.smk"
    config:
        config


use rule * from assessment as assessment_*
