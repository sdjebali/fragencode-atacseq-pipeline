import pandas as pd
from snakemake.utils import validate

include: "rules/common.smk"

# path to snakefile
snakedir = os.path.dirname(workflow.snakefile)

##### load config and read file #####
configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["sample_file"]).set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

reads = pd.read_table(config["read_file"]).set_index("sample", drop=False)
#validate(reads, schema="schemas/reads.schema.yaml")


wildcard_constraints:
    sample="\w+"

##--------------------------------------
##  Global variables
##--------------------------------------
WORKDIR = config["workdir"]
workdir: WORKDIR
message("The current working directory is "+WORKDIR)


###############
# main rule   #
###############
rule all:
    input:
        expand("{tiss}/{anim}/peaks/{sample}.q10.noMT.nodup.macs_peaks.narrowPeak", sample=samples.index, tiss=samples["tissue"], anim=samples["animal"])

#### load rules #####
include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/filtering.smk"
include: "rules/peakcalling.smk"
