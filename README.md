# FR-AgENCODE ATAC-seq pipeline
This is the pipeline to process ATAC-seq data from fastq to peaks used in the FR-AgENCODE project http://www.fragencode.org/

## Download the code
<pre>
git clone https://github.com/sdjebali/fragencode-atacseq-pipeline.git
</pre>

## Installation
Use conda
<pre>
conda env create --name atacseq --file environment.yaml
</pre>

If everything went fine activate the environment
<pre>
conda activate atacseq
conda activate base
</pre>

And make sure you have java8 and bowtie2-2.3.5.1 installed, for example using module load (genologin)
<pre>
module load system/Java8
module load bioinfo/bowtie2-2.3.5.1
</pre>


## Usage
**************************************
Make sure samples.tsv, reads.tsv and config.yaml files have the correct content

## Caution
The genome files must be indexed with bowtie2 before running the pipeline

## Running snakemake
<pre>
export DRMAA_LIBRARY_PATH="/tools/libraries/slurm-drmaa/slurm-drmaa-1.0.7/lib/libdrmaa.so"
snakemake --use-conda --debug-dag --jobs 30 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" --configfile config.yaml -n -p
</pre>
