

#################################
# Trim reads to remove adaptors #
#################################
# Use cutadapt version 2.1
# module load bioinfo/cutadapt-2.1
# Use latest version of trim_galore modified on July 2019
# Takes 2h30 with 1 cpu and 0.5 G used
rule trimming:
    input:
        read1=get_fastq1,
        read2=get_fastq2
    output:
        #outdir="{tiss}/{anim}/trimming/{sample}",
        outfile="{tiss}/{anim}/trimming/{sample}/{sample}_trim_galore.out",
        read1="{tiss}/{anim}/trimming/{sample}/{sample}_R1_val_1.fq.gz",
        read2="{tiss}/{anim}/trimming/{sample}/{sample}_R2_val_2.fq.gz"
    log:
        "logs/trimming/{tiss}/{anim}/{sample}.log"
    shell:
        "outdir=`dirname {output.outfile}`;"
        "trim_galore --basename {wildcards.sample} --stringency 3 -q 20 --paired --nextera"
        " -o $outdir {input.read1} {input.read2} > {output.outfile} 2> {log}"
