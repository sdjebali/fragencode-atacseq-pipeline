

#####################
# Map trimmed reads #
#####################
# Use bowtie2 version 2.3.5.1
# module load bioinfo/bowtie2-2.3.5.1
# Takes several hours with 4 cpus and 25.5G of ram
rule map:
    input:
        read1="{sp}/{tiss}/{anim}/trimming/{sample}/{sample}_R1_val_1.fq.gz",
        read2="{sp}/{tiss}/{anim}/trimming/{sample}/{sample}_R2_val_2.fq.gz"
    output:
        "{sp}/{tiss}/{anim}/mapping/{sample}.sam"
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.map.log"
    threads:
        4
    params:
    # note that the following is not generic since we want to apply it to any species not just to sus_scrofa
        genome=config["genome"]["sus_scrofa"]
    shell:
        "bowtie2 -t -p {threads} -X 2000 -x {params.genome}"
        " -1 {input.read1} -2 {input.read2} --met 60 --met-file {log}"
        " -S {output}"

###################
# Convert to bam  #
###################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule samtobam:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.sam"
    output:
        temp("{sp}/{tiss}/{anim}/mapping/{sample}.unsorted.bam")
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.tobam.log"
    threads:
        4
    shell:
        "samtools view {input} -Sbh -@ {threads} -o {output} 2> {log}"

#####################
# Sort the bam file #
#####################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule bamsort:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.unsorted.bam"
    output:
        "{sp}/{tiss}/{anim}/mapping/{sample}.bam"
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.bamsort.log"
    threads:
        4
    shell:
        "samtools sort -@ {threads} -m 4G {input} -O bam -o {output} 2> {log}"


######################
# Index the bam file #
######################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule bamindex:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.bam"
    output:
        "{sp}/{tiss}/{anim}/mapping/{sample}.bam.bai"
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.bamidx.log"
    shell:
        "samtools index {input} 2> {log}"


#########################
# Flagstat the bam file #
#########################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule bamflag:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.bam"
    output:
        "{sp}/{tiss}/{anim}/mapping/{sample}.bam.flagstat"
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.bamflag.log"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"

############################################################
# Extract the mapped reads with MAPQ>=10 and mapped mate   #
############################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule mapq10:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.bam"
    output:
        temp("{sp}/{tiss}/{anim}/mapping/{sample}.mono.bam")
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.mapq10.log"
    threads:
        4
    shell:
        "samtools view -@ {threads} -F 12 -q 10 -bh {input}"
        " > {output} 2> {log}"

##############################################################################################
# Sort the reads by name to keep pairs with both reads mapped with mapq10 on the same chrom  #
##############################################################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule mapq10paired:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.mono.bam"
    output:
        temp("{sp}/{tiss}/{anim}/mapping/{sample}.paired.unsorted.bam")
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.mapq10.paired.log"
    threads:
        4
    shell:
        "samtools sort -n -@ {threads} -m 4G {input} -o - |"
        " samtools view -h - | awk -f scripts/keep.paired.reads.awk |"
        " samtools view -Sbh - > {output} 2> {log}"

#################################################
# Sort the reads after the mapq10 paired filter #
#################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule bamq10sort:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.paired.unsorted.bam"
    output:
        "{sp}/{tiss}/{anim}/mapping/{sample}.q10.bam"
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.bamq10.sort.log"
    threads:
        4
    shell:
        "samtools sort -@ {threads} -m 4G {input} -O bam -o {output} 2> {log}"

################################################
# Make index of the resulting final bam file   #
################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule bamq10index:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.q10.bam"
    output:
        "{sp}/{tiss}/{anim}/mapping/{sample}.q10.bam.bai"
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.bamq10.index.log"
    shell:
        "samtools index {input} 2> {log}"

###################################################
# Make flagstat of the resulting final bam file   #
###################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule bamq10flag:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.q10.bam"
    output:
        "{sp}/{tiss}/{anim}/mapping/{sample}.q10.bam.flagstat"
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.bamq10.flag.log"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"

###################################################
# Make idxstats of the resulting final bam file   #
###################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule bamq10_idxstats:
    input:
        "{sp}/{tiss}/{anim}/mapping/{sample}.q10.bam"
    output:
        "{sp}/{tiss}/{anim}/mapping/{sample}.q10.bam.idxstats"
    log:
        "logs/mapping/{sp}/{tiss}/{anim}/{sample}.bamq10.idxstats.log"
    shell:
        "samtools idxstats {input} > {output} 2> {log}"
