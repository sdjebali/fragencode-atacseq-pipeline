

#####################
# Map trimmed reads #
#####################
# Takes several hours with 4 cpus and 25.5G of ram
rule map:
    input:
        read1="{tiss}/{anim}/trimming/{sample}/{sample}_R1_val_1.fq.gz",
        read2="{tiss}/{anim}/trimming/{sample}/{sample}_R2_val_2.fq.gz"
    output:
        "{tiss}/{anim}/mapping/{sample}.sam"
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.map.log"
    threads:
        4
    params:
    # note that the following is not generic since we want to apply it to any species not just to sus_scrofa
        genome=config["genome"]
    shell:
        "bowtie2 -t -p {threads} -X 2000 -x {params.genome}"
        " -1 {input.read1} -2 {input.read2} --met 60 --met-file {log}"
        " -S {output}"

###################
# Convert to bam  #
###################
rule samtobam:
    input:
        "{tiss}/{anim}/mapping/{sample}.sam"
    output:
        temp("{tiss}/{anim}/mapping/{sample}.unsorted.bam")
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.tobam.log"
    threads:
        4
    shell:
        "samtools view {input} -Sbh -@ {threads} -o {output} 2> {log}"

#####################
# Sort the bam file #
#####################
rule bamsort:
    input:
        "{tiss}/{anim}/mapping/{sample}.unsorted.bam"
    output:
        "{tiss}/{anim}/mapping/{sample}.bam"
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.bamsort.log"
    threads:
        4
    shell:
        "samtools sort -@ {threads} -m 4G {input} -O bam -o {output} 2> {log}"


######################
# Index the bam file #
######################
rule bamindex:
    input:
        "{tiss}/{anim}/mapping/{sample}.bam"
    output:
        "{tiss}/{anim}/mapping/{sample}.bam.bai"
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.bamidx.log"
    shell:
        "samtools index {input} 2> {log}"


#########################
# Flagstat the bam file #
#########################
rule bamflag:
    input:
        "{tiss}/{anim}/mapping/{sample}.bam"
    output:
        "{tiss}/{anim}/mapping/{sample}.bam.flagstat"
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.bamflag.log"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"

############################################################
# Extract the mapped reads with MAPQ>=10 and mapped mate   #
############################################################
rule mapq10:
    input:
        "{tiss}/{anim}/mapping/{sample}.bam"
    output:
        temp("{tiss}/{anim}/mapping/{sample}.mono.bam")
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.mapq10.log"
    threads:
        4
    shell:
        "samtools view -@ {threads} -F 12 -q 10 -bh {input}"
        " > {output} 2> {log}"

##############################################################################################
# Sort the reads by name to keep pairs with both reads mapped with mapq10 on the same chrom  #
##############################################################################################
rule mapq10paired:
    input:
        "{tiss}/{anim}/mapping/{sample}.mono.bam"
    output:
        temp("{tiss}/{anim}/mapping/{sample}.paired.unsorted.bam")
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.mapq10.paired.log"
    params:
        path=snakedir
    threads:
        4
    shell:
        "samtools sort -n -@ {threads} -m 4G {input} -o - |"
        " samtools view -h - | awk -f {params.path}/scripts/keep.paired.reads.awk |"
        " samtools view -Sbh - > {output} 2> {log}"

#################################################
# Sort the reads after the mapq10 paired filter #
#################################################
rule bamq10sort:
    input:
        "{tiss}/{anim}/mapping/{sample}.paired.unsorted.bam"
    output:
        "{tiss}/{anim}/mapping/{sample}.q10.bam"
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.bamq10.sort.log"
    threads:
        4
    shell:
        "samtools sort -@ {threads} -m 4G {input} -O bam -o {output} 2> {log}"

################################################
# Make index of the resulting final bam file   #
################################################
rule bamq10index:
    input:
        "{tiss}/{anim}/mapping/{sample}.q10.bam"
    output:
        "{tiss}/{anim}/mapping/{sample}.q10.bam.bai"
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.bamq10.index.log"
    shell:
        "samtools index {input} 2> {log}"

###################################################
# Make flagstat of the resulting final bam file   #
###################################################
rule bamq10flag:
    input:
        "{tiss}/{anim}/mapping/{sample}.q10.bam"
    output:
        "{tiss}/{anim}/mapping/{sample}.q10.bam.flagstat"
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.bamq10.flag.log"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"

###################################################
# Make idxstats of the resulting final bam file   #
###################################################
rule bamq10_idxstats:
    input:
        "{tiss}/{anim}/mapping/{sample}.q10.bam"
    output:
        "{tiss}/{anim}/mapping/{sample}.q10.bam.idxstats"
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.bamq10.idxstats.log"
    shell:
        "samtools idxstats {input} > {output} 2> {log}"
