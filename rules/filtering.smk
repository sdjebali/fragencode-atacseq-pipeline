##############################
# Remove mitochondrial reads #
##############################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule remove_mito:
    input:
        bam="{sp}/{tiss}/{anim}/mapping/{sample}.q10.bam",
        bai="{sp}/{tiss}/{anim}/mapping/{sample}.q10.bam.bai"
    output:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    log:
        "logs/filtering/{sp}/{tiss}/{anim}/{sample}.noMT.log"
    params:
        mito=config["mito"]
    shell:
        "samtools idxstats {input.bam} | awk '{{print $1}}' | grep -v {params.mito} |"
        " xargs samtools view -bh {input} > {output} 2> {log}"


###############################################################
# Making index of bam file after mitochondrial read removal   #
###############################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule nomt_index:
    input:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    output:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.bam.bai"
    log:
        "logs/filtering/{sp}/{tiss}/{anim}/{sample}.noMTindex.log"
    shell:
        "samtools index {input} 2> {log}"


######################################################################
# Making idxstats of the bam file after mitochondrial read removal   #
######################################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule nomt_idxstats:
    input:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    output:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.bam.idxstats"
    log:
        "logs/filtering/{sp}/{tiss}/{anim}/{sample}.noMTidxstats.log"
    shell:
        "samtools idxstats {input} > {output} 2> {log}"

######################################################################
# Making flagstat of the bam file after mitochondrial read removal   #
######################################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule nomt_flag:
    input:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    output:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.bam.flagstat"
    log:
        "logs/filtering/{sp}/{tiss}/{anim}/{sample}.noMTflagstat.log"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"


##########################
# Remove PCR duplicates  #
##########################
# Use Java8
# module load system/Java8
# Use picard 2.20.4 see config.yaml
rule remove_pcrdup:
    input:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    output:
        bam="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam",
        metrics="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.metrics.txt"
    log:
        "logs/filtering/{sp}/{tiss}/{anim}/{sample}.noMTremovePCRdup.log"
    params:
        picard=config["picard"]
    shell:
        "java -jar -Xmx4g {params.picard} MarkDuplicates I={input} O={output.bam}"
        " M={output.metrics} REMOVE_DUPLICATES=true 2> {log}"


#######################################################
# Make index of bam file after PCR duplicate removal  #
#######################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule nomtdup_index:
    input:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam"
    output:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.bai"
    log:
        "logs/filtering/{sp}/{tiss}/{anim}/{sample}.nomtdupindex.log"
    shell:
        "samtools index {input} 2> {log}"

##############################################################
# Make idxstats of the bam file after PCR duplicate removal  #
##############################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule nomtdup_idxstats:
    input:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam"
    output:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.idxstats"
    log:
        "logs/filtering/{sp}/{tiss}/{anim}/{sample}.nomtdup_idxstats.log"
    shell:
        "samtools idxstats {input} > {output} 2> {log}"

##############################################################
# Make flagstat of the bam file after PCR duplicate removal  #
##############################################################
# Use samtools version 1.9
# module load bioinfo/samtools-1.9
rule nomtdup_flagstat:
    input:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam"
    output:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.flagstat"
    log:
        "logs/filtering/{sp}/{tiss}/{anim}/{sample}.nomtdup_flag.log"
    shell:
        "samtools flagstat {input} > $nodup.flagstat 2> {log}"

############################################################################
# Colleting multiple metrix from the bam file after PCR duplicate removal  #
############################################################################
# Use Java8
# module load system/Java8
rule collect_metrics:
    input:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam"
    output:
        insize1="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.insert_size_metrics",
        insize2="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.insert_size_histogram.pdf",
        align="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.alignment_summary_metrics",
        qual1="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.quality_distribution_metrics",
        qual2="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.quality_distribution.pdf",
        qual3="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.quality_by_cycle_metrics",
        qual4="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.quality_by_cycle.pdf"
    log:
        "logs/filtering/{sp}/{tiss}/{anim}/{sample}.nomtdup_collectmetrics.log"
    params:
    # note that the following is not generic since we want to apply it to any species not just to sus_scrofa
        picard=config["picard"],
        genome=config["genome"]["sus_scrofa"],
        prefix="{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM"
    shell:
        "java -jar -Xmx4g {params.picard} CollectMultipleMetrics"
        " I={input} O={params.prefix} R={params.genome}"
