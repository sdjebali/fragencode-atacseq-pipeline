##############################
# Remove mitochondrial reads #
##############################
rule remove_mito:
    input:
        bam="{tiss}/{anim}/mapping/{sample}.q10.bam",
        bai="{tiss}/{anim}/mapping/{sample}.q10.bam.bai"
    output:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    log:
        "logs/filtering/{tiss}/{anim}/{sample}.noMT.log"
    params:
        mito=config["mito"]
    shell:
        "samtools idxstats {input.bam} | awk '{{print $1}}' | grep -v {params.mito} |"
        " xargs samtools view -bh {input} > {output} 2> {log}"


###############################################################
# Making index of bam file after mitochondrial read removal   #
###############################################################
rule nomt_index:
    input:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    output:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.bam.bai"
    log:
        "logs/filtering/{tiss}/{anim}/{sample}.noMTindex.log"
    shell:
        "samtools index {input} 2> {log}"


######################################################################
# Making idxstats of the bam file after mitochondrial read removal   #
######################################################################
rule nomt_idxstats:
    input:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    output:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.bam.idxstats"
    log:
        "logs/filtering/{tiss}/{anim}/{sample}.noMTidxstats.log"
    shell:
        "samtools idxstats {input} > {output} 2> {log}"

######################################################################
# Making flagstat of the bam file after mitochondrial read removal   #
######################################################################
rule nomt_flag:
    input:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    output:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.bam.flagstat"
    log:
        "logs/filtering/{tiss}/{anim}/{sample}.noMTflagstat.log"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"


##########################
# Remove PCR duplicates  #
##########################
rule remove_pcrdup:
    input:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.bam"
    output:
        bam="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam",
        metrics="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.metrics.txt"
    log:
        "logs/filtering/{tiss}/{anim}/{sample}.noMTremovePCRdup.log"
    params:
        picard=config["picard"]
    shell:
        "java -jar -Xmx4g {params.picard} MarkDuplicates I={input} O={output.bam}"
        " M={output.metrics} REMOVE_DUPLICATES=true 2> {log}"


#######################################################
# Make index of bam file after PCR duplicate removal  #
#######################################################
rule nomtdup_index:
    input:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam"
    output:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.bai"
    log:
        "logs/filtering/{tiss}/{anim}/{sample}.nomtdupindex.log"
    shell:
        "samtools index {input} 2> {log}"

##############################################################
# Make idxstats of the bam file after PCR duplicate removal  #
##############################################################
rule nomtdup_idxstats:
    input:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam"
    output:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.idxstats"
    log:
        "logs/filtering/{tiss}/{anim}/{sample}.nomtdup_idxstats.log"
    shell:
        "samtools idxstats {input} > {output} 2> {log}"

##############################################################
# Make flagstat of the bam file after PCR duplicate removal  #
##############################################################
rule nomtdup_flagstat:
    input:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam"
    output:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.flagstat"
    log:
        "logs/filtering/{tiss}/{anim}/{sample}.nomtdup_flag.log"
    shell:
        "samtools flagstat {input} > $nodup.flagstat 2> {log}"

############################################################################
# Colleting multiple metrix from the bam file after PCR duplicate removal  #
############################################################################
rule collect_metrics:
    input:
        "{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam"
    output:
        insize1="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.insert_size_metrics",
        insize2="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.insert_size_histogram.pdf",
        align="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.alignment_summary_metrics",
        qual1="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.quality_distribution_metrics",
        qual2="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.quality_distribution.pdf",
        qual3="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.quality_by_cycle_metrics",
        qual4="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM.quality_by_cycle.pdf"
    log:
        "logs/filtering/{tiss}/{anim}/{sample}.nomtdup_collectmetrics.log"
    params:
    # note that the following is not generic since we want to apply it to any species not just to sus_scrofa
        picard=config["picard"],
        genome=config["genome"],
        prefix="{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam.CMM"
    shell:
        "java -jar -Xmx4g {params.picard} CollectMultipleMetrics"
        " I={input} O={params.prefix} R={params.genome}"
