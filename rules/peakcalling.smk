
############################
# Calling peaks with MACS2 #
############################
# Use macs2 version 2.1.1.20160309 in /tools/python/2.7.2/bin/macs2
rule peakcall:
    input:
        "{sp}/{tiss}/{anim}/filtering/{sample}.q10.noMT.nodup.bam"
    output:
        "{sp}/{tiss}/{anim}/peaks/{sample}.q10.noMT.nodup.macs_peaks.narrowPeak"
    log:
        "logs/peakcalling/{sp}/{tiss}/{anim}/{sample}.peakcall.log"
    params:
        prefix="{sp}/{tiss}/{anim}/peaks/{sample}.q10.noMT.nodup.macs"
    conda:
        "envs/macs2.yaml"
    shell:
        "macs2 callpeak --nomodel -f BAMPE -t {input} -n {params.prefix}"
        " --keep-dup all --verbose 0 2> {log}"
