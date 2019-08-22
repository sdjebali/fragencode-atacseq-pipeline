#!/bin/sh
# properties = {"type": "single", "rule": "map", "local": false, "input": ["sus_scrofa/liver/pig1/trimming/pig1_liver/pig1_liver_R1_val_1.fq.gz", "sus_scrofa/liver/pig1/trimming/pig1_liver/pig1_liver_R2_val_2.fq.gz"], "output": ["sus_scrofa/liver/pig1/mapping/pig1_liver.sam"], "wildcards": {"sp": "sus_scrofa", "tiss": "liver", "anim": "pig1", "sample": "pig1_liver"}, "params": {"genome": "genomes[sus_scrofa]"}, "log": ["logs/mapping/sus_scrofa/liver/pig1/pig1_liver.map.log"], "threads": 4, "resources": {}, "jobid": 9, "cluster": {"mem": 8, "name": "map.anim=pig1,sample=pig1_liver,sp=sus_scrofa,tiss=liver", "time": "24:00:00", "output": "logs/cluster/map.anim=pig1,sample=pig1_liver,sp=sus_scrofa,tiss=liver.out", "error": "logs/cluster/map.anim=pig1,sample=pig1_liver,sp=sus_scrofa,tiss=liver.err"}}
cd /work/project/fragencode/tools/atacseq/snakemake/test && \
/home/sdjebali/.conda/envs/atacseq/bin/python3.6 \
-m snakemake sus_scrofa/liver/pig1/mapping/pig1_liver.sam --snakefile /work/project/fragencode/tools/atacseq/snakemake/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/project/fragencode/tools/atacseq/snakemake/test/.snakemake/tmp.q7rkd9u1 sus_scrofa/liver/pig1/trimming/pig1_liver/pig1_liver_R1_val_1.fq.gz sus_scrofa/liver/pig1/trimming/pig1_liver/pig1_liver_R2_val_2.fq.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /work/project/fragencode/tools/atacseq/snakemake/test/config.yaml -p --allowed-rules map --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda 

