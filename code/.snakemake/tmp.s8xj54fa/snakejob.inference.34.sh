#!/bin/sh
# properties = {"type": "single", "rule": "inference", "local": false, "input": ["../data/roo/Skin-Melanoma.mucosal_nucleotidesubstitution1_ROO.RDS"], "output": ["../data/inference/Skin-Melanoma.mucosal_nucleotidesubstitution1_15000_LNMROO.RData"], "wildcards": {"cancer_type": "Skin-Melanoma.mucosal", "feature_type": "nucleotidesubstitution1", "nits": "15000", "model": "LNM"}, "params": {"cancer_type": "Skin-Melanoma.mucosal", "feature_type": "nucleotidesubstitution1", "model": ["LNM"], "nits": ["15000"]}, "log": ["logs/inference/Skin-Melanoma.mucosal_nucleotidesubstitution1_15000_LNM.log"], "threads": 1, "resources": {}, "jobid": 34, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/inference/Skin-Melanoma.mucosal_nucleotidesubstitution1_15000_LNMROO.RData --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.s8xj54fa ../data/roo/Skin-Melanoma.mucosal_nucleotidesubstitution1_ROO.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.s8xj54fa/34.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.s8xj54fa/34.jobfailed; exit 1)

