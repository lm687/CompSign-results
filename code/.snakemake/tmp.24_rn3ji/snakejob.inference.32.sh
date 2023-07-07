#!/bin/sh
# properties = {"type": "single", "rule": "inference", "local": false, "input": ["../data/roo/Skin-Melanoma.acral_signatures_ROO.RDS"], "output": ["../data/inference/Skin-Melanoma.acral_signatures_2000_LNMROO.RData"], "wildcards": {"cancer_type": "Skin-Melanoma.acral", "feature_type": "signatures", "nits": "2000", "model": "LNM"}, "params": {"cancer_type": "Skin-Melanoma.acral", "feature_type": "signatures", "model": ["LNM"], "nits": ["2000"]}, "log": ["logs/inference/Skin-Melanoma.acral_signatures_2000_LNM.log"], "threads": 1, "resources": {}, "jobid": 32, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/inference/Skin-Melanoma.acral_signatures_2000_LNMROO.RData --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.24_rn3ji ../data/roo/Skin-Melanoma.acral_signatures_ROO.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.24_rn3ji/32.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.24_rn3ji/32.jobfailed; exit 1)

