#!/bin/sh
# properties = {"type": "single", "rule": "inference", "local": false, "input": ["../data/roo/Liver-HCC_signatures_ROO.RDS"], "output": ["../data/inference/Liver-HCC_signatures_15000_LNMROO.RData"], "wildcards": {"cancer_type": "Liver-HCC", "feature_type": "signatures", "nits": "15000", "model": "LNM"}, "params": {"cancer_type": "Liver-HCC", "feature_type": "signatures", "model": ["LNM"], "nits": ["15000"]}, "log": ["logs/inference/Liver-HCC_signatures_15000_LNM.log"], "threads": 1, "resources": {}, "jobid": 21, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/inference/Liver-HCC_signatures_15000_LNMROO.RData --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp._r6_ks6o ../data/roo/Liver-HCC_signatures_ROO.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp._r6_ks6o/21.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp._r6_ks6o/21.jobfailed; exit 1)

