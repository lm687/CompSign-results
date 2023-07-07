#!/bin/sh
# properties = {"type": "single", "rule": "inference", "local": false, "input": ["../data/roo/Breast-AdenoCA_nucleotidesubstitution3_ROO.RDS"], "output": ["../data/inference/Breast-AdenoCA_nucleotidesubstitution3_2000_LNMROO.RData"], "wildcards": {"cancer_type": "Breast-AdenoCA", "feature_type": "nucleotidesubstitution3", "nits": "2000", "model": "LNM"}, "params": {"cancer_type": "Breast-AdenoCA", "feature_type": "nucleotidesubstitution3", "model": ["LNM"], "nits": ["2000"]}, "log": ["logs/inference/Breast-AdenoCA_nucleotidesubstitution3_2000_LNM.log"], "threads": 1, "resources": {}, "jobid": 84, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/inference/Breast-AdenoCA_nucleotidesubstitution3_2000_LNMROO.RData --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.53ad2xse ../data/roo/Breast-AdenoCA_nucleotidesubstitution3_ROO.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.53ad2xse/84.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.53ad2xse/84.jobfailed; exit 1)

