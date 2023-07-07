#!/bin/sh
# properties = {"type": "single", "rule": "simulate_for_comparison_multiple", "local": false, "input": ["../data/assessing_models_simulation/GenerationJnorm3"], "output": ["../data/assessing_models_simulation/datasets/multiple_GenerationJnorm3_100_100_80_7_0.6_NA_NA_NA_dataset3.RDS"], "wildcards": {"datasetgeneration": "GenerationJnorm3", "n": "100", "nlambda": "100", "lmbda": "80", "d": "7", "beta_intensity": "0.6", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "3"}, "params": {"datasetgeneration": "GenerationJnorm3", "n": "100", "nlambda": "100", "lmbda": "80", "d": "7", "beta_intensity": "0.6", "itnum": "3", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA"}, "log": [], "threads": 1, "resources": {}, "jobid": 43200, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/datasets/multiple_GenerationJnorm3_100_100_80_7_0.6_NA_NA_NA_dataset3.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.zx8fjxen ../data/assessing_models_simulation/GenerationJnorm3 --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules simulate_for_comparison_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.zx8fjxen/43200.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.zx8fjxen/43200.jobfailed; exit 1)

