#!/bin/sh
# properties = {"type": "single", "rule": "simulate_for_comparison_multiple", "local": false, "input": ["../data/assessing_models_simulation/GenerationInoRE"], "output": ["../data/assessing_models_simulation/datasets/multiple_GenerationInoRE_50_100_10_5_0.4_NA_NA_NA_dataset1.RDS"], "wildcards": {"datasetgeneration": "GenerationInoRE", "n": "50", "nlambda": "100", "lmbda": "10", "d": "5", "beta_intensity": "0.4", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "1"}, "params": {"datasetgeneration": "GenerationInoRE", "n": "50", "nlambda": "100", "lmbda": "10", "d": "5", "beta_intensity": "0.4", "itnum": "1", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA"}, "log": [], "threads": 1, "resources": {}, "jobid": 4388, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/datasets/multiple_GenerationInoRE_50_100_10_5_0.4_NA_NA_NA_dataset1.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.p8uy4k5w ../data/assessing_models_simulation/GenerationInoRE --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules simulate_for_comparison_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.p8uy4k5w/4388.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.p8uy4k5w/4388.jobfailed; exit 1)

