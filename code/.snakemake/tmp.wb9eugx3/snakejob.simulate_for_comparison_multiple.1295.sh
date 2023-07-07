#!/bin/sh
# properties = {"type": "single", "rule": "simulate_for_comparison_multiple", "local": false, "input": ["../data/assessing_models_simulation/GenerationMGnorm"], "output": ["../data/assessing_models_simulation/datasets/multiple_GenerationMGnorm_200_180_100_6_0_betaintercept3_betaslope3_sdRE1_dataset294.RDS"], "wildcards": {"datasetgeneration": "GenerationMGnorm", "n": "200", "nlambda": "180", "lmbda": "100", "d": "6", "beta_intensity": "0", "fixed_beta_intercept": "betaintercept3", "fixed_beta_slope": "betaslope3", "sdRE_input": "sdRE1", "itnum": "294"}, "params": {"datasetgeneration": "GenerationMGnorm", "n": "200", "nlambda": "180", "lmbda": "100", "d": "6", "beta_intensity": "0", "itnum": "294", "fixed_beta_intercept": "betaintercept3", "fixed_beta_slope": "betaslope3", "sdRE_input": "sdRE1"}, "log": [], "threads": 1, "resources": {}, "jobid": 1295, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/datasets/multiple_GenerationMGnorm_200_180_100_6_0_betaintercept3_betaslope3_sdRE1_dataset294.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.wb9eugx3 ../data/assessing_models_simulation/GenerationMGnorm --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules simulate_for_comparison_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.wb9eugx3/1295.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.wb9eugx3/1295.jobfailed; exit 1)

