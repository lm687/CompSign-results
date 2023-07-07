#!/bin/sh
# properties = {"type": "single", "rule": "simulate_for_comparison_multiple", "local": false, "input": ["../data/assessing_models_simulation/GenerationJnorm"], "output": ["../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_200_180_100_5_0_betaintercept1d4_betaslopecov1d4_cov1d4_dataset507.RDS"], "wildcards": {"datasetgeneration": "GenerationJnorm", "n": "200", "nlambda": "180", "lmbda": "100", "d": "5", "beta_intensity": "0", "fixed_beta_intercept": "betaintercept1d4", "fixed_beta_slope": "betaslopecov1d4", "sdRE_input": "cov1d4", "itnum": "507"}, "params": {"datasetgeneration": "GenerationJnorm", "n": "200", "nlambda": "180", "lmbda": "100", "d": "5", "beta_intensity": "0", "itnum": "507", "fixed_beta_intercept": "betaintercept1d4", "fixed_beta_slope": "betaslopecov1d4", "sdRE_input": "cov1d4"}, "log": [], "threads": 1, "resources": {}, "jobid": 1508, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_200_180_100_5_0_betaintercept1d4_betaslopecov1d4_cov1d4_dataset507.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.8n_d2xv1 ../data/assessing_models_simulation/GenerationJnorm --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules simulate_for_comparison_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.8n_d2xv1/1508.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.8n_d2xv1/1508.jobfailed; exit 1)

