#!/bin/sh
# properties = {"type": "single", "rule": "simulate_for_comparison_multiple", "local": false, "input": ["../data/assessing_models_simulation/GenerationCnorm"], "output": ["../data/assessing_models_simulation/datasets/multiple_GenerationCnorm_80_180_9_6_0_betaintercept2_betaslope2_dataset389.RDS"], "wildcards": {"datasetgeneration": "GenerationCnorm", "n": "80", "nlambda": "180", "lmbda": "9", "d": "6", "beta_intensity": "0", "fixed_beta_intercept": "betaintercept2", "fixed_beta_slope": "betaslope2", "itnum": "389"}, "params": {"datasetgeneration": "GenerationCnorm", "n": "80", "nlambda": "180", "lmbda": "9", "d": "6", "beta_intensity": "0", "itnum": "389", "fixed_beta_intercept": "betaintercept2", "fixed_beta_slope": "betaslope2"}, "log": [], "threads": 1, "resources": {}, "jobid": 3780, "cluster": {}}
 cd /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code && \
PATH='/scratcha/fmlab/morril01/software/miniconda3/envs/snakemake-globalDA/bin':$PATH /scratcha/fmlab/morril01/software/miniconda3/envs/snakemake-globalDA/bin/python3.9 \
-m snakemake ../data/assessing_models_simulation/datasets/multiple_GenerationCnorm_80_180_9_6_0_betaintercept2_betaslope2_dataset389.RDS --snakefile /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.4k91h1as ../data/assessing_models_simulation/GenerationCnorm --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules simulate_for_comparison_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.4k91h1as/3780.jobfinished || (touch /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.4k91h1as/3780.jobfailed; exit 1)

