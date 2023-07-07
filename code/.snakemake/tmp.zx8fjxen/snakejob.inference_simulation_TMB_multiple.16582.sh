#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationJnorm2_100_100_10_3_0.0002_NA_NA_NA_dataset1.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm2_100_100_10_3_0.0002_diagREDM_NA_NA_NA_dataset1.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationJnorm2", "n": "100", "nlambda": "100", "lmbda": "10", "d": "3", "beta_intensity": "0.0002", "model": "diagREDM", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "1"}, "params": {"model": "diagREDM", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationJnorm2_100_100_10_3_0.0002_ModeldiagREDM_NA_NA_NA_dataset1.log"], "threads": 1, "resources": {}, "jobid": 16582, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm2_100_100_10_3_0.0002_diagREDM_NA_NA_NA_dataset1.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.zx8fjxen ../data/assessing_models_simulation/datasets/multiple_GenerationJnorm2_100_100_10_3_0.0002_NA_NA_NA_dataset1.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.zx8fjxen/16582.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.zx8fjxen/16582.jobfailed; exit 1)

