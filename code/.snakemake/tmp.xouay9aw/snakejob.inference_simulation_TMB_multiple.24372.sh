#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationJnorm3_20_100_80_6_0.001_NA_NA_NA_dataset3.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm3_20_100_80_6_0.001_fullREM_NA_NA_NA_dataset3.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationJnorm3", "n": "20", "nlambda": "100", "lmbda": "80", "d": "6", "beta_intensity": "0.001", "model": "fullREM", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "3"}, "params": {"model": "fullREM", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationJnorm3_20_100_80_6_0.001_ModelfullREM_NA_NA_NA_dataset3.log"], "threads": 1, "resources": {}, "jobid": 24372, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm3_20_100_80_6_0.001_fullREM_NA_NA_NA_dataset3.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.xouay9aw ../data/assessing_models_simulation/datasets/multiple_GenerationJnorm3_20_100_80_6_0.001_NA_NA_NA_dataset3.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.xouay9aw/24372.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.xouay9aw/24372.jobfailed; exit 1)

