#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationJnorm2_20_100_10_6_0_NA_NA_NA_dataset3.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm2_20_100_10_6_0_fullREDMsinglelambda_NA_NA_NA_dataset3.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationJnorm2", "n": "20", "nlambda": "100", "lmbda": "10", "d": "6", "beta_intensity": "0", "model": "fullREDMsinglelambda", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "3"}, "params": {"model": "fullREDMsinglelambda", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationJnorm2_20_100_10_6_0_ModelfullREDMsinglelambda_NA_NA_NA_dataset3.log"], "threads": 1, "resources": {}, "jobid": 18232, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm2_20_100_10_6_0_fullREDMsinglelambda_NA_NA_NA_dataset3.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.zx8fjxen ../data/assessing_models_simulation/datasets/multiple_GenerationJnorm2_20_100_10_6_0_NA_NA_NA_dataset3.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.zx8fjxen/18232.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.zx8fjxen/18232.jobfailed; exit 1)

