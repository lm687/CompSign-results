#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationMGnorm_200_180_100_6_0_betaintercept3_betaslope3_sdRE1_dataset189.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMGnorm_200_180_100_6_0_fullREDM_betaintercept3_betaslope3_sdRE1_dataset189.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationMGnorm", "n": "200", "nlambda": "180", "lmbda": "100", "d": "6", "beta_intensity": "0", "model": "fullREDM", "fixed_beta_intercept": "betaintercept3", "fixed_beta_slope": "betaslope3", "sdRE_input": "sdRE1", "itnum": "189"}, "params": {"model": "fullREDM", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationMGnorm_200_180_100_6_0_ModelfullREDM_betaintercept3_betaslope3_sdRE1_dataset189.log"], "threads": 1, "resources": {}, "jobid": 190, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMGnorm_200_180_100_6_0_fullREDM_betaintercept3_betaslope3_sdRE1_dataset189.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.6d_v42et ../data/assessing_models_simulation/datasets/multiple_GenerationMGnorm_200_180_100_6_0_betaintercept3_betaslope3_sdRE1_dataset189.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.6d_v42et/190.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.6d_v42et/190.jobfailed; exit 1)

