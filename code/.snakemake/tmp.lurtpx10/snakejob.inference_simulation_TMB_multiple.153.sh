#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationMGnorm_80_180_100_6_0_betaintercept3_betaslope3_sdRE1_dataset152.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMGnorm_80_180_100_6_0_diagREDM_betaintercept3_betaslope3_sdRE1_dataset152.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationMGnorm", "n": "80", "nlambda": "180", "lmbda": "100", "d": "6", "beta_intensity": "0", "model": "diagREDM", "fixed_beta_intercept": "betaintercept3", "fixed_beta_slope": "betaslope3", "sdRE_input": "sdRE1", "itnum": "152"}, "params": {"model": "diagREDM", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationMGnorm_80_180_100_6_0_ModeldiagREDM_betaintercept3_betaslope3_sdRE1_dataset152.log"], "threads": 1, "resources": {}, "jobid": 153, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMGnorm_80_180_100_6_0_diagREDM_betaintercept3_betaslope3_sdRE1_dataset152.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.lurtpx10 ../data/assessing_models_simulation/datasets/multiple_GenerationMGnorm_80_180_100_6_0_betaintercept3_betaslope3_sdRE1_dataset152.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.lurtpx10/153.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.lurtpx10/153.jobfailed; exit 1)

