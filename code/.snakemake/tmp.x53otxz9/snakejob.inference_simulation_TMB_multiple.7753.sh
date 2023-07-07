#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationMixturefewersignaturespairedThyAdenoCAPCAWG_200_200_80_4_-12_NA_NA_NA_dataset12.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturefewersignaturespairedThyAdenoCAPCAWG_200_200_80_4_-12_diagREDM_NA_NA_NA_dataset12.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationMixturefewersignaturespairedThyAdenoCAPCAWG", "n": "200", "nlambda": "200", "lmbda": "80", "d": "4", "beta_intensity": "-12", "model": "diagREDM", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "12"}, "params": {"model": "diagREDM", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationMixturefewersignaturespairedThyAdenoCAPCAWG_200_200_80_4_-12_ModeldiagREDM_NA_NA_NA_dataset12.log"], "threads": 1, "resources": {}, "jobid": 7753, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturefewersignaturespairedThyAdenoCAPCAWG_200_200_80_4_-12_diagREDM_NA_NA_NA_dataset12.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.x53otxz9 ../data/assessing_models_simulation/datasets/multiple_GenerationMixturefewersignaturespairedThyAdenoCAPCAWG_200_200_80_4_-12_NA_NA_NA_dataset12.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.x53otxz9/7753.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.x53otxz9/7753.jobfailed; exit 1)

