#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationMixturefewersignaturespairedSkinMelanomacutaneousPCAWG_200_200_80_4_-12_NA_NA_NA_dataset14.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturefewersignaturespairedSkinMelanomacutaneousPCAWG_200_200_80_4_-12_fullREDMsinglelambda_NA_NA_NA_dataset14.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationMixturefewersignaturespairedSkinMelanomacutaneousPCAWG", "n": "200", "nlambda": "200", "lmbda": "80", "d": "4", "beta_intensity": "-12", "model": "fullREDMsinglelambda", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "14"}, "params": {"model": "fullREDMsinglelambda", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationMixturefewersignaturespairedSkinMelanomacutaneousPCAWG_200_200_80_4_-12_ModelfullREDMsinglelambda_NA_NA_NA_dataset14.log"], "threads": 1, "resources": {}, "jobid": 3795, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturefewersignaturespairedSkinMelanomacutaneousPCAWG_200_200_80_4_-12_fullREDMsinglelambda_NA_NA_NA_dataset14.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.x53otxz9 ../data/assessing_models_simulation/datasets/multiple_GenerationMixturefewersignaturespairedSkinMelanomacutaneousPCAWG_200_200_80_4_-12_NA_NA_NA_dataset14.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.x53otxz9/3795.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.x53otxz9/3795.jobfailed; exit 1)

