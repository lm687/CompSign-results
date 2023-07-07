#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationMixturefewersignaturespairedUterusAdenoCAPCAWG_100_200_80_4_-999_NA_NA_NA_dataset2.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturefewersignaturespairedUterusAdenoCAPCAWG_100_200_80_4_-999_diagREDMsinglelambda_NA_NA_NA_dataset2.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationMixturefewersignaturespairedUterusAdenoCAPCAWG", "n": "100", "nlambda": "200", "lmbda": "80", "d": "4", "beta_intensity": "-999", "model": "diagREDMsinglelambda", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "2"}, "params": {"model": "diagREDMsinglelambda", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationMixturefewersignaturespairedUterusAdenoCAPCAWG_100_200_80_4_-999_ModeldiagREDMsinglelambda_NA_NA_NA_dataset2.log"], "threads": 1, "resources": {}, "jobid": 9483, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturefewersignaturespairedUterusAdenoCAPCAWG_100_200_80_4_-999_diagREDMsinglelambda_NA_NA_NA_dataset2.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.x53otxz9 ../data/assessing_models_simulation/datasets/multiple_GenerationMixturefewersignaturespairedUterusAdenoCAPCAWG_100_200_80_4_-999_NA_NA_NA_dataset2.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.x53otxz9/9483.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.x53otxz9/9483.jobfailed; exit 1)

