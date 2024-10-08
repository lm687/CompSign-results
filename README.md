## Introduction

In this repository we apply the methods from the package `CompSign` to answer questions on the differential abundance of mutational signatures between clonal and subclonal samples in the PCAWG cohort.

We refer the reader to the [CompSign](https://github.com/lm687/CompSign) github repository for

- example input data
- minimal examples showing how to estimate the parameters
- different variants of the model

all of this information can be accessed in the Vignettes:

```
browseVignettes("CompSign")
```

## Ignored folders
Some elements inside the data folder `data/` is restricted and therefore not available. Specifically, count PCAWG data are indentifiable and with restricted access.

The summarised data which is not identifiable and available in this github repository is, however, enough to run the code:

- `data/roo/` (see [here](https://github.com/lm687/CompSign-results/tree/main/data/roo)) includes the exposures for each of the PCAWG datasets

## Environment
The results from CompSign applied to PCAWG can be reproduced in a conda environment.

To set up the environment, run

```
cd code/
conda env create -f environment.yaml
```

alternatively, to use additional functionality such as being able to run competing models, a more complex environment is found in `environment-extended.yaml`.

To enter the environment, type

```
source activate snakemake-globalDA
...
conda deactivate
```

## Running in cluster

If running in a cluster, you might have to specify additional parameters (see the example below, for a slurm system). Note that for Snakemake >=8 the syntax has changed substantially. Here Snakemake < 8 is used.

```
source activate snakemake-globalDA
snakemake --cluster "sbatch -t 0:10:00 --cores 1" --jobs 40 --printshellcmds
conda deactivate
```

<!-- ## Creating Snakemake's config file
The snakemake pipeline needs an input file, `config_PCAWG.yaml`. Y -->

<!-- The file `config_PCAWG.yaml` is created by running
```
sh make_config.sh
```

`make_config` contains the arguments (i.e. parameters, for the most part) for the different Simulation Generations. Note using `bash make_config.sh` will throw an error - use sh instead. -->

## Other considerations
Check `text/` for more information.

## Availability of simulated data and of inference results

Large files have been deposited in https://zenodo.org/records/10546525.

Inference results
- TMB_PCAWG: results from the CompSign models applied to PCAWG cohorts
- TMB: results from the CompSign models applied to simulated datasets in which matrices of exposures have been simulated, or in which first mutation substitutions have been generated and then exposures have been extracted from ground truth signatures
- TMB_with_sAgextrast: results from the CompSign models applied to simulated datasets in which first mutation substitutions have been generated and then exposures have been extracted from signatures selected using a variety of active signatures
- TCSM: results from the TCSM model applied to simulated datasets in which mutation substitutions have been generated
- HiLDA: results from the HiLDA model applied to simulated datasets in which mutation substitutions have been generated

Datasets
- Datasets_simulation: input datasets for the models above (all except for MB_PCANG)
- Datasets_PCAWG: input datasets of observed exposures in the PCAWG cohorts (input data for TMB PCAWG)
- Datasets _simulation_paraneters: paraneters used for the synthetic datasets (e.g. betal parameters)
G bay

An example CompSign run name is:
```
  multiple_GenerationMixtureSimulationTwoCT_100_50_NA_NA_-3.75_diagREDM_NA_NA_NA
```

where the information is as follows
```
  {multiple = fixed}_{GenerationMixtureSimulationTwoCT = simulation framework}_{100 = number of simulated patients}_{50 = number of simulated mutations per subsample}_{NA
  Lambda}_{NA = d number of mutational signatures}_{-3.75 - beta1intensity,
  differential abundance parameter used to simulate beta1}_{diagREDM = model}_{NA
  beta0, path to file with beta0 parameters}_{NA = beta1, path to file with beta1 parameters, incompatible with beta1intensity}_{NA = cov}
```

Names of datasets
```
  A1    GenerationJnorm_200_180_2_5_0_betaintercept1d4_betaslope1d4_covmat1d4
  A2    GenerationJnorm_200_180_20_5_0_betaintercept1d4_betaslope1d4_covmat1d4
  A3    GenerationJnorm_200_180_100_5_0_betaintercept1d4_betaslope1d4_covmat1d4
  
  B1    GenerationJnorm_200_3401_18_6_0_betainterceptPCAWG2_betaslopePCAWG2_covmatPCAWG2
  B2    GenerationJnorm_200_14072_80_6_0_betainterceptPCAWG4_betaslopePCAWG4_covmatPCAWG4
  B3    GenerationJnorm_200_14072_87_6_0_betainterceptPCAWG4_betaslopePCAWG4_covmatFULLPCAWG4
  B4    GenerationJnorm_nPCAWG6_{nlambdaPCAWG6,lownlambdaPCAWG6,low2nlambdaPCAWG6}_lambdaPCAWG6_dPCAWG6_0_betainterceptPCAWG6_betaslopePCAWG6_covmatFULLPCAWG6
  
  C1    GenerationMixtureSimulationTwoCT
  C2    GenerationMixtureSimulationv4
  C3    GenerationMixtureSimulation
  C3    GenerationMixtureSimulationv7
```
``
