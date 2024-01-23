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

If running in a cluster, you might have to specify additional parameters (see the example below, for a slurm system)

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
### PCAWG samples for which there is no VCF file
These files appear in the metadata and for which we might have the mutccf file, but not their VCF, which is the only file that contains details about the mutation (the mutccf files give information about position and CCF, but not mutation type).

```
../data/restricted/pcawg/pcawg_restricted_snv_counts/f8467ec8-2d61-ba21-e040-11ac0c483584
../data/restricted/pcawg/pcawg_restricted_snv_counts/f856fa85-fdb8-c0b0-e040-11ac0d480b4e
../data/restricted/pcawg/pcawg_restricted_snv_counts/f8696c79-b165-92a6-e040-11ac0c4804bf
../data/restricted/pcawg/pcawg_restricted_snv_counts/f393bb00-888d-710f-e040-11ac0d484518
../data/restricted/pcawg/pcawg_restricted_snv_counts/fc8130df-1e8f-c879-e040-11ac0d485df4
../data/restricted/pcawg/pcawg_restricted_snv_counts/fc8130e0-0dcf-b558-e040-11ac0c483285
../data/restricted/pcawg/pcawg_restricted_snv_counts/2ce48f01-2f61-49d9-a56a-7438bf4a37d7
../data/restricted/pcawg/pcawg_restricted_snv_counts/841eb82b-347d-4d7f-805f-3f3701a2983d
../data/restricted/pcawg/pcawg_restricted_snv_counts/4e361622-f9a8-4e9b-a89e-19bafebe1d6a
../data/restricted/pcawg/pcawg_restricted_snv_counts/5ab6a1d3-76f8-45d4-a430-d9831daa9ec4
../data/restricted/pcawg/pcawg_restricted_snv_counts/86ae34f9-e16a-4593-8e55-b1296782bc1f
../data/restricted/pcawg/pcawg_restricted_snv_counts/924bcc4a-c982-43bf-8bbb-641dc983d65e
../data/restricted/pcawg/pcawg_restricted_snv_counts/00493087-9d9d-40ca-86d5-936f1b951c93
../data/restricted/pcawg/pcawg_restricted_snv_counts/2182ce2c-5941-4b65-9419-fc7966d5e6d5
../data/restricted/pcawg/pcawg_restricted_snv_counts/303abbe5-4155-4a0d-bc3b-f8995261ca52
../data/restricted/pcawg/pcawg_restricted_snv_counts/31e63f89-a6a9-40fb-823d-f41587bd73d8
../data/restricted/pcawg/pcawg_restricted_snv_counts/711c8a16-3cf8-42d8-b29e-fd1e9ef1c82b
../data/restricted/pcawg/pcawg_restricted_snv_counts/8853cbee-7931-49a6-b063-a806943a10ad
../data/restricted/pcawg/pcawg_restricted_snv_counts/d3aff5d3-23c0-43ae-9c01-8ddd776b530b
../data/restricted/pcawg/pcawg_restricted_snv_counts/ec43c4b5-fb72-4a4a-af03-10c2d05ff159
../data/restricted/pcawg/pcawg_restricted_snv_counts/f4b9d98f-7b76-4eaa-9595-10b0973d5ff7
../data/restricted/pcawg/pcawg_restricted_snv_counts/526b3796-2cbd-4eec-8273-064b41456279
../data/restricted/pcawg/pcawg_restricted_snv_counts/d25a4c65-9cb4-4611-909e-e68f93408d84
../data/restricted/pcawg/pcawg_restricted_snv_counts/322f0b01-2118-4dbe-aba1-3875a54ee71b
../data/restricted/pcawg/pcawg_restricted_snv_counts/9d29543e-8601-4fd0-8e76-3df3de465cab
../data/restricted/pcawg/pcawg_restricted_snv_counts/005794f1-5a87-45b5-9811-83ddf6924568
../data/restricted/pcawg/pcawg_restricted_snv_counts/0040b1b6-b07a-4b6e-90ef-133523eaf412
../data/restricted/pcawg/pcawg_restricted_snv_counts/84a0bc36-9f29-4b23-aee1-bf5ff71f697b
../data/restricted/pcawg/pcawg_restricted_snv_counts/cb573c96-f6d4-4897-8919-9827f623b6a7
../data/restricted/pcawg/pcawg_restricted_snv_counts/ef78f09c-c622-11e3-bf01-24c6515278c0
../data/restricted/pcawg/pcawg_restricted_snv_counts/f5a97315-1906-4774-980e-0879c6ad368e
../data/restricted/pcawg/pcawg_restricted_snv_counts/692dfa4f-45e5-4183-b5da-6650a1fbcabd
../data/restricted/pcawg/pcawg_restricted_snv_counts/2f324d8b-ec7b-4d3b-9d64-65f9fc6630a2
../data/restricted/pcawg/pcawg_restricted_snv_counts/0ead45d8-d785-4404-8319-2ef951e02e03
../data/restricted/pcawg/pcawg_restricted_snv_counts/37522f18-77b2-4414-8df8-3c2c8048adba
../data/restricted/pcawg/pcawg_restricted_snv_counts/536dedba-46c4-4a21-b112-13c030b13069
../data/restricted/pcawg/pcawg_restricted_snv_counts/541f91bd-7e9d-4348-9e78-45b948d8967e
../data/restricted/pcawg/pcawg_restricted_snv_counts/92dc0e0c-842f-40de-9c39-486b491ea80a
../data/restricted/pcawg/pcawg_restricted_snv_counts/a6ebe0c0-8aab-4b9f-8328-4b795895a77d
../data/restricted/pcawg/pcawg_restricted_snv_counts/bc395326-1656-4ef2-bb19-0cb29194b91c
../data/restricted/pcawg/pcawg_restricted_snv_counts/dd4fdb6a-8067-4b64-ab74-bbb0fec34ca9
../data/restricted/pcawg/pcawg_restricted_snv_counts/ec77847e-48fd-4ba5-bc3e-3cd1b149b552
../data/restricted/pcawg/pcawg_restricted_snv_counts/4d2204f1-be84-4f58-b7e7-61ae9fbf6d25
../data/restricted/pcawg/pcawg_restricted_snv_counts/94b5dc5a-701a-45e3-8f63-8231031a055a
../data/restricted/pcawg/pcawg_restricted_snv_counts/c307688c-b1fa-47f6-a9e2-1ea41f7645b6
../data/restricted/pcawg/pcawg_restricted_snv_counts/005e85a3-3571-462d-8dc9-2babfc7ace21
../data/restricted/pcawg/pcawg_restricted_snv_counts/07531318-87e8-4db8-aa61-9b93597d063b
../data/restricted/pcawg/pcawg_restricted_snv_counts/33de44a2-bec1-402d-872c-d78c1f2d52b3
../data/restricted/pcawg/pcawg_restricted_snv_counts/4d11d7da-1204-437e-87b1-e8337a67c9a8
../data/restricted/pcawg/pcawg_restricted_snv_counts/4e596add-a7c5-4617-9649-b4ac6612e39c
../data/restricted/pcawg/pcawg_restricted_snv_counts/db198301-6c69-4d56-88d1-c650406423dd
../data/restricted/pcawg/pcawg_restricted_snv_counts/f9854144-d92c-46da-ac87-9d1fd7efe67d
../data/restricted/pcawg/pcawg_restricted_snv_counts/f9c52187-2e82-d58a-e040-11ac0d484fc4
../data/restricted/pcawg/pcawg_restricted_snv_counts/f9c52414-385d-8cc7-e040-11ac0d485037
../data/restricted/pcawg/pcawg_restricted_snv_counts/f9c6f4ca-4bb8-26b4-e040-11ac0d485600
../data/restricted/pcawg/pcawg_restricted_snv_counts/f9c70e38-dd99-3fe2-e040-11ac0d4862f2
../data/restricted/pcawg/pcawg_restricted_snv_counts/1d0617e8-2725-4411-b50f-e46ea1d43242
../data/restricted/pcawg/pcawg_restricted_snv_counts/c73f3f82-3091-46dd-b667-e96a1d8c501c
../data/restricted/pcawg/pcawg_restricted_snv_counts/dbd834cb-b14f-4380-9741-f96551268447
../data/restricted/pcawg/pcawg_restricted_snv_counts/dc107863-2c7d-4b19-8afb-666c7798f0da
../data/restricted/pcawg/pcawg_restricted_snv_counts/dd9d2e9e-02dc-40fc-842b-c5b9707fca56
../data/restricted/pcawg/pcawg_restricted_snv_counts/ffb4f42b-58e9-40c3-8963-11804f041375
../data/restricted/pcawg/pcawg_restricted_snv_counts/577d5c9e-fbda-41d5-b0b3-cdb733453ea5
../data/restricted/pcawg/pcawg_restricted_snv_counts/638e80c7-9a6e-4a32-a621-fc4168e72343
../data/restricted/pcawg/pcawg_restricted_snv_counts/950486ad-14f8-480a-b079-9cc3cd842090
../data/restricted/pcawg/pcawg_restricted_snv_counts/9d2671b9-bd30-4e3c-aa74-01e31dd2531e
../data/restricted/pcawg/pcawg_restricted_snv_counts/accfc45b-eae0-4991-a488-e217cdb46655
../data/restricted/pcawg/pcawg_restricted_snv_counts/bcf858fd-cc3b-4fde-ab10-eb96216f4366
../data/restricted/pcawg/pcawg_restricted_snv_counts/c082dc34-457e-40ec-8258-e11e8ed362c2
```

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

multiple_GenerationMixtureSimulationTwoCT_100_50_NA_NA_-3.75_diagREDM_NA_NA_NA

where the information is as follows

{multiple = fixed}_{GenerationMixtureSimulationTwoCT = simulation framework}_{100 = number of simulated patients}_{50 = number of simulated mutations per subsample}_{NA
Lambda}_{NA = d number of mutational signatures}_{-3.75 - beta1intensity,
differential abundance parameter used to simulate beta1}_{diagREDM = model}_{NA
beta0, path to file with beta0 parameters}_{NA = beta1, path to file with beta1 parameters, incompatible with beta1intensity}_{NA = cov}

Names of datasets
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

