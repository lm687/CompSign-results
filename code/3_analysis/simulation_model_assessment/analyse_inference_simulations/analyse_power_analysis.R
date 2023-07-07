## Power analysis for the models

rm(list = ls())

debugging = T
if(debugging){
  ## debugging
  setwd("/Users/morril01/Documents/PhD/GlobalDA/code/")
}

library(optparse)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(TMB)
library(wesanderson)
source("1_create_ROO/roo_functions.R")
source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("2_inference_TMB/mm_multinomial/helper_functions.R")
source("2_inference_TMB/helper_TMB.R")


debugging = T
if(debugging){
  opt = list()
  opt$dataset_generation = 'GenerationJnorm'
  opt$input_list = list.files("../data/assessing_models_simulation/inference_results/TMB/nlminb/", full.names = T)
  opt$input_list = opt$input_list[grep(paste0(opt$dataset_generation, '_') , opt$input_list)]
  opt$input_list = opt$input_list[grep(opt$model , opt$input_list)]
  opt$input_list = opt$input_list[grepl("multiple_", opt$input_list)]
  opt$output_string = paste0(opt$dataset_generation , '_', opt$model, '_manual')
  # opt$output_folder_name = paste0("../results/results_TMB/simulated_datasets/mixed_effects_models_multiple/",
  #                                 opt$dataset_generation , "/", opt$dataset_generation , "_", opt$model, "_manual/")
  
  opt$multiple_runs = T
}else{
  opt$input_list = strsplit(opt$input_list, " ")[[1]]
}

if(opt$multiple_runs){
  flder_save= "../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries_multiple/"
}else{
  flder_save= "../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries/"
}

opt$input_list = sort(opt$input_list)
print(opt$input_list)

# folder_output = opt$output_folder_name
# system(paste0("mkdir -p ", folder_output))

datasets_files = list.files("../data/assessing_models_simulation/datasets/", full.names = TRUE)
datasets_files = datasets_files[grep(pattern = opt$dataset_generation, datasets_files)]

if(opt$multiple_runs){
  ## if using multiple, only select multiple
  datasets_files = datasets_files[grep('multiple_', datasets_files)]
}else{
  ## if using single dataset+run, remove all multiple datasets
  datasets_files = datasets_files[-grep('multiple_', datasets_files)]
}

if(opt$dataset_generation == 'GenerationCnorm'){
  datasets_files = datasets_files[-grep(pattern = 'GenerationCnormsimpler', datasets_files)]
}

mdels_input <- sapply(opt$input_list, function(i) strsplit(basename(i), '_')[[1]][8])
table(mdels_input)
fles_split_by_model <- split(opt$input_list, mdels_input)
names(fles_split_by_model) <- unique(mdels_input)

datasets = lapply(datasets_files, readRDS)
names(datasets) = gsub(".RDS", "", basename(datasets_files))
DA_bool = ( sapply(datasets, function(i) i$beta_gamma_shape) > 0 )


datasets_files[[1]]


## our effect is going to be the perturbation, averaged by signatures, as in the function below
give_totalperturbation_TMBobj_sigaverage

## compute this for each of the datasets, and compare it to the beta_slope that I have used to simulate them
datasets_files

table(unlist(sapply(datasets, `[`, 'beta_gamma_shape')))


total_perts_datasets <- sapply(datasets, function(dataset_it){
  give_totalperturbation_TMBobj_sigaverage(list(Y=do.call('rbind', slot(dataset_it$objects_counts,
                                                                        'count_matrices_all')),
                                                X=dataset_it$X_sim,
                                                z=t(dataset_it$Z_sim)),
                                           addone = T)})
effectsize1 <- sapply(datasets, function(dataset_it){
  compute_effect_size_1(slot(dataset_it$objects_counts,'count_matrices_all')[[1]],
                        slot(dataset_it$objects_counts,'count_matrices_all')[[2]])})


plot(unlist(sapply(datasets, `[`, 'beta_gamma_shape')), total_perts_datasets)

effect_sizes <- data.frame(beta_gamma_shape=unlist(sapply(datasets, `[`, 'beta_gamma_shape')),
                           total_perturbation= total_perts_datasets,
                           effectsize1=effectsize1,
                           d=unlist(sapply(datasets, `[`, 'd')),
                           n=unlist(sapply(datasets, `[`, 'n')),
                           lambda=unlist(sapply(sapply(datasets, `[`, 'lambda'), function(i) i[1])))


##------------------------------------------------------------------------------------##
## not sensitive to small values of beta_gamma_shape

## not sensitive to d (good)
ggplot(effect_sizes,
       aes(x=beta_gamma_shape+0.001, group=interaction(beta_gamma_shape, d), y=total_perturbation, col=factor(d)))+geom_violin()+
  theme_bw()+scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+geom_jitter(alpha=0.5)

## sensitive to the number of samples (not very good!)
ggplot(effect_sizes,
       aes(x=beta_gamma_shape+0.001, group=interaction(beta_gamma_shape, n), y=total_perturbation, col=factor(n)))+
  geom_violin()+ theme_bw()+scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+
  geom_jitter(alpha=0.5)

ggplot(effect_sizes,
       aes(x=beta_gamma_shape+0.001, group=interaction(beta_gamma_shape, lambda), y=total_perturbation, col=factor(lambda)))+
  geom_violin()+ theme_bw()+scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+
  geom_jitter(alpha=0.04)+
  scale_color_manual(values=wes_palette(n=3, name="Darjeeling2"))+
  theme(legend.position = "bottom")
ggsave(paste0("../results/results_TMB/simulated_datasets/effect_sizes/total_perturbation_simulation_", opt$dataset_generation,
              "_lambda.pdf"), width = 3, height = 3)

##------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------##
## bimodal when beta_gamma_shape is low, due to the overdispersion (see below)
ggplot(effect_sizes,
       aes(x=beta_gamma_shape+0.001, group=beta_gamma_shape, y=effectsize1))+geom_violin()+
  theme_bw()+scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+geom_jitter(alpha=0.05)

## not sensitive to d (good)
ggplot(effect_sizes,
       aes(x=beta_gamma_shape+0.001, group=interaction(beta_gamma_shape, d), y=effectsize1, col=factor(d)))+geom_violin()+
  theme_bw()+scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+geom_jitter(alpha=0.5)

## not sensitive to n (good)
ggplot(effect_sizes,
       aes(x=beta_gamma_shape+0.001, group=interaction(beta_gamma_shape, n), y=effectsize1, col=factor(n)))+
  geom_violin()+ theme_bw()+scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+
  geom_jitter(alpha=0.5)

## why the bimodality?

## looking at lambda
ggplot(effect_sizes,
       aes(x=beta_gamma_shape+0.001, group=interaction(beta_gamma_shape, lambda), y=effectsize1, col=factor(lambda)))+
  geom_violin()+ theme_bw()+scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+
  geom_jitter(alpha=0.04)+
  scale_color_manual(values=wes_palette(n=3, name="Darjeeling2"))+
  theme(legend.position = "bottom")
ggsave(paste0("../results/results_TMB/simulated_datasets/effect_sizes/effectsize1_simulation_", opt$dataset_generation,
              "_lambda.pdf"), width = 3, height = 3)
## clearly, it's due to lambda. This is a cool result!

##------------------------------------------------------------------------------------##
