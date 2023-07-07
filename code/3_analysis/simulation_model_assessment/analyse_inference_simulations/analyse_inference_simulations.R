rm(list = ls())

debugging = F
if(debugging){
  ## debugging
  # setwd("/Users/morril01/Documents/PhD/GlobalDA/code/")
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
  setwd("../../../")
  getwd()
  try(library(TMB))
  library(TMB, lib.loc = "/mnt/scratcha/fmlab/morril01/software/miniconda3/lib/R/library/")
  opt <- list()
  opt$input_list = '../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_50_200_80_7_0_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_50_200_80_7_0.01_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_50_200_80_7_0.1_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_50_200_80_7_0.6_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_100_200_80_7_0_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_100_200_80_7_0.01_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_100_200_80_7_0.1_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_100_200_80_7_0.6_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_200_200_80_7_0_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_200_200_80_7_0.01_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_200_200_80_7_0.1_fullREM_NA_NA_NA_dataset0.RDS ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturePCAWG_200_200_80_7_0.6_fullREM_NA_NA_NA_dataset0.RDS'
  opt$input_list <- strsplit(opt$input_list, ' ')[[1]]
  opt$output_folder_name = '../results/results_TMB/simulated_datasets/mixed_effects_models_multiple/GenerationMixturePCAWG/GenerationMixturePCAWG_fullREM/'
  opt$output_string = 'GenerationMixturePCAWG_fullREM'
  opt$model = 'fullREM'
  opt$dataset_generation = 'GenerationMixturePCAWG'
  opt$multiple_runs = T
  
}

library(optparse)
# library(ROCR)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(TMB)
# library(xtable)
# library(mvtnorm)
source("1_create_ROO/roo_functions.R")
source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("2_inference_TMB/mm_multinomial/helper_functions.R")
source("2_inference_TMB/helper_TMB.R")

option_list = list(
  make_option(c("--input_list"), type="character", default=NA,
              help="Text with small description of the type of simulation being carried out", metavar="character"),
  make_option(c("--output_folder_name"), type="character", default=NA,
              help="Name of the folder for the output files", metavar="character"),
  make_option(c("--output_string"), type="character", default=NA,
              help="String to be attached to the output files", metavar="character"),
  make_option(c("--model"), type="character", default=NA,
              help="Name of model", metavar="character"),
  make_option(c("--dataset_generation"), type="character", default=NA,
              help="Generation name", metavar="character"),
  make_option(c("--multiple_runs"), type="logical", default=F,
              help="Boolean: are we analysing mutiple runs?", metavar="logical"),
  make_option(c("--read_input_from_file"), type="logical", default=F,
              help="Boolean: are we reading the input files from a file? (set to T if the number of input files is very large)", metavar="logical")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(debugging){
  ## debugging
  opt = list()
  opt$dataset_generation = 'GenerationJnorm'
  opt$dataset_generation = 'generationFnorm'
  opt$model = 'fullREDMsinglelambda'
  opt$model = 'fullREM'
  opt$model = 'diagREDM'
  # opt$model = 'diagREDMsinglelambda'
  opt$input_list = list.files("../data/assessing_models_simulation/inference_results/TMB/nlminb/", full.names = T)
  opt$input_list = opt$input_list[grep(paste0(opt$dataset_generation, '_') , opt$input_list)]
  opt$input_list = opt$input_list[grep(opt$model , opt$input_list)]
  opt$input_list = opt$input_list[grepl("multiple_", opt$input_list)]
  opt$output_string = paste0(opt$dataset_generation , '_', opt$model, '_manual')
  opt$output_folder_name = paste0("../results/results_TMB/simulated_datasets/mixed_effects_models_multiple/",
                                  opt$dataset_generation , "/", opt$dataset_generation , "_", opt$model, "_manual/")
  
  opt$multiple_runs = T
  # opt$input_list = python_like_select(python_like_select(list.files("../data/assessing_models_simulation/inference_results/TMB/", full.names = T),
  #                                                        "GenerationCnorm_"), "fullREDM")
  # opt$output_folder_name= 'GenerationCnorm_fullREDM'
  # opt$output_string = 'GenerationCnorm_fullREDM'
  # opt$model = 'fullREDM'
  # opt$dataset_generation = 'GenerationCnorm'
}else{
  if(opt$read_input_from_file){
    ## large number of input files: read list from file
    opt$input_list = readLines(opt$input_list)
  }else{
    opt$input_list = strsplit(opt$input_list, " ")[[1]]
  }
}

if(opt$multiple_runs){
  flder_save= "../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries_multiple/"
}else{
  flder_save= "../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries/"
}

opt$input_list = sort(opt$input_list)
print(opt$input_list)

folder_output = opt$output_folder_name
system(paste0("mkdir -p ", folder_output))

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
  if(length(grep(pattern = 'GenerationCnormsimpler', datasets_files)) > 0){
    ## if there are any files from GenerationCnormsimpler (if there are none the line below gives an empty list)
    datasets_files = datasets_files[-grep(pattern = 'GenerationCnormsimpler', datasets_files)]
  }
}

cat('Matching datasets and runs\n')
# match
if(opt$multiple_runs){
  datasets_files = datasets_files[match(gsub(paste0("_", opt$model), "", basename(opt$input_list)),
                                        basename(datasets_files))]
}else{
  datasets_files = datasets_files[match(gsub(paste0("_", opt$model, ".RDS"), "", basename(opt$input_list)),
        gsub("_dataset.RDS", "", basename(datasets_files)))]
}


datasets = lapply(datasets_files, readRDS)

if((opt$dataset_generation %in% c("GenerationMixturePCAWG", "GenerationMixturefewersignaturesPCAWG", "GenerationMixturefewersignaturespairedPCAWG")) | grepl('GenerationMixturefewersignaturespaired', opt$dataset_generation) | grepl('signaturespaired', opt$dataset_generation) ){
  cat('Transforming beta gamma shape from logR to probability')
  datasets <- lapply(datasets, function(i){
    if(i$beta_gamma_shape == -999){
      i$beta_gamma_shape = 0
    }else{
      i$beta_gamma_shape = softmax(c(i$beta_gamma_shape, 0))[1]
    }
    i
  })
}else{
  if(grepl('Mixture', opt$dataset_generation)){
    stop('Are you sure you are using probabilities beta_gamma_shape for and not log-ratios?\n')
  }
}
names(datasets) = gsub(".RDS", "", basename(datasets_files))

print(length(datasets))
print(length(opt$input_list))

DA_bool = ( sapply(datasets, function(i) i$beta_gamma_shape) > 0 )

runs = lapply(opt$input_list, readRDS)
names(runs) <- opt$input_list
# sapply(runs, function(i) i$pdHess)

#' ## assess if there was good convergence
# sapply(runs, typeof)

#' Convert betas to NA whenever the convergence was not successful
bool_has_not_converged <- sapply(runs, function(i) try(give_summary_per_sample(i))) != "Good"
runs[bool_has_not_converged] = NA

cat('Adding slope and intercept to failed runs\n')
print(which(sapply(runs, typeof) %in% c("logical", "character")))

for(j in which(sapply(runs, typeof) %in% c("logical", "character"))){
  # runs[[j]] = list(par.fixed=c(beta=rep(NA, 2*(datasets[[j]]$d-1)))) ## *2 for slope and intercept
  runs[[j]] = list(par.fixed=c(beta=rep(NA, 2*(ncol(datasets[[j]]$W)-1)))) ## *2 for slope and intercept
}
pvals = as.numeric(sapply(runs, function(i) try(wald_TMB_wrapper(i, verbatim=FALSE))))

pvals_adj = pvals#*length(pvals_M)

table(DA_bool, M=pvals_adj <= 0.05)

if(opt$dataset_generation == 'GenerationPois'){
  cat('Betas are not in log-ratio space')
  ## here the betas are not in log-ratios
  for(j in 1:length(datasets)){
    datasets[[j]]$effective_d <- datasets[[j]]$d + 1
  }
}else{
  cat('Datasets created with logratios\n')
  
  ## the betas are in log-ratios
  for(j in 1:length(datasets)){
    datasets[[j]]$effective_d <- datasets[[j]]$d
  }
}

if(!is.null(dim(datasets[[1]]$beta))){
  ## if the datasets are created from the DM or M models
  df_beta_recovery = cbind.data.frame(beta_true = unlist(sapply(datasets, function(i) i$beta[2,])),
                                      beta_intercept_true = unlist(sapply(datasets, function(i) i$beta[1,])),
                                      idx = rep(1:length(datasets) , unlist(sapply(datasets, function(i) i$effective_d))-1),
                                      d =  rep(unlist(sapply(datasets, function(i) i$d)), unlist(sapply(datasets, function(i) i$effective_d))-1),
                                      n =  rep(unlist(sapply(datasets, function(i) i$n)), unlist(sapply(datasets, function(i) i$effective_d))-1),
                                      beta_gamma_shape =  rep(unlist(sapply(datasets, function(i) i$beta_gamma_shape)), unlist(sapply(datasets, function(i) i$effective_d))-1),
                                      beta_est = unlist(sapply(runs, function(i) select_slope_2(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE))),
                                      beta_intercept_est = unlist(sapply(runs, function(i) select_intercept(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE))),
                                      beta_stderr = unlist(sapply(runs, give_stderr)),
                                      beta_intercept_stderr = unlist(sapply(runs, give_stderr, only_slopes=F)[c(T,F)]),
                                      pvals_adj=rep(pvals_adj, unlist(sapply(datasets, function(i) i$effective_d))-1),
                                      DA_bool=rep(DA_bool, unlist(sapply(datasets, function(i) i$effective_d))-1),
                                      idx_within_dataset=unlist(sapply(datasets, function(i) 1:(i$effective_d-1))))
  
  df_beta_recovery$bool_zero_true_beta = factor(df_beta_recovery$beta_true == 0, levels=c(TRUE, FALSE))
  df_beta_recovery$converged = sapply(sapply(runs, '[', 'pdHess'), function(i) if(is.null((i))){FALSE}else{i})[df_beta_recovery$idx]
  
  first_entries_runs = sapply(unique(df_beta_recovery$idx), function(i) which(df_beta_recovery$idx == i)[1])
  table(Sim=df_beta_recovery$DA_bool[first_entries_runs], est_M=df_beta_recovery$pvals_adj[first_entries_runs] <= 0.05)
  
}else{
  ## if they are created with perturbation, or similar methods
  cat('Model simulated non-parametrically\n')

  cat('d-1\n')
  print(as.vector(sapply(datasets, function(i) length(rep(NA, i$effective_d-1)))))
  cat('give_stderr\n')
  print(as.vector(sapply(runs, function(i) length(give_stderr(i)))))
  print(names(runs)[(as.vector(sapply(datasets, function(i) length(rep(NA, i$effective_d-1))))) != (as.vector(sapply(runs, function(i) length(give_stderr(i)))))])

  cat('Number of datasets: ', length(datasets), '\n')
  cat('Number of runs: ', length(runs), '\n')
  print(sapply(runs, function(i) length(i$par.fixed)))
  print(names(runs)[sapply(runs, function(i) length(i$par.fixed)) == 6 ])
  print((runs)[sapply(runs, function(i) length(i$par.fixed)) == 6 ])
  
  print(datasets[[1]]$d)
  print(dim(datasets[[1]]$W))
  print(runs[[1]])
  
  
  ### i$d changed to ncol(i$W) because otherwise there could be some problems in datasets that were incorrectly created
  #' in df_beta_recovery, as vector() added to unlist(), because with  opt$dataset_generation
  #' [1] "GenerationMixturePCAWG" I got a df_beta_recovery with matrices inside, i.e. many columns
  df_beta_recovery = cbind.data.frame(beta_true = as.vector(unlist(sapply(datasets, function(i) rep(NA, ncol(i$W)-1)))),
                                      idx = rep(1:length(datasets) , unlist(sapply(datasets, function(i) ncol(i$W)))-1),
                                      d =  rep(unlist(sapply(datasets, function(i) ncol(i$W))), unlist(sapply(datasets, function(i) ncol(i$W)))-1),
                                      n =  rep(unlist(sapply(datasets, function(i) i$n)), unlist(sapply(datasets, function(i) ncol(i$W)))-1),
                                      beta_gamma_shape =  rep(unlist(sapply(datasets, function(i) i$beta_gamma_shape)), unlist(sapply(datasets, function(i) ncol(i$W)))-1),
                                      beta_est = as.vector(unlist(sapply(runs, function(i) select_slope_2(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE)))),
                                      beta_stderr = as.vector(unlist(sapply(runs, give_stderr))),
                                      beta_intercept_est = unlist(sapply(runs, function(i) select_intercept(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE))), ## new 20220210
                                      beta_intercept_stderr = unlist(sapply(runs, give_stderr, only_slopes=F)[c(T,F)]), ## new 20220210
                                      pvals_adj=rep(pvals_adj, unlist(sapply(datasets, function(i) ncol(i$W)))-1),
                                      DA_bool=rep(DA_bool, unlist(sapply(datasets, function(i) ncol(i$W)))-1),
                                      idx_within_dataset=as.vector(unlist(sapply(datasets, function(i) 1:(ncol(i$W)-1)))))
  rownames(df_beta_recovery) <- make.names_allowing_hyphen(names(datasets)[df_beta_recovery$idx]) ## added 12 jan 2022
  df_beta_recovery$bool_zero_true_beta = factor(df_beta_recovery$beta_gamma_shape == 0, levels=c(TRUE, FALSE))
  df_beta_recovery$converged = sapply(sapply(runs, '[', 'pdHess'), function(i) if(is.null((i))){FALSE}else{i})[df_beta_recovery$idx]
  
  first_entries_runs = sapply(unique(df_beta_recovery$idx), function(i) which(df_beta_recovery$idx == i)[1])
  table(Sim=df_beta_recovery$DA_bool[first_entries_runs], est_M=df_beta_recovery$pvals_adj[first_entries_runs] <= 0.05)
  
}

saveRDS(df_beta_recovery,
        paste0(flder_save, opt$output_string, ".RDS"))

#' Remove non-converged
df_beta_recovery = df_beta_recovery[df_beta_recovery$converged,]

#' Remove outliers
df_beta_recovery[which.max(df_beta_recovery$beta_est),]
if(opt$dataset_generation == 'GenerationCnorm' & opt$model == 'fullREM'){
  df_beta_recovery = df_beta_recovery[!(df_beta_recovery$idx == 98),]
}

#' ## zero slopes are not well detected
try({
ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est), col=beta_gamma_shape))+geom_point()+
  geom_abline(intercept = 0, slope = 1)
ggsave(paste0(folder_output, "recovery_betaslope_scatter.pdf"))
})
# ggplot(df_beta_recovery,
#        aes(x=(beta_true), y=(beta_est), col=beta_gamma_shape))+geom_point()+
#   geom_abline(intercept = 0, slope = 1)+facet_wrap(.~idx,scales = "free")
# ggsave(paste0(folder_output, "recovery_betaslope_scatter_facets.pdf"), width = 10, height = 10)


## plotting separately the differentially abundant (right) and non-differentially abundant (left) datasets
try({ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est), col=n))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")
ggsave(paste0(folder_output, "recovery_betaslope_scatter_sepDA.pdf"))
})

#' ## Multinomial
filename_betarecovery1 = paste0(folder_output, "recovery_betaslope2.pdf")
try({
  p <- ggplot(df_beta_recovery,
            aes(x=(beta_true), y=(beta_est), col=(pvals_adj<0.05)))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")
gp <- ggplotGrob(p); facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]; gp$widths[facet.columns] <- gp$widths[facet.columns] * c(1,4); grid::grid.draw(gp)
ggsave(filename_betarecovery1)
})

#' ## Looking at true zeros
#' ## i.e. runs where all beta slopes are zero
filename_betarecovery2 = paste0(folder_output, "betaslopes_nonDA.pdf")
try({
  ggplot(df_beta_recovery[!(df_beta_recovery$DA_bool),],
       aes(x=(beta_true), y=(beta_est), col=(pvals_adj<0.05)))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~(pvals_adj<0.05), scales = "free_x")
ggsave(filename_betarecovery2, height = 3, width = 8)
})

# ggplot(df_beta_recovery[!(df_beta_recovery$DA_bool),],
#        aes(x=idx_within_dataset, y=(beta_est), col=(pvals_adj<0.05)))+
#   geom_point(data = df_beta_recovery[!(df_beta_recovery$DA_bool),], aes(x=idx_within_dataset, y=beta_true), shape=4)+
#   geom_abline(slope = 0, intercept = 0, alpha=0.2)+
#   geom_point()+
#   geom_errorbar(aes(ymin=beta_est-beta_stderr, ymax=beta_est+beta_stderr), width=.2,
#                 position=position_dodge(.9))+theme_bw()+
#   facet_wrap(.~idx, scales='free_x')+theme(legend.position = "bottom")
filename_betarecovery3 = paste0(folder_output, "betaslopes_stderr_nonDA.pdf")
# ggsave(filename_betarecovery3, height = 8, width = 8)

# ggplot(df_beta_recovery[(df_beta_recovery$DA_bool),],
#        aes(x=idx_within_dataset, y=(beta_est), col=(pvals_adj<0.05)))+
#   geom_point(data = df_beta_recovery[df_beta_recovery$DA_bool,], aes(x=idx_within_dataset, y=beta_true), shape=4)+
#   geom_abline(slope = 0, intercept = 0, alpha=0.2)+geom_point()+
#   geom_errorbar(aes(ymin=beta_est-beta_stderr, ymax=beta_est+beta_stderr), width=.2,
#                 position=position_dodge(.9))+
#   facet_wrap(.~interaction(idx, beta_gamma_shape), scales='free_x', nrow=length(unique(df_beta_recovery$beta_gamma_shape)))+
#   theme_bw()+theme(legend.position = "bottom")
filename_betarecovery4 = paste0(folder_output, "betaslopes_stderr_DA.pdf")
# ggsave(filename_betarecovery4, height = 8, width = 8)

print(paste0('Saved files:', c(filename_betarecovery1, filename_betarecovery2, filename_betarecovery3, filename_betarecovery4)))


## How does the FP rate change as several parameters vary?
table(df_beta_recovery$beta_true)

## The correlation between the parameters
## This should be corrected somehow as larger beta_gamma_shape leads to higher betas

try({
pdf(paste0(opt$output_folder_name, "/", opt$output_string, "_facets_distance_estimates.pdf"), width = 7, height = 2)
grid.arrange(ggplot(df_beta_recovery %>% group_by(idx) %>% mutate(dist=dist(rbind(beta_est, beta_true)), .groups="drop"),
       aes(x=log(beta_gamma_shape), y=dist))+geom_point()+geom_smooth()+labs(y='Distance of estimates to true values'),
ggplot(df_beta_recovery %>% group_by(idx) %>% mutate(dist=dist(rbind(beta_est, beta_true)), .groups="drop"),
       aes(x=d, y=dist))+geom_point()+geom_smooth()+labs(y='Distance of estimates to true values'),
ggplot(df_beta_recovery %>% group_by(idx) %>% mutate(dist=dist(rbind(beta_est, beta_true)), .groups="drop"),
       aes(x=n, y=dist))+geom_point()+geom_smooth()+labs(y='Distance of estimates to true values'),
ncol=3)
dev.off()
})

## Compute FP, FN, for each beta
cat('Compute FP, FN, for each beta\n')
sum(is.na(pvals))
df_beta_recovery_accuracy = df_beta_recovery %>% group_by(idx) %>% mutate(acc= paste0(c(unique(pvals_adj<0.05), unique(DA_bool)), collapse='-'))
table(df_beta_recovery_accuracy$acc)

change_da_acc = function(i){
  if(i == 'FALSE-FALSE'){
    'TN'
  }else if(i == 'FALSE-TRUE'){
    'FN'
  }else if(i == 'TRUE-FALSE'){
    'FP'
  }else if(i == 'TRUE-TRUE'){
    'TP'
  }else{
    NA
  }
}


df_beta_recovery_accuracy$acc = sapply(df_beta_recovery_accuracy$acc, change_da_acc)
df_beta_recovery_accuracy$acc
ggplot(df_beta_recovery_accuracy, aes(x=beta_gamma_shape, fill=acc))+geom_bar()

df_beta_recovery_accuracy = df_beta_recovery_accuracy[!is.na(df_beta_recovery_accuracy$acc),]
df_beta_recovery_accuracy_grouped = df_beta_recovery_accuracy %>% group_by(idx) %>% select('idx', 'beta_gamma_shape', 'acc') %>% group_by(beta_gamma_shape) %>% 
  summarise(sensitivity=(sum(acc == 'TP')/(sum(acc == 'TP')+sum(acc == 'FN'))),
            specificity=(sum(acc == 'TN')/(sum(acc == 'TN')+sum(acc == 'FP'))))
df_beta_recovery_accuracy_grouped

cat('Plotting accuracy\n')
pdf(paste0(opt$output_folder_name, "/", opt$output_string, "_specificity_sensitivity.pdf"), width = 7, height = 2)
grid.arrange(ggplot(df_beta_recovery_accuracy_grouped,
       aes(x=beta_gamma_shape, y=sensitivity, group=beta_gamma_shape))+geom_point()+geom_violin(),
ggplot(df_beta_recovery_accuracy %>% group_by(beta_gamma_shape) %>% summarise(sensitivity=((acc == 'TP')/((acc == 'TP')+(acc == 'FN'))), specificity=((acc == 'TN')/((acc == 'TN')+(acc == 'FP')))),
       aes(x=beta_gamma_shape, y=specificity, group=beta_gamma_shape))+geom_point()+geom_violin(),
ncol=2)
dev.off()

ggplot(df_beta_recovery_accuracy_grouped,
       aes(x=beta_gamma_shape, y=sensitivity, group=beta_gamma_shape))+geom_point()+geom_violin()

df_beta_recovery_accuracy2 = df_beta_recovery_accuracy
# df_beta_recovery_accuracy2 = df_beta_recovery_accuracy2[!(is.na(df_beta_recovery_accuracy2$acc)),]

cat('Plotting accuracy 2\n')

try({
ggplot(df_beta_recovery_accuracy2, aes(x=beta_gamma_shape, fill=acc))+geom_bar(position='stack')+
  facet_wrap(.~interaction(beta_gamma_shape,df_beta_recovery_accuracy2$n),ncol=7 )

ggplot(df_beta_recovery_accuracy2, aes(x=beta_gamma_shape, group=interaction(n,beta_gamma_shape), fill=acc))+
  geom_bar(position='stack')+
  facet_wrap(.~(n),ncol=1 )


ggplot(df_beta_recovery_accuracy2, aes(x=factor(beta_gamma_shape), fill=acc))+
  # geom_bar(position='stack')+facet_wrap(.~n)
  geom_bar(position='fill')+facet_wrap(.~n)
ggsave(paste0(opt$output_folder_name, "/", opt$output_string, "_specificity_sensitivity_2.pdf"), height = 3, width = 8)
})

outfile_text <- paste0(opt$output_folder_name, "/", opt$output_string, "_results_info.txt")
cat('Creating output file for snakemake: ', outfile_text, '\n')
write(Sys.time(), file = outfile_text)

