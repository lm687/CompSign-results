#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd("~/Documents/PhD/GlobalDA/code/2_inference_TMB/mm_multinomial/")
library(TMB)
library(scales)
library(uuid)
library(ROCR)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(xtable)
library(mvtnorm)
source("../1_create_ROO/roo_functions.R")
source("../2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("mm_multinomial/helper_functions.R")
source("helper_TMB.R")
uuid = uuid::UUIDgenerate()
re_run_test = FALSE ## use cache
re_make_plots = TRUE
# set.seed(1245)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
# TMB::compile("mm_multinomial/ME_LNM.cpp", "-std=gnu++17")
# dyn.load(dynlib("mm_multinomial/ME_LNM"))
# TMB::compile("mm_multinomial/ME_multinomial.cpp", "-std=gnu++17")
# dyn.load(dynlib("mm_multinomial/ME_multinomial"))
# TMB::compile("mm_multinomial/ME_dirichletmultinomial.cpp", "-std=gnu++17")
# dyn.load(dynlib("mm_multinomial/ME_dirichletmultinomial"))

TMB::compile("mm_multinomial/fullRE_ME_multinomial.cpp")
dyn.load(dynlib("mm_multinomial/fullRE_ME_multinomial"))
TMB::compile("mm_multinomial/fullRE_ME_dirichletmultinomial.cpp")
dyn.load(dynlib("mm_multinomial/fullRE_ME_dirichletmultinomial"))

#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
##' Which dataset generation to analyse?
##' 20200625: first generation
##' (B) (not run) same as 20200625, but with non-zero intercept for beta
##' GenerationC: trying to create data that resembles better the PCAWG cohort data
dataset_generation='GenerationE'
dataset_generation=20200625
dataset_generation='GenerationD'
dataset_generation='GenerationC'
dataset_generation='GenerationC_norm'
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
datasets_files = list.files("../../data/assessing_models_simulation/datasets/", full.names = TRUE)
datasets_files = datasets_files[grep(pattern = dataset_generation, datasets_files)]
datasets = lapply(datasets_files, readRDS)
names(datasets) = gsub(".RDS", "", basename(datasets_files))
DA_bool = ( sapply(datasets, function(i) i$beta_gamma_shape) > 0 )
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
if(re_run_test){
  runs_M = mclapply(datasets_files, function(i)  try(wrapper_run_TMB(i, model = "fullRE_M", typedata = "simulation", simulation = TRUE)))
  names(runs_M) = gsub(".RDS", "", basename(datasets_files))
  saveRDS(object = runs_M, file = paste0("../../data/robjects_cache/tmb_results_simulations/", dataset_generation, "_simulation_runs_M_", uuid, ".RDS"))
}else{
  runs_M_names = list("../../data/robjects_cache/tmb_results_simulations/singleRE_20200625_simulation_runs_M_970f6533-5025-4eef-a598-956422846ddf.RDS",
                      "../../data/robjects_cache/tmb_results_simulations/GenerationC_simulation_runs_M_9c1d179c-8003-4398-923a-e08df2f4a439.RDS",
                      "../../data/robjects_cache/tmb_results_simulations/GenerationD_simulation_runs_M_96f86b9d-c9e5-4ce0-bc46-18503a851944.RDS",
                      "../../data/robjects_cache/tmb_results_simulations/GenerationE_simulation_runs_M_4c8d8557-8369-4c67-871d-ecea7ae2aa30.RDS",
                      "../../data/robjects_cache/tmb_results_simulations/GenerationC_norm_simulation_runs_M_ba7b1ea0-3ef3-4775-9a58-246e85627503.RDS")
  names(runs_M_names) = c(20200625, 'GenerationC', 'GenerationD', 'GenerationE', 'GenerationC_norm')
  runs_M = readRDS(runs_M_names[[as.character(dataset_generation)]])
}
if(re_run_test){
  runs_DM = mclapply(datasets_files, function(i)  try(wrapper_run_TMB(i, model = "fullRE_DM", typedata = "simulation", simulation = TRUE)))
  names(runs_DM) = gsub(".RDS", "", basename(datasets_files))
  saveRDS(object = runs_DM, file = paste0("../../data/robjects_cache/tmb_results_simulations/", dataset_generation, "_simulation_runs_DM_", uuid, ".RDS"))
}else{
  runs_DM_names = list("../../data/robjects_cache/tmb_results_simulations/singleRE_20200625_simulation_runs_DM_d8ea36c9-d0ed-4673-b86b-071de763dc65.RDS", #"../../data/robjects_cache/tmb_results_simulations/20200625_simulation_runs_DM_5c268769-daa7-4418-a41b-08c6ed9cf5d3.RDS",
                      "../../data/robjects_cache/tmb_results_simulations/GenerationC_simulation_runs_DM_9c1d179c-8003-4398-923a-e08df2f4a439.RDS",
                      "../../data/robjects_cache/tmb_results_simulations/GenerationD_simulation_runs_DM_96f86b9d-c9e5-4ce0-bc46-18503a851944.RDS",
                      "../../data/robjects_cache/tmb_results_simulations/GenerationE_simulation_runs_DM_4c8d8557-8369-4c67-871d-ecea7ae2aa30.RDS",
                      "../../data/robjects_cache/tmb_results_simulations/GenerationC_norm_simulation_runs_DM_1e2eaa10-e45a-4803-b827-16c95594dcf3.RDS")
  names(runs_DM_names) = c(20200625, 'GenerationC', 'GenerationD', 'GenerationE', 'GenerationC_norm')
  runs_DM = readRDS(runs_DM_names[[as.character(dataset_generation)]])
}
# if(re_run_test){
#   runs_LNM = lapply(datasets_files, function(i)  try( wrapper_run_TMB(i, model = "LNM", typedata = "simulation", simulation = TRUE)))
#   names(runs_LNM) = gsub(".RDS", "", basename(datasets_files))
#   saveRDS(object = runs_LNM, file = paste0("../../data/robjects_cache/tmb_results_simulations/", dataset_generation, "_simulation_runs_LNM_", uuid, ".RDS"))
# }else{
#   runs_LNM_names = list("../../data/robjects_cache/tmb_results_simulations/20200625_simulation_runs_LNM_5c268769-daa7-4418-a41b-08c6ed9cf5d3.RDS",
#                        "../../data/robjects_cache/tmb_results_simulations/GenerationC_simulation_runs_LNM_ae8b28fd-de16-4e2e-906e-9c9d49ab896c.RDS",
#                        "../../data/robjects_cache/tmb_results_simulations/GenerationD_simulation_runs_LNM_96f86b9d-c9e5-4ce0-bc46-18503a851944.RDS")
#   names(runs_LNM_names) = c(20200625, 'GenerationC', 'GenerationD')
#   runs_LNM = readRDS(runs_LNM_names[[as.character(dataset_generation)]])
# }
#-------------------------------------------------------------------------------------------------#

## How many failures have there been?
table(sapply(runs_DM, function(i) all(is.na(i$par.fixed))))
idx_to_repeat = 31
idx_to_repeat = 15
(sapply(runs_DM, function(i) all(is.na(i$par.fixed))))[idx_to_repeat]
# try(wrapper_run_TMB(datasets_files[[idx_to_repeat]], model = "fullRE_M", typedata = "simulation", simulation = TRUE))
grep('GenerationE_80_200_NA_4_0_dataset', datasets_files)
get_count_object_file(datasets_files[[idx_to_repeat]])$objects_counts@count_matrices_all
if(re_make_plots) createbarplot_object(datasets_files[[idx_to_repeat]])
#-------------------------------------------------------------------------------------------------#
# Keep only good runs
runs_M[sapply(runs_M, function(i) try(give_summary_per_sample(i))) != "Good"] = NA
runs_DM[sapply(runs_DM, function(i) try(give_summary_per_sample(i))) != "Good"] = NA
# runs_LNM[sapply(runs_LNM, function(i) try(give_summary_per_sample(i))) != "Good"] = NA
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
wrapper_run_ttest_ilr = function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  props = sapply(x, normalise_rw, simplify = FALSE)
  return(Compositional::hotel2T2(x1 = compositions::ilr(props[[1]]), x2 = compositions::ilr(props[[2]]))$pvalue)
}
wrapper_run_ttest_props = function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  props = sapply(x, normalise_rw, simplify = FALSE)
  return(Compositional::hotel2T2(x1 =props[[1]][,-1], x2 = props[[2]][,-1])$pvalue)
}
runs_ttest_irl = lapply(datasets_files, function(i)  try(wrapper_run_ttest_ilr(i)))
runs_ttest_props = lapply(datasets_files, function(i)  try(wrapper_run_ttest_props(i)))
pvals_ttest_ilr = as.numeric(unlist(runs_ttest_irl))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## assess if there was good convergence
sapply(runs_M, typeof)

## check if the names are the same
all(names(runs_M) == names(datasets))
all(names(runs_DM) == names(datasets))
# all(names(runs_LNM) == names(datasets))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
pvals_M = as.numeric(sapply(runs_M, function(i) try(wald_TMB_wrapper(i, verbatim=FALSE))))
pvals_DM = as.numeric(sapply(runs_DM, function(i) try(wald_TMB_wrapper(i, verbatim=FALSE))))
# pvals_LNM = as.numeric(sapply(runs_LNM, function(i) try(wald_TMB_wrapper(i, verbatim=FALSE))))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## adjust p-values
## TO DO
pvals_M_adj = pvals_M#*length(pvals_M)
pvals_DM_adj = pvals_DM#*length(pvals_DM)
# pvals_LNM_adj = pvals_LNM*length(pvals_LNM)
pvals_ttest_ilr_adj = pvals_ttest_ilr#*length(pvals_ttest_ilr)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## Basic analysis
## How many are said to be differentially abundant, how many actually are?
sapply(list(M=pvals_M_adj <= 0.05, DM=pvals_DM_adj <= 0.05,
            # LNM=pvals_LNM_adj <= 0.05,
            ttest=pvals_ttest_ilr_adj <= 0.05, true=DA_bool), table)

## Contingency table
table(DA_bool, M=pvals_M_adj <= 0.05)
table(DA_bool, DM=pvals_DM_adj <= 0.05)
table(DA_bool, LNM=pvals_LNM_adj <= 0.05)

## which datasets are TP, TN, FP and FN?
which_contingency = lapply(list(pvals_M_adj, pvals_DM_adj#, pvals_LNM_adj
                                ), function(i) list(TP = as.vector(which(DA_bool & (i <= 0.05))),
                                                                        TN = as.vector(which(!DA_bool & (i > 0.05))),
                                                                        FP = as.vector(which(!DA_bool & (i <= 0.05))),
                                                                        FN = as.vector(which(DA_bool & (i > 0.05)))))
names(which_contingency) = c('M', 'DM')#, 'LNM')

## It doesn't seem to have to do with the number of patients in the datasets
lapply(which_contingency, function(i) sapply(i, function(j) sapply(j, function(k) datasets[[k]]$n)))

#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## compare estimated betas to actual betas
for(j in which(sapply(runs_M, typeof) %in% c("logical", "character"))){
  runs_M[[j]] = list(par.fixed=c(beta=rep(NA, 2*(datasets[[j]]$d-1)))) ## *2 for slope and intercept
}

for(j in which(sapply(runs_DM, typeof) %in% c("logical", "character"))){
  runs_DM[[j]] = list(par.fixed=c(beta=rep(NA, 2*(datasets[[j]]$d-1)))) ## *2 for slope and intercept
}

# for(j in which(sapply(runs_LNM, typeof) %in% c("logical", "character"))){
#   runs_LNM[[j]] = list(par.fixed=c(beta=rep(NA, 2*(datasets[[j]]$d-1)))) ## *2 for slope and intercept
# }

## beta slopes are not well recovered, even though in fullinterceptRE_ME_multinomial.R it's done well!

true_slope_dataset = (sapply(datasets, function(i) i$beta[2,]))
inferred_slope_M = sapply(runs_M, function(i) try(select_slope_2(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE)))
par(mfrow=c(1,1))
plot(unlist(true_slope_dataset),
     unlist(inferred_slope_M))
abline(coef = c(0, 1))

#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## Analysing the characteristics of the DA datasets according to models
#' (to see if there's anything obsviously wrong about the way that data are generated and that could
#' explain poor accuracy of the tests)
par(mfrow=c(1,1))
boxplot(split(log(pvals_M), DA_bool)) ## p-values are generally lower in the truly-DA class.

## how to the betas look in each group? if they are not too different then simply THE SIMULATION IS DONE WRONG
boxplot(sapply(split(sapply(datasets, function(i) i$beta[2,]), DA_bool), function(j) as.vector(unlist(j))))
plot(sapply(split(sapply(datasets, function(i) i$beta[2,]), DA_bool), function(j) as.vector(unlist(j))))
boxplot(log(.2+sapply(split(sapply(datasets, function(i) i$beta[2,]), DA_bool), function(j) as.vector(unlist(j)))))


# some metric indicating change
metric_change = sapply(datasets, function(i){
  sum(abs(normalise_rw(i$W[1:i$n,]) - normalise_rw(i$W[(i$n+1):(2*(i$n)),])))/(i$n*i$d)
})
## here it shows that the "DA_bool" doesn't really correlate with how much of an effect there is
par(mfrow=c(1,2))
plot(metric_change, pvals_M_adj <= 0.05, col=(1+as.numeric(DA_bool)), pch=9)
plot(metric_change, pvals_DM_adj <= 0.05, col=(1+as.numeric(DA_bool)), pch=9)

## there are some betas which are VERY low in the differential abundant case
## that doesn't explain the false positives, though
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
if(re_make_plots){
  pdf(paste0("../../results/TMB/", dataset_generation, "_pvals_simulation.pdf"), height = 3, width = 10)
  par(mfrow=c(1,3))
  plot(sort(na.omit(pvals_M)), type = "h", main="M")
  plot(sort(na.omit(pvals_DM)), type = "h", main='DM')
  # plot(sort(na.omit(pvals_LNM)), type = "h", main="LNM")
  dev.off()
}
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
# How do the FP datasets look?
lapply(which_contingency, function(i) i$FP)
lapply(which_contingency, function(i) i$FN)

par(mfrow=c(1,2))
## example of a FN in all datasets
sapply(datasets[[6]]$objects_counts@count_matrices_all, function(i) image(t(i)))
createbarplot_ROOSigs(datasets[[6]]$objects_counts, slot = 'count_matrices_all', pre_path = "../../../CDA_in_Cancer/code/")
datasets[[6]]$beta_gamma_shape
datasets[[6]]$beta
matrix(python_like_select_name(runs_M[[6]]$par.fixed, 'beta'), nrow=2) 
pvals_DM[[6]]
which_beta_slope = select_slope_2(grep('beta', names(runs_M[[6]]$par.fixed)), verbatim = FALSE)

stat = t(matrix(select_slope_2(python_like_select_name(runs_M[[6]]$par.fixed, 'beta'), verbatim = FALSE))) %*% solve(runs_M[[6]]$cov.fixed[which_beta_slope,which_beta_slope]) %*% matrix(select_slope_2(python_like_select_name(runs_M[[6]]$par.fixed, 'beta'), verbatim = FALSE))
pchisq(stat, 2, lower.tail = FALSE)

## example of a FP in all
sapply(datasets[[3]]$objects_counts@count_matrices_all, function(i) image(t(i)))
createbarplot_ROOSigs(datasets[[3]]$objects_counts, slot = 'count_matrices_all', pre_path = "../../../CDA_in_Cancer/code/")
datasets[[3]]$beta
matrix(python_like_select_name(runs_M[[3]]$par.fixed, 'beta'), nrow=2) 

## plotting all datasets truly DA or not DA, to see if the distinction is obvious
plts = sapply(datasets, function(i) createbarplot_ROOSigs(i$objects_counts, slot = 'count_matrices_all', pre_path = "../../../CDA_in_Cancer/code/"))

## Possibilities of why the accuracy is so low
## a) I am not computing the p-value properly
## b) DA_bool is wrong, perhaps because beta_shape... is not the best thing to look at

if(re_make_plots){
  pdf(paste0("../../results/assessing_models/barplots_datasets", dataset_generation, "_DA.pdf"), width = 20, height = 15)
  do.call('grid.arrange', plts[DA_bool]) ## differentially abundant
  dev.off()
  pdf(paste0("../../results/assessing_models/barplots_datasets", dataset_generation, "_nonDA.pdf"), width = 20, height = 15)
  do.call('grid.arrange', plts[!DA_bool]) ## not differentially abundant
  dev.off()
}

## the ones that DM think are DA or not
# do.call('grid.arrange', plts[which(pvals_DM_adj <= 0.05)]) ## differentially abundant
# do.call('grid.arrange', plts[-which(pvals_DM_adj <= 0.05)]) ## not differentially abundant
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## Remove outlier datasets
remove_outliers=F
if(remove_outliers){
  plot(density(sapply(datasets, function(i) max(i$beta[2,]))))
  remove_idx_datasets = which(sapply(datasets, function(i) max(i$beta[2,])) > 5)
  if(length(remove_idx_datasets) > 0){
    ## there is some problematic dataset
    datasets = datasets[-remove_idx_datasets]
    pvals_M_adj = pvals_M_adj[-remove_idx_datasets]
    pvals_DM_adj = pvals_DM_adj[-remove_idx_datasets]
    pvals_LNM_adj = pvals_LNM_adj[-remove_idx_datasets]
    pvals_ttest_ilr_adj = pvals_ttest_ilr_adj[-remove_idx_datasets]
    runs_M = runs_M[-remove_idx_datasets]
    runs_DM = runs_DM[-remove_idx_datasets]
    runs_LNM = runs_LNM[-remove_idx_datasets]
    DA_bool = DA_bool[-remove_idx_datasets]
  }
}

summarise_DA_detection = function(true, predicted){
  ## remove NAs
  which_na = which(is.na(predicted))
  if(length(which_na) > 0){ ## some NA
    true = true[-which_na]
    predicted = predicted[-which_na]
  }
  
  FPs = sum(!true & predicted)/sum(predicted)
  TPs = sum(true & predicted)/sum(predicted)
  TNs = sum(!true & !predicted)/sum(!predicted)
  FNs = sum(true & !predicted)/sum(!predicted)
  total_pos = sum(true | predicted)
  Power = TPs/total_pos
  Sensitivity = TPs / (TPs + FNs)
  Specificity = TNs / (TNs + FPs)
  pred <- prediction(as.numeric(true), as.numeric(predicted))
  AUC = try(performance(pred, "auc")@y.values[[1]])
  return(c(FP=FPs, TP=TPs, Power=Power, AUC=AUC, Specificity=Specificity, Sensitivity=Sensitivity))
}

## they are all pretty horrendous in terms of FP
summarise_DA_detection(true = DA_bool, predicted = pvals_M_adj <= 0.05)
summarise_DA_detection(true = DA_bool, predicted = pvals_DM_adj <= 0.05)
summarise_DA_detection(true = DA_bool, predicted = pvals_LNM_adj <= 0.05)
summarise_DA_detection(true = DA_bool, predicted = pvals_ttest_ilr_adj <= 0.05)

all_summary = do.call('rbind', list(M=summarise_DA_detection(true = DA_bool, predicted = pvals_M_adj <= 0.05),
                                    DM=summarise_DA_detection(true = DA_bool, predicted = pvals_DM_adj <= 0.05),
                                    #LNM=summarise_DA_detection(true = DA_bool, predicted = pvals_LNM_adj <= 0.05),
                                    ttest=summarise_DA_detection(true = DA_bool, predicted = pvals_ttest_ilr_adj <= 0.05)))
xtable::xtable(all_summary, digits=c(0,4,4,4,4,4,4))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## Creating tables of prediction accuracy with balanced datasets
## I need to split by beta_gamma_shape

table(sapply(datasets, function(i) i$beta_gamma_shape))
splits = lapply(list(DA_bool, pvals_M_adj, pvals_DM_adj, pvals_ttest_ilr_adj), function(vec)
  split(vec,f = sapply(datasets, function(i) i$beta_gamma_shape))
  )
names(splits) = c('DA_bool', 'pvals_M_adj', 'pvals_DM_adj', 'pvals_ttest_ilr_adj')

## get groups of beta_gamma_shape=0 with any of the others

all_summary_balanced = lapply(2:length(splits[[1]]), function(comparison_idx){
  DA_bool_comparison = unlist(splits$DA_bool[c(1, comparison_idx)])
  do.call('rbind', list(M=summarise_DA_detection(true = DA_bool_comparison,
                                                 predicted = unlist(splits$pvals_M_adj[c(1, comparison_idx)]) <= 0.05),
                        # DM=summarise_DA_detection(true = DA_bool_comparison,
                        #                           predicted =  unlist(splits$pvals_DM_adj[c(1, comparison_idx)]) <= 0.05),
                        # LNM=summarise_DA_detection(true = DA_bool_comparison,
                        #                            predicted =  unlist(splits$pvals_LNM_adj[c(1, comparison_idx)]) <= 0.05),
                        ttest=summarise_DA_detection(true = DA_bool_comparison,
                                                     predicted =  unlist(splits$pvals_ttest_ilr_adj[c(1, comparison_idx)]) <= 0.05)))
})
all_summary_balanced
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## it should be the case that the number of TRUE increases
sapply(splits$pvals_DM_adj, function(i) table(i <= 0.05))
## problem: still too many false positives! there are 

## fraction of false positives in the first group
sapply(splits, function(i) sum(na.omit(i[[1]] <= 0.05)) /(sum(!is.na(i[[1]]))))

## here there should be a negative correlation (higher effect size <=> lower p-value) (which there is to some extend)
par(mfrow=c(1,3))
plot( sapply(datasets, function(i) i$beta_gamma_shape), pvals_M_adj)
plot( sapply(datasets, function(i) i$beta_gamma_shape), pvals_DM_adj)
plot( sapply(datasets, function(i) i$beta_gamma_shape), pvals_LNM_adj)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## here there is only the intercept
df_beta_recovery = cbind.data.frame(beta_true = unlist(sapply(datasets, function(i) i$beta[2,])),
                                    idx = rep(1:length(datasets) , unlist(sapply(datasets, function(i) i$d))-1),
                                    d =  rep(unlist(sapply(datasets, function(i) i$d)), unlist(sapply(datasets, function(i) i$d))-1),
                                    n =  rep(unlist(sapply(datasets, function(i) i$n)), unlist(sapply(datasets, function(i) i$d))-1),
                                    beta_gamma_shape =  rep(unlist(sapply(datasets, function(i) i$beta_gamma_shape)), unlist(sapply(datasets, function(i) i$d))-1),
                                    beta_est_M = unlist(sapply(runs_M, function(i) select_slope_2(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE))),
                                    beta_stderr_M = unlist(sapply(runs_M, give_stderr)),
                                    beta_est_DM = unlist(sapply(runs_DM, function(i) select_slope_2(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE))),
                                    beta_stderr_DM = unlist(sapply(runs_DM, give_stderr)),
                                    # beta_est_LNM = unlist(sapply(runs_LNM, function(i) select_slope_2(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE))),
                                    pvals_M_adj=rep(pvals_M_adj, unlist(sapply(datasets, function(i) i$d))-1),
                                    pvals_DM_adj=rep(pvals_DM_adj, unlist(sapply(datasets, function(i) i$d))-1),
                                    # pvals_LNM_adj=rep(pvals_LNM_adj, unlist(sapply(datasets, function(i) i$d))-1)
                                    DA_bool=rep(DA_bool, unlist(sapply(datasets, function(i) i$d))-1),
                                    idx_within_dataset=unlist(sapply(datasets, function(i) 1:(i$d-1))))
df_beta_recovery$bool_zero_true_beta = factor(df_beta_recovery$beta_true == 0, levels=c(TRUE, FALSE))
pairs(df_beta_recovery[,c('pvals_M_adj', 'pvals_DM_adj'#, 'pvals_LNM_adj'
                          )], main='Pairs plot of p-values')
pairs(df_beta_recovery[,c('beta_est_M', 'beta_est_DM'#, 'beta_est_LNM'
                          )], main='Pairs plot of betas')

table(na.omit(df_beta_recovery$pvals_M_adj <= 0.05))
table(na.omit(df_beta_recovery$pvals_DM_adj <= 0.05))
table(na.omit(df_beta_recovery$pvals_LNM_adj <= 0.05))

#-------------------------------------------------------------------------------------------------#
## zero slopes are not well detected
ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_M), col=beta_gamma_shape))+geom_point()+
  geom_abline(intercept = 0, slope = 1)

# ggplot(df_beta_recovery,
#        aes(x=(beta_true), y=(beta_est_M), col=beta_gamma_shape))+geom_point()+
#   geom_abline(intercept = 0, slope = 1)+facet_wrap(.~idx,scales = "free")
# 
ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_DM), col=beta_gamma_shape))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~idx,scales = "free")

# ggplot(df_beta_recovery,
#        aes(x=(beta_true), y=(beta_est_LNM), col=beta_gamma_shape))+geom_point()+
#   geom_abline(intercept = 0, slope = 1)+facet_wrap(.~idx,scales = "free")
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## plotting separately the differentially abundant (right) and non-differentially abundant (left) datasets
ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_M), col=n))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")

## Multinomial
p <- ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_M), col=(pvals_M_adj<0.05)))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")
gp <- ggplotGrob(p); facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]; gp$widths[facet.columns] <- gp$widths[facet.columns] * c(1,4); grid::grid.draw(gp)
ggsave(paste0("../../results/assessing_models/recovery_betaslope_M_", dataset_generation, ".pdf"))

ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_DM), col=(pvals_DM_adj<0.05)))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")
ggsave(paste0("../../results/assessing_models/recovery_betaslope_DM_", dataset_generation, ".pdf"))

## Looking at true zeros
## i.e. runs where all beta slopes are zero
if(re_make_plots){
  ggplot(df_beta_recovery[!(df_beta_recovery$DA_bool),],
         aes(x=(beta_true), y=(beta_est_M), col=(pvals_M_adj<0.05)))+geom_point()+
    geom_abline(intercept = 0, slope = 1)+facet_wrap(.~(pvals_M_adj<0.05), scales = "free_x")
  ggsave(paste0("../../results/TMB/", dataset_generation, "_M_betaslopes_nonDA.pdf"), height = 3, width = 8)
  
  ggplot(df_beta_recovery[!(df_beta_recovery$DA_bool),],
         aes(x=idx_within_dataset, y=(beta_est_M), col=(pvals_M_adj<0.05)))+
    geom_point(data = df_beta_recovery[!(df_beta_recovery$DA_bool),], aes(x=idx_within_dataset, y=beta_true), shape=4)+
    geom_abline(slope = 0, intercept = 0, alpha=0.2)+
    geom_point()+
    geom_errorbar(aes(ymin=beta_est_M-beta_stderr_M, ymax=beta_est_M+beta_stderr_M), width=.2,
                  position=position_dodge(.9))+theme_bw()+
    facet_wrap(.~idx, scales='free_x')+theme(legend.position = "bottom")
  ggsave(paste0("../../results/TMB/", dataset_generation, "_M_betaslopes_stderr_nonDA.pdf"), height = 8, width = 8)
  ggplot(df_beta_recovery[(df_beta_recovery$DA_bool),],
         aes(x=idx_within_dataset, y=(beta_est_M), col=(pvals_M_adj<0.05)))+
    geom_point(data = df_beta_recovery[df_beta_recovery$DA_bool,], aes(x=idx_within_dataset, y=beta_true), shape=4)+
    geom_abline(slope = 0, intercept = 0, alpha=0.2)+geom_point()+
    geom_errorbar(aes(ymin=beta_est_M-beta_stderr_M, ymax=beta_est_M+beta_stderr_M), width=.2,
                  position=position_dodge(.9))+
    facet_wrap(.~interaction(idx, beta_gamma_shape), scales='free_x', nrow=length(unique(df_beta_recovery$beta_gamma_shape)))+
    theme_bw()+theme(legend.position = "bottom")
  ggsave(paste0("../../results/TMB/", dataset_generation, "_M_betaslopes_stderr_DA.pdf"), height = 8, width = 8)
}

## Dirichlet-Multinomial
## there are some very extreme values of beta
plot(density(log(abs(na.omit(df_beta_recovery$beta_est_DM))))) ## all betas
plot(density(na.omit(sapply(runs_DM, function(i) max(python_like_select_name(i$par.fixed, "beta"))))))
df_beta_recovery[df_beta_recovery$idx %in% which(sapply(runs_DM, function(i) max(abs(python_like_select_name(i$par.fixed, "beta")))) > 10),'beta_est_DM'] = NA
p <- ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_DM), col=(pvals_DM_adj)<0.05))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")
gp <- ggplotGrob(p); facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]; gp$widths[facet.columns] <- gp$widths[facet.columns] * c(1,4); grid::grid.draw(gp)

if(re_make_plots){
  ggplot(df_beta_recovery[!(df_beta_recovery$DA_bool),],
         aes(x=(beta_true), y=(beta_est_DM), col=(pvals_DM_adj<0.05)))+geom_point()+
    geom_abline(intercept = 0, slope = 1)+facet_wrap(.~(pvals_DM_adj<0.05), scales = "free_x")
  ggsave(paste0("../../results/TMB/", dataset_generation, "_DM_betaslopes_nonDA.pdf"), height = 3, width = 8)
  
  ggplot(df_beta_recovery[!(df_beta_recovery$DA_bool),],
         aes(x=idx_within_dataset, y=(beta_est_DM), col=(pvals_DM_adj<0.05)))+
    geom_point(data = df_beta_recovery[!(df_beta_recovery$DA_bool),], aes(x=idx_within_dataset, y=beta_true), shape=4)+
    geom_abline(slope = 0, intercept = 0, alpha=0.2)+
    geom_point()+
    geom_errorbar(aes(ymin=beta_est_DM-beta_stderr_DM, ymax=beta_est_DM+beta_stderr_DM), width=.2,
                  position=position_dodge(.9))+theme_bw()+
    facet_wrap(.~idx, scales='free_x')+theme(legend.position = "bottom")
  ggsave(paste0("../../results/TMB/", dataset_generation, "_DM_betaslopes_stderr_nonDA.pdf"), height = 8, width = 8)
  ggplot(df_beta_recovery[(df_beta_recovery$DA_bool),],
         aes(x=idx_within_dataset, y=(beta_est_DM), col=(pvals_DM_adj<0.05)))+
    geom_point(data = df_beta_recovery[df_beta_recovery$DA_bool,], aes(x=idx_within_dataset, y=beta_true), shape=4)+
    geom_abline(slope = 0, intercept = 0, alpha=0.2)+geom_point()+
    geom_errorbar(aes(ymin=beta_est_DM-beta_stderr_DM, ymax=beta_est_DM+beta_stderr_DM), width=.2,
                  position=position_dodge(.9))+
    facet_wrap(.~interaction(idx, beta_gamma_shape), scales='free_x', nrow=length(unique(df_beta_recovery$beta_gamma_shape)))+
    theme_bw()+theme(legend.position = "bottom")
  ggsave(paste0("../../results/TMB/", dataset_generation, "_DM_betaslopes_stderr_DA.pdf"), height = 8, width = 8)
}

## LNM
# df_beta_recovery[df_beta_recovery$idx %in% which(sapply(runs_LNM, function(i) max(abs(python_like_select_name(i$par.fixed, "beta")))) > 20),'beta_est_LNM'] = NA
# p <- ggplot(df_beta_recovery,
#        aes(x=(beta_true), y=(beta_est_LNM), col=(pvals_LNM_adj)<0.05))+geom_point()+
#   geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")
# gp <- ggplotGrob(p); facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]; gp$widths[facet.columns] <- gp$widths[facet.columns] * c(1,4); grid::grid.draw(gp)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## Do betas compensate for each other? i.e. if there is an overestimate in slope, is there an overestimate
## in the intersect corresponding to the same log-ratio?

beta_gamma_shapes_df = cbind.data.frame(idx=1:length(datasets), beta_gamma_shape= sapply(datasets, function(i) i$beta_gamma_shape))
cut_beta_gamma_shapes = split(1:length(datasets), f = sapply(datasets, function(i) i$beta_gamma_shape))

df_beta_recovery$beta_gamma_shape = beta_gamma_shapes_df[match(df_beta_recovery$idx, beta_gamma_shapes_df$idx),'beta_gamma_shape']

lapply(cut_beta_gamma_shapes, function(i)  df_beta_recovery[match(i, df_beta_recovery$idx),c('beta_true', 'beta_est_M')])

if(re_make_plots){
  ggplot(df_beta_recovery,
         aes(x=(beta_true), y=(beta_est_M), col=n))+geom_point()+
    geom_abline(intercept = 0, slope = 1)+facet_wrap(.~beta_gamma_shape, scales = "free_x", ncol=5)
  ggsave(paste0("../../results/TMB/", dataset_generation, "_recoveries_betashape_facet.pdf"), height = 3, width = 11)
  
  ggplot(df_beta_recovery,
         aes(x=(beta_true), y=(beta_est_M), col=n))+geom_point()+
    geom_abline(intercept = 0, slope = 1)+facet_wrap(.~beta_gamma_shape, scales = "free", ncol=5)
  ggsave(paste0("../../results/TMB/", dataset_generation, "_recoveries_betashape_facet_freescales.pdf"), height = 3, width = 11)
}  
## Problem when beta shape is 0! (why?)

#-------------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------------#
## colour code per patient
ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_M), col=factor(idx)))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+guides(color=FALSE)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## there are some simulated betas which are extreme and I have removed those

## a clear problem is some datasets with no zero or near-zero differential abundance having a very small p-value
df_beta_recovery[df_beta_recovery$beta_true == 0,][1:10,]
## e.g. idx=6
df_beta_recovery[df_beta_recovery$idx == 6,]


simulate_theta_from_TMB_M = function(tmb_object, dataset_object){
  ## there is a column of zeros because the intercepts have been absorbed by the random effects intercept
  theta_sim = softmax(
    cbind(t(dataset_object$X_sim) %*% t(matrix(python_like_select_name(tmb_object$par.fixed, 'beta'), ncol = 2)) +
    replicate(n = dataset_object$d-1, tmb_object$par.random %*% dataset_object$Z_sim, simplify = TRUE), 0)
    )
  list(theta_sim[1:dataset_object$n,], theta_sim[(dataset_object$n+1):(2*dataset_object$n),])
}

simulate_theta_from_TMB_DM = function(tmb_object, dataset_object){
  ## there is a column of zeros because the intercepts have been absorbed by the random effects intercept
  alphabar_sim = 
    t(dataset_object$X_sim) %*% t(matrix(python_like_select_name(tmb_object$par.fixed, 'beta'), ncol = 2)) +
            replicate(n = dataset_object$d-1, tmb_object$par.random %*% dataset_object$Z_sim, simplify = TRUE)
  alpha = softmax(cbind(alphabar_sim * exp(rep(python_like_select_name(tmb_object$par.fixed, 'log_lambda'),
                                         each=nrow(dataset_object$Z_sim))), 0))
  theta_sim = t(apply(alpha, 1, function(i) MCMCpack::rdirichlet(1, i)))
  list(theta_sim[1:dataset_object$n,], theta_sim[(dataset_object$n+1):(2*dataset_object$n),])
}

simulate_theta_from_TMB_DM_fullRE = function(tmb_object, dataset_object){
  ## there is a column of zeros because the intercepts have been absorbed by the random effects intercept
  # alphabar_sim = 
  #   t(dataset_object$X_sim) %*% t(matrix(python_like_select_name(tmb_object$par.fixed, 'beta'), ncol = 2)) +
  #   replicate(n = dataset_object$d-1, tmb_object$par.random %*% dataset_object$Z_sim, simplify = TRUE)
  # alpha = softmax(cbind(alphabar_sim * exp(rep(python_like_select_name(tmb_object$par.fixed, 'log_lambda'),
  #                                              each=nrow(dataset_object$Z_sim))), 0))
  # theta_sim = t(apply(alpha, 1, function(i) MCMCpack::rdirichlet(1, i)))
  # list(theta_sim[1:dataset_object$n,], theta_sim[(dataset_object$n+1):(2*dataset_object$n),])
}

simulate_theta_from_TMB_LNM = function(tmb_object, dataset_object){
  ## there is a column of zeros because the intercepts have been absorbed by the random effects intercept
  mu_sim = t(dataset_object$X_sim) %*% t(matrix(python_like_select_name(tmb_object$par.fixed, 'beta'), ncol = 2)) +
            replicate(n = dataset_object$d-1, tmb_object$par.random %*% dataset_object$Z_sim, simplify = TRUE)
  
  cov_matrix = matrix(NA, nrow=(dataset_object$d-1), ncol=(dataset_object$d-1))
  fill_cov_matrix = matrix(python_like_select_name(tmb_object$par.fixed, 'cov_par'))
  cov_matrix = give_UNSTRUCTURED_CORR_t_matrix(vec = fill_cov_matrix, dim_mat = (dataset_object$d-1))
  ## fill in the order that TMB's UNSTRUCTURED_CORR_t saves the covariances
  # cov_matrix[unlist(sapply(2:nrow(cov_matrix), function(rw) seq(from = rw,length.out = (rw-1), by = nrow(cov_matrix) )))] = fill_cov_matrix
  # cov_matrix[unlist(sapply(2:nrow(cov_matrix), function(cl) seq(from = (((cl-1)*nrow(cov_matrix))+1),length.out = (cl-1), by = 1 )))] = fill_cov_matrix
  # diag(cov_matrix) = 1 ## because this is how it was simulated
  theta_sim2 = softmax(cbind(t(apply(mu_sim, 1, function(mu_it) mvtnorm::rmvnorm(n = 1, mean = mu_it, sigma = cov_matrix))),
                     0))
  list(theta_sim2[1:dataset_object$n,], theta_sim2[(dataset_object$n+1):(2*dataset_object$n),])
}

simulate_theta_from_TMB_LNM_fullRE = function(tmb_object, dataset_object){
  # ## there is a column of zeros because the intercepts have been absorbed by the random effects intercept
  # mu_sim = t(dataset_object$X_sim) %*% t(matrix(python_like_select_name(tmb_object$par.fixed, 'beta'), ncol = 2)) +
  #   replicate(n = dataset_object$d-1, tmb_object$par.random %*% dataset_object$Z_sim, simplify = TRUE)
  # 
  # cov_matrix = matrix(NA, nrow=(dataset_object$d-1), ncol=(dataset_object$d-1))
  # fill_cov_matrix = matrix(python_like_select_name(tmb_object$par.fixed, 'cov_par'))
  # cov_matrix = give_UNSTRUCTURED_CORR_t_matrix(vec = fill_cov_matrix, dim_mat = (dataset_object$d-1))
  # theta_sim2 = softmax(cbind(t(apply(mu_sim, 1, function(mu_it) mvtnorm::rmvnorm(n = 1, mean = mu_it, sigma = cov_matrix))),
  #                            0))
  # list(theta_sim2[1:dataset_object$n,], theta_sim2[(dataset_object$n+1):(2*dataset_object$n),])
}

par(mfrow=c(4,2))
sapply(datasets[[6]][[1]]@count_matrices_all, function(i) image(t(i), main='true'))
sapply(simulate_theta_from_TMB_M(runs_M[[6]], datasets[[6]]) , function(i) image(t(i), main='sim M est'))
sapply(simulate_theta_from_TMB_DM(runs_DM[[6]], datasets[[6]]) , function(i) image(t(i), main='sim DM est'))
sapply(simulate_theta_from_TMB_LNM(runs_LNM[[6]], datasets[[6]]) , function(i) image(t(i), main='sim LNM est'))

simulate_theta_from_TMB_DM(runs_M[[6]], datasets[[6]])

## compare the true and simulated parameters

## (option A) compare theta and normalised counts directly
plot(as.vector(normalise_rw(datasets[[6]]$W)),
     as.vector(do.call('rbind', simulate_theta_from_TMB_M(runs_M[[6]], datasets[[6]]))))



par(mfrow=c(1,3))
## no variability
plot(as.vector(replicate(100, as.vector(do.call('rbind', simulate_theta_from_TMB_M(runs_M[[6]], datasets[[6]]))))),
     as.vector(replicate(100, as.vector(normalise_rw(datasets[[6]]$W)))))
abline(coef = c(0,1), col='blue', lty='dashed')

## overdispersed
sim_multiple_DM = list((replicate(100, as.vector(do.call('rbind', simulate_theta_from_TMB_DM(runs_DM[[6]], datasets[[6]]))))),
                       (replicate(100, as.vector(normalise_rw(datasets[[6]]$W)))))
plot(as.vector(sim_multiple_DM[[1]]),
     as.vector(sim_multiple_DM[[2]]))
abline(coef = c(0,1), col='blue', lty='dashed')

## overdispersed
sim_multiple_LNM = list((replicate(100, as.vector(do.call('rbind', simulate_theta_from_TMB_LNM(runs_LNM[[6]], datasets[[6]]))))),
                        (replicate(100, as.vector(normalise_rw(datasets[[6]]$W)))))
plot(as.vector(sim_multiple_LNM[[1]]),
     as.vector(sim_multiple_LNM[[2]]))
abline(coef = c(0,1), col='blue', lty='dashed')

## with same data, create the confidence intervals
quantiles_DM = apply(sim_multiple_DM[[1]], 1, function(i) quantile(i, probs=c(0.025, 0.975), na.rm = TRUE) )
quantiles_LNM = apply(sim_multiple_LNM[[1]], 1, function(i) quantile(i, probs=c(0.025, 0.975), na.rm = TRUE) )

## many NAs!
sum(sapply(1:nrow(sim_multiple_DM[[1]]), function(j){
  (sim_multiple_DM[[2]][j,1] > quantiles_DM[1,j] ) & (sim_multiple_DM[[2]][j,1] < quantiles_DM[2,j])
}), na.rm = TRUE)/nrow(sim_multiple_DM[[1]])

cred_int_theta_sim = function(model, idx_dataset){
  if(typeof(runs_DM[[idx_dataset]]) == "character"){
    return(NA)
  }else{
    if(model == 'DM'){
      sim_multiple = list((replicate(100, as.vector(do.call('rbind', simulate_theta_from_TMB_DM(runs_DM[[idx_dataset]], datasets[[idx_dataset]]))))),
                             (replicate(100, as.vector(normalise_rw(datasets[[idx_dataset]]$W)))))
    }else if(model == 'LNM'){
      sim_multiple = list((replicate(100, as.vector(do.call('rbind', simulate_theta_from_TMB_LNM(runs_LNM[[idx_dataset]], datasets[[idx_dataset]]))))),
                              (replicate(100, as.vector(normalise_rw(datasets[[idx_dataset]]$W)))))
    }else{
      stop('Incorrect <model> specification')
    }
  
    quantiles_sim = apply(sim_multiple[[1]], 1, function(i) quantile(i, probs=c(0.025, 0.975), na.rm = TRUE) )
    
    return(sum(sapply(1:nrow(sim_multiple[[1]]), function(j){
      (sim_multiple[[2]][j,1] > quantiles_sim[1,j] ) & (sim_multiple[[2]][j,1] < quantiles_sim[2,j])
    }), na.rm = TRUE)/nrow(sim_multiple[[1]]))
  }
}

# cred_int_theta_sims_DM = sapply(1:length(datasets), cred_int_theta_sim, model = 'DM')
# saveRDS(cred_int_theta_sims_DM, paste0("../../data/robjects_cache/cred_int_theta_sims_DM_", dataset_generation, ".RDS"))
# cred_int_theta_sims_LNM = sapply(1:length(datasets), cred_int_theta_sim, model = 'LNM')
# saveRDS(cred_int_theta_sims_LNM, paste0("../../data/robjects_cache/cred_int_theta_sims_LNM_", dataset_generation, ".RDS"))
cred_int_theta_sims_DM = readRDS("../../data/robjects_cache/cred_int_theta_sims_DM.RDS")
cred_int_theta_sims_LNM = readRDS("../../data/robjects_cache/cred_int_theta_sims_LNM.RDS")
cred_int_theta_sims_DM = readRDS("../../data/robjects_cache/cred_int_theta_sims_DM_GenerationD.RDS")
cred_int_theta_sims_LNM = readRDS("../../data/robjects_cache/cred_int_theta_sims_LNM_GenerationD.RDS")

plot(cred_int_theta_sims_DM, cred_int_theta_sims_LNM)

df_cred_int_theta = data.frame(cred_int_theta_sims_DM=cred_int_theta_sims_DM,
           cred_int_theta_sims_LNM=cred_int_theta_sims_LNM,
           d=unlist(sapply(datasets, `[`, "d")),
           n=unlist(sapply(datasets, `[`, "n")),
           DM_convergence=as.character(sapply(runs_DM, `[`, "pdHess")),
           LNM_convergence=as.character(sapply(runs_LNM, `[`, "pdHess")),
           names=names(datasets)
)
which(df_cred_int_theta$LNM_convergence == "NULL")
which(df_cred_int_theta$DM_convergence == "NULL")

df_cred_int_theta = df_cred_int_theta[!(df_cred_int_theta$LNM_convergence == "NULL" | df_cred_int_theta$DM_convergence == "NULL"),]

if(re_make_plots){
  ggplot(df_cred_int_theta,
         aes(x=cred_int_theta_sims_DM, y=cred_int_theta_sims_LNM, col=d))+ geom_point()
  ggplot(df_cred_int_theta,
         aes(x=cred_int_theta_sims_DM, y=cred_int_theta_sims_LNM, col=n))+ geom_point()
  ggplot(df_cred_int_theta,
         aes(x=cred_int_theta_sims_DM, y=cred_int_theta_sims_LNM, col=DM_convergence, shape=LNM_convergence))+ geom_point()+
    geom_abline(slope = 1, intercept = 0)
  ggsave("../../results/assessing_models/theta_in_credint_DM_and_LNM.pdf", width = 7, height = 6)
  ## pairing the same cancer types
  ggplot(melt(df_cred_int_theta, id.vars = c('d', 'n', 'DM_convergence', 'LNM_convergence', 'names')),
         aes(x=variable, y=as.numeric(value), col=DM_convergence, group=names, col=DM_convergence))+ geom_line()+facet_wrap(.~factor(DM_convergence, levels=c('TRUE', 'FALSE')))
  ggsave("../../results/assessing_models/theta_in_credint_DM_and_LNM_2.pdf", width = 9, height = 5)
  ggplot(df_cred_int_theta,
         aes(x=cred_int_theta_sims_DM, y=cred_int_theta_sims_LNM, col=LNM_convergence, shape=DM_convergence))+ geom_point()
}

hist(cred_int_theta_sims_DM)

## (option B) we simulate theta, now we have to draw from multinomial
sim_M = t(sapply(1:nrow(datasets[[6]]$W), function(i){
  rmultinom(n = 1,
            prob = do.call('rbind', simulate_theta_from_TMB_M(runs_M[[6]], datasets[[6]]))[i,],
            size = rowSums(datasets[[6]]$W)[i])
}))
par(mfrow=c(1,1))
plot(as.vector(datasets[[6]]$W), as.vector(sim_M))

## another problem is clearly differential abundance runs having a very large p-value
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## Overdispersion coefficients for DM runs
plot(do.call('rbind', sapply(runs_DM, function(i) python_like_select_name(i$par.fixed, 'log_lambda'))))
if(re_make_plots){
  ggplot(data.frame(t(sapply(runs_DM, function(i){
    if(length(python_like_select_name(i$par.fixed, 'log_lambda')) == 0){
    c(NA, NA)
    }else{
      python_like_select_name(i$par.fixed, 'log_lambda')
    }
    })),
    converged=unlist(sapply(runs_DM, function(i)if(!is.null(i$pdHess)){i$pdHess}else{NA}))  ), aes(x=log_lambda, y=log_lambda.1, col=converged))+geom_point()
}
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
# Compare to stan runs
folder_stan = "../../data/assessing_models_simulation/inference_results/"
stanfiles = list.files(folder_stan, full.names = TRUE)
load_stan = function(stanfile){
  load(stanfile)
  return(list(W=W, X=X, Z=Z, d=d, fit_stan=fit_stan))
}
# names(stan_simA) = gsub(".RData", "", basename(stanfiles))
# rhats = sapply(stan_simA, function(i) rstan::Rhat(as.matrix(i$fit_stan)))

stan_simA_M = gsub(".RData", "", basename(stanfiles[grepl("ModelMROO", stanfiles)]))
stan_simA_DM = stanfiles[grepl("ModelDMROO", stanfiles)]
stan_simA_DM = gsub(".RData", "", basename(stan_simA_DM[!grepl("Nits1000_ModelDMROO", stan_simA_DM)]))
stan_simA_LNM = gsub(".RData", "", basename(stanfiles[grepl("ModelLNMROO", stanfiles)]))

# min(rhats)
## supossing this is well computed, none of these runs converged

modify_M_name = function(df){
  df %>% gsub("_nlambda", "_", .) %>% gsub("_d", "_", .) %>% gsub("_n", "_", .) %>% 
  gsub("beta_intensity", "", .) %>% 
  gsub("_Nits1000_ModelMROO", "", .) 
}
modify_DM_name = function(df){
  df %>% gsub("_nlambda", "_", .) %>% gsub("_d", "_", .) %>% gsub("_n", "_", .) %>% 
  gsub("beta_intensity", "", .) %>% gsub("lambda", "", .) %>% 
  gsub("_Nits2000_ModelDMROO", "", .) 
}
modify_LNM_name = function(df){
  df %>% gsub("_nlambda", "_", .) %>% gsub("_d", "_", .) %>% gsub("_n", "_", .) %>% 
  gsub("beta_intensity", "", .) %>% 
  gsub("_Nits2000_ModelLNMROO", "", .) 
}
mod_names_stan_M = modify_M_name(stan_simA_M)
mod_names_stan_DM = modify_DM_name(stan_simA_DM)
mod_names_stan_LNM = modify_DM_name(stan_simA_LNM)
mod_names_M = gsub("_dataset", "", names(runs_M))
mod_names_DM = gsub("_dataset", "", names(runs_DM))
mod_names_LNM = gsub("_dataset", "", names(runs_LNM))

modify_M_name_2 = function(i){
  paste0(paste0(strsplit(i, '_')[[1]][1:3], collapse='_'),
         '_10_',
         paste0(strsplit(i, '_')[[1]][4:5], collapse='_'))}
modify_DM_name_2 = function(i){
  paste0(paste0(strsplit(i, '_')[[1]][1:3], collapse='_'),
         '_10_',
         paste0(strsplit(i, '_')[[1]][4:5], collapse='_'))}
if(dataset_generation == 20200625){
  ## there is a bit missing in mod_names_stan_M compared to mod_names_M, which is always a 10. I am
  ## not sure if it's the coefficient for overdispersion but I add it to mod_names_stan_M anyway
  mod_names_stan_M = sapply(mod_names_stan_M, modify_M_name_2)
  mod_names_stan_DM = sapply(mod_names_stan_DM, modify_DM_name_2)
}

# dim(python_like_select_colnames(as.matrix(stan_simA[[1]]$fit_stan), "beta\\[2,"))

##' I am not sure why the TMB (mod_names_DM) contains the number of features but the stan 
##' runs do not. for the match I remove the features from the mod_names_DM names
stan_simA_M = stan_simA_M[match(mod_names_M, mod_names_stan_M)]
stan_simA_DM = stan_simA_DM[match(mod_names_DM, mod_names_stan_DM)]
# stan_simA_DM = stan_simA_DM[match(sapply(mod_names_DM, function(i){.x = strsplit(i, '_')[[1]]; return(paste0(c(paste0(.x[1:4], collapse="_"), '_', .x[6]), collapse=''))}),
#                                   mod_names_stan_DM)]
stan_simA_LNM = stan_simA_LNM[match(mod_names_LNM, mod_names_stan_LNM)]

sapply(list(stan_simA_M, stan_simA_DM, stan_simA_LNM), names)

idx = (1:length(stan_simA_M))[51]

give_plot_stan_vs_TMB_density = function(stan_list, TMB_list, nameout){
  pdf(paste0("../../results/assessing_models/", nameout), width = 10, height = 3)
  sapply(1:length(stan_list), function(idx){
    if(!is.na(stan_list[[idx]])){
      .xx = load_stan(paste0(folder_stan, stan_list[[idx]], ".RData"))
      post_betaslope = python_like_select_colnames(as.matrix(.xx$fit_stan), "beta\\[2,")
      ML_betaslope = select_slope_2(python_like_select_name(runs_M[[idx]]$par.fixed, "beta"))
      
      par(mfrow=c(1,ncol(post_betaslope)))
      for(i in 1:ncol(post_betaslope)){
        plot(density(post_betaslope[,i]))
        abline(v=ML_betaslope[i])
        # abline(v=ML_betaslope_2[i], col='red')
      }
    }
  })
  dev.off()
}

if(re_make_plots){
  give_plot_stan_vs_TMB_density(stan_list = stan_simA_M, TMB_list = runs_M, nameout = "simA_M_stan_TMB_betaslope.pdf")
  give_plot_stan_vs_TMB_density(stan_list = stan_simA_DM, TMB_list = runs_DM, nameout = "simA_DM_stan_TMB_betaslope.pdf")
  give_plot_stan_vs_TMB_density(stan_list = stan_simA_LNM, TMB_list = runs_LNM, nameout = "simA_LNM_stan_TMB_betaslope.pdf")
}

## scatterplot for betas for mean of posterior and ML from TMB

means_stan = function(stanfiles_list){
  lapply(stanfiles_list, function(stanfile){
  if(!is.na(stanfile)){
    stanfile = paste0(folder_stan, stanfile, '.RData')
    cat(stanfile, '\n')
    load(stanfile)
    betas = (python_like_select_colnames(as.matrix(fit_stan), 'beta\\[2'))
    if(is.null(dim(betas))){ ## one column
      beta_means=mean(betas)
    }else{
      beta_means=colMeans(betas)
    }
    return(beta_means)
  }else{return(NA)}
  })
}

means_stan_M = means_stan(stan_simA_M)
means_stan_DM = means_stan(stan_simA_DM)
means_stan_LNM = means_stan(stan_simA_LNM)
## add as many NAs as necessary for runs for which we don't have data
for(j in which(sapply(means_stan_M, typeof) == "logical")){ ## they are NA
  means_stan_M[[j]] = rep(NA,  (datasets[[j]]$d-1))
}
for(j in which(sapply(means_stan_DM, typeof) == "logical")){ ## they are NA
  means_stan_DM[[j]] = rep(NA,  (datasets[[j]]$d-1))
}
for(j in which(sapply(means_stan_LNM, typeof) == "logical")){ ## they are NA
  means_stan_LNM[[j]] = rep(NA,  (datasets[[j]]$d-1))
}


as.vector(which(table(df_beta_recovery$idx) != sapply(means_stan_DM, length)))
df_beta_recovery[df_beta_recovery$idx == 1,]
means_stan_DM[[1]]
## why do they all have 7 dimensions in the stan versions?
## I have removed the feature from the name
means_stan_DM


df_beta_recovery = cbind(df_beta_recovery,
                         means_stan_M = as.vector(unlist(means_stan_M)),
                         means_stan_DM = as.vector(unlist(means_stan_DM))
                         #means_stan_LNM = as.vector(unlist(means_stan_LNM))
                         )
df_beta_recovery$converged_M_TMB = sapply(gsub(".{1}$", "", rownames(df_beta_recovery)), function(i) 
  runs_M[[i]]$pdHess)
# df_beta_recovery$converged_LNM_TMB = sapply(gsub(".{1}$", "", rownames(df_beta_recovery)), function(i) 
#                                             runs_LNM[[i]]$pdHess)
df_beta_recovery$converged_DM_TMB = as.character(sapply(gsub(".{1}$", "", rownames(df_beta_recovery)), function(i) 
  runs_DM[[i]]$pdHess))

if(re_make_plots){
  ggplot(df_beta_recovery, aes(x=beta_est_M, means_stan_M))+geom_point()#+geom_abline(coef = c(0,1))
  ggplot(df_beta_recovery[df_beta_recovery$converged_M_TMB,], aes(x=beta_est_M, means_stan_M))+geom_point()#+geom_abline(coef = c(0,1))
  #ggplot(df_beta_recovery[!is.na(df_beta_recovery$means_stan_M),], aes(x=beta_est_M, means_stan_M))+geom_point()#+geom_abline(coef = c(0,1))
  ggsave("../../results/assessing_models/simA_M_scatter.pdf")
  ggplot(df_beta_recovery, aes(x=beta_est_DM, y=means_stan_DM))+
    geom_point()+geom_abline(coef = c(0,1))+
    guides(shape=FALSE)
  ggsave("../../results/assessing_models/simA_DM_scatter.pdf")
  ggplot(df_beta_recovery, aes(x=beta_est_DM, y=means_stan_DM, col=converged_DM_TMB))+
    geom_point()+geom_abline(coef = c(0,1))+
    guides(shape=FALSE)
  ggsave("../../results/assessing_models/simA_DM_scatter_col.pdf")
  ggplot(df_beta_recovery[df_beta_recovery$converged_DM_TMB == 'TRUE',], aes(x=beta_est_DM, y=means_stan_DM))+
    geom_point()+geom_abline(coef = c(0,1))+
    guides(shape=FALSE)
  ggsave("../../results/assessing_models/simA_DM_scatter_col_nodivergent.pdf", width = 5, height = 5)
  ggplot(df_beta_recovery %>% filter(!is.na(means_stan_LNM)),
         aes(x=(beta_est_LNM), y=means_stan_LNM, col=converged_LNM_TMB))+
    geom_point()+
    geom_abline(coef = c(0,1))
  ggsave("../../results/assessing_models/simA_LNM_scatter.pdf")
}
## even in cases where LNM hasn't converged (which is many) there seems to be an okay correlation

#-------------------------------------------------------------------------------------------------#

## Are zeros recovered with Stan?
ggplot(droplevels(df_beta_recovery) %>% filter(!is.na(means_stan_M)),
       aes(x=(beta_true), y=(means_stan_M), col=n))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~beta_gamma_shape, scales = "free_x", ncol=5)
ggsave(paste0("../../results/TMB/", dataset_generation, "_stan_recoveries_betashape_facet.pdf"), height = 3, width = 11)

ggplot(droplevels(df_beta_recovery) %>% filter(!is.na(means_stan_M)),
       aes(x=(beta_true), y=(means_stan_M), col=n))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~beta_gamma_shape, scales = "free", ncol=5)
ggsave(paste0("../../results/TMB/", dataset_generation, "_stan_recoveries_betashape_facet_freescales.pdf"), height = 3, width = 11)

## Looking at the posteriors of these betas that are not well recovered
## find the datasets where there are true zeros
stan_simA_M
with_zeros_name = sapply(rownames(df_beta_recovery[df_beta_recovery$beta_true == 0,]), function(i){
  .x = strsplit(i, '_')[[1]]
  return(paste0(.x[1:(length(.x)-1)], collapse = '_'))
})

posteriors_truezeros = lapply(paste0(folder_stan, python_like_select(stan_simA_M, 'beta_intensity0_'), ".RData"), load_stan)

ggplot(melt(lapply(posteriors_truezeros, function(i) python_like_select_colnames(as(i$fit_stan, 'matrix'), 'beta\\[2,'))),
       aes(x=value, col=interaction(parameters, L1), group=interaction(parameters, L1)))+geom_density()

library(ggridges)

give_names = function(vec, names){
  names(vec) = names
  vec
}

give_names_beta2 = function(vec, names){
  names(vec) =  paste0('beta[2,', 1:length(vec), ']')
  vec
}
names_as_row = function(i){
  .x = rbind.data.frame(names(i), i)
  names(.x) = NULL
  .x
}

ML_est = (melt(lapply(paste0(sapply(modify_M_name(python_like_select(stan_simA_M, 'beta_intensity0_')), modify_M_name_2), '_'),
       function(i) t((give_names_beta2(df_beta_recovery[grepl(i, rownames(df_beta_recovery)),'beta_est_M']) )))))
colnames(ML_est) = c('iterations', 'parameters', 'value', 'L1')
head(melt(lapply(posteriors_truezeros, function(i) python_like_select_colnames(as(i$fit_stan, 'matrix'), 'beta\\[2,'))))
ML_est$group=apply(ML_est[,c('parameters', 'L1')], 1, paste0, collapse='.')

posterior_est = melt(lapply(posteriors_truezeros, function(i) python_like_select_colnames(as(i$fit_stan, 'matrix'), 'beta\\[2,')))
posterior_est$group=apply(posterior_est[,c('parameters', 'L1')], 1, paste0, collapse='.')
posterior_est$group = factor(posterior_est$group, levels=apply(ML_est[order(ML_est$value),c('parameters', 'L1')], 1, paste0, collapse='.'))
ML_est$group = factor(ML_est$group, levels=apply(ML_est[order(ML_est$value),c('parameters', 'L1')], 1, paste0, collapse='.'))

ggplot(posterior_est,
       aes(x=value, y=1,#interaction(parameters, L1),
           group=group, fill=L1)#, fill=interaction(parameters, L1))
       )+
  geom_density_ridges()+
  facet_wrap(.~group, ncol=4)+
  geom_vline(data = ML_est, aes(xintercept=value, group=group),
             #, col=interaction(parameters, L1)))#+facet_wrap(.~interaction(L1,Var2))
             col='red')+
  geom_vline(aes(xintercept=0),
             #, col=interaction(parameters, L1)))#+facet_wrap(.~interaction(L1,Var2))
             col='gray', lty='dashed')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
ggsave(paste0("../../results/assessing_models/zero_recovery_", dataset_generation, "_stan_recoveries_truezero.png"), height = 2.5, width = 5)

## another view of the same data
ggplot(posterior_est,
       aes(x=value, y=parameters,
           group=interaction(parameters, L1), fill=L1)#, fill=interaction(parameters, L1))
)+
  geom_density_ridges()+
  facet_wrap(.~L1, ncol=4)+
  geom_vline(aes(xintercept=0),
             #, col=interaction(parameters, L1)))#+facet_wrap(.~interaction(L1,Var2))
             col='gray', lty='dashed')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

## Trying to find how we could avoid this inflation in type I errors
## looking at the correlation betweeen estimates, to see if it could be a cause

## subsample the data
sapply(1:length(posteriors_truezeros), function(idx_posteriors_truezeros){
  pairs(base::subset(python_like_select_colnames(as(posteriors_truezeros[[idx_posteriors_truezeros]]$fit_stan, 'matrix'), 'beta\\[2'),
               sample(c(T,F), size = 2000, replace = T, prob = c(0.1, 0.9))))
})

## is the error in zero coefficients larger than for other coefficients?
## plot in which in the x axis are the sorted coefficients, and in the y axis the errors (|est-true|)

pdf(paste0("../../results/TMB/", dataset_generation, "_ordered_betas_and_error.pdf"), height = 3, width = 5)
plot((df_beta_recovery$beta_est_M - df_beta_recovery$beta_true)[order(df_beta_recovery$beta_true, decreasing = F)], type='h',
     xlab='Rank of betas', ylab='Error of betas')
dev.off()

pdf(paste0("../../results/TMB/", dataset_generation, "_ordered_betas_and_error_zoomed_in.pdf"), height = 3, width = 5)
plot((df_beta_recovery$beta_est_M - df_beta_recovery$beta_true)[order(df_beta_recovery$beta_true, decreasing = F)], type='h',
     xlab='Rank of betas', ylab='Error of betas', ylim=c(-2,2))
dev.off()

#-------------------------------------------------------------------------------------------------#

## Two M runs for some reason
runs_M1 = readRDS("../../data/robjects_cache/tmb_results_simulations/20200625_simulation_runs_M_5c268769-daa7-4418-a41b-08c6ed9cf5d3.RDS")
runs_M2 = readRDS("../../data/robjects_cache/tmb_results_simulations/20200625_simulation_runs_M_970f6533-5025-4eef-a598-956422846ddf.RDS")

all(sapply(sapply(runs_M1, function(i) i$par.fixed), length) == sapply(sapply(runs_M2, function(i) i$par.fixed), length))

plot(sapply(runs_M1, function(i) i$par.fixed) %>% unlist,
     sapply(runs_M2, function(i) i$par.fixed) %>% unlist)
