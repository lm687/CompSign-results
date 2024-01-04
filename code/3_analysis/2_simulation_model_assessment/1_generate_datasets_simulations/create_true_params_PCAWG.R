## create true params based on PCAWG data
## ----------------------------------------------------------------------------------------- ##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##
source("../../../2_inference_TMB/helper_TMB.R")
flder_save <- '../../../../data/assessing_models_simulation/additional_files/'
system(paste0('mkdir -p ', flder_save))
## ----------------------------------------------------------------------------------------- ##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../../../data/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
enough_samples

read_info <- function(ct){
  .x <- list(diagRE_DMDL_SP = try(readRDS(paste0("../../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS"))),
             fullRE_DMDL_SP = try(readRDS(paste0("../../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS"))),
             dataset_active_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../../data/", override_warning_X_Z = T))
  
  .x
}
read_info_list <- lapply(enough_samples, function(ct){
  read_info(ct)
}); names(read_info_list) <- enough_samples

fullDM_extra <- list()
##-----------------------------------------------------------------------------------------------------##

#---------------------------------------------------------------------------#
TMB::compile("../../../2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("../../../2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial"))
TMB::compile("../../../2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("../../../2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial"))

## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##

read_diagRE <- function(ct){
  try(readRDS(paste0("../../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS")))
}

read_fullRE <- function(ct){
  try(readRDS(paste0("../../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS")))
}

diagRE_params <- lapply(enough_samples, read_diagRE)
fullRE_params <- lapply(enough_samples, read_fullRE)
names(diagRE_params) <- names(fullRE_params) <- enough_samples
## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##
ct = 'Kidney-ChRCC'
betas1 <- matrix(python_like_select_name(diagRE_params[[ct]]$par.fixed, 'beta'), nrow=2)
cov1 <- diag(exp(python_like_select_name(diagRE_params[[ct]]$par.fixed, 'logs_sd_RE')))
fullDM_extra[[ct]] <- wrapper_run_TMB(read_info_list[[ct]]$dataset_active_sigs,
                                      model = "fullRE_DM", use_nlminb=T, smart_init_vals=T)
fullDM_extra[[ct]]
# cov1_full <- L_to_cov(python_like_select_name(fullRE_params[[ct]]$par.fixed, 'cov_par_RE'), d = ncol(betas1))
lambdas <- python_like_select_name(diagRE_params[[ct]]$par.fixed, 'log_lambda')
## lambdas: we don't save them but add them in the Snakemake code
exp(lambdas)
saveRDS(cov1, paste0(flder_save,"multiple_fixed_covmatPCAWG1.RDS"))
cov1 <- readRDS(paste0(flder_save,"multiple_fixed_covmatPCAWG1.RDS"))
saveRDS(betas1[1,], paste0(flder_save,"multiple_fixed_betainterceptPCAWG1.RDS"))
saveRDS(betas1[2,], paste0(flder_save,"multiple_fixed_betaslopePCAWG1.RDS"))
## 7 log-ratios
median(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) ## median num of mutations
## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##
ct <- 'CNS-GBM'
fullDM_extra <- list()
betas2 <- matrix(python_like_select_name(diagRE_params$`CNS-GBM`$par.fixed, 'beta'), nrow=2)
cov2 <- diag(exp(python_like_select_name(diagRE_params$`CNS-GBM`$par.fixed, 'logs_sd_RE')))
fullDM_extra[[ct]] <- wrapper_run_TMB(read_info_list[[ct]]$dataset_active_sigs,
                model = "fullRE_DM", use_nlminb=T, smart_init_vals=T)
fullDM_extra[[ct]]
## save run that worked!
saveRDS(object = fullDM_extra[[ct]],
        file = paste0("../../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS"))
cov2
fullDM_extra[[ct]]
cov2_full <- L_to_cov(python_like_select_name(fullDM_extra[[ct]]$par.fixed, 'cov_par_RE'), d = ncol(betas2))
diag(cov2_full) <- (exp(python_like_select_name(fullDM_extra[[ct]]$par.fixed, 'logs_sd_RE')))

cov2_full ## cov (full) for full
cov2 ## cov (diag) for diag

python_like_select_name(diagRE_params[[ct]]$par.fixed, 'log_lambda')
python_like_select_name(fullDM_extra[[ct]]$par.fixed, 'log_lambda')

python_like_select_name(diagRE_params[[ct]]$par.fixed, 'beta')
python_like_select_name(fullDM_extra[[ct]]$par.fixed, 'beta')

par(mfrow=c(1,1))
plot(unlist(cov2_full),
     unlist(cov(matrix(fullDM_extra[[ct]]$par.random, ncol=ncol(betas2)))))
abline(coef=c(0,1))

# xxx <- readRDS(paste0("../../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS"))
# xxx <- readRDS("../../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_CNS-GBM_signaturesMSE.RDS")
lambdas <- python_like_select_name(diagRE_params$`CNS-GBM`$par.fixed, 'log_lambda')
## lambdas: we don't save them but add them in the Snakemake code
exp(lambdas) ## from diag
mean(exp(python_like_select_name(fullDM_extra[[ct]]$par.fixed, 'log_lambda'))) ## from full. using 28

saveRDS(cov2, paste0(flder_save,"multiple_fixed_covmatPCAWG2.RDS"))
saveRDS(cov2full,  paste0(flder_save, "multiple_fixed_covmatFULLPCAWG2.RDS"))
saveRDS(betas2[1,],  paste0(flder_save,"multiple_fixed_betainterceptPCAWG2.RDS"))
saveRDS(betas2[2,],  paste0(flder_save,"multiple_fixed_betaslopePCAWG2.RDS"))
median(rowSums(read_info_list$`CNS-GBM`$dataset_active_sigs$Y))
# [1] 3401
read_info_list$`CNS-GBM`$diagRE_DMDL_SP
## log_lambda  -3.38334738   0.17484245
# log_lambda  -4.73237424   0.17066885
exp(-4)*1000 ## for diagRE, where lambdas are scaled
saveRDS(diag(1,5),  paste0(flder_save,"multiple_fixed_covdiagd6.RDS"))
## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##
sapply(read_info_list, function(i) dim(i$dataset_active_sigs$Y))
sapply(read_info_list, function(i) mean(rowSums(i$dataset_active_sigs$Y)))
## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##
ct <- 'Lung-SCC'
ct <- 'Prost-AdenoCA'

betas2 <- matrix(python_like_select_name(diagRE_params[[ct]]$par.fixed, 'beta'), nrow=2)
cov2 <- diag(exp(python_like_select_name(diagRE_params[[ct]]$par.fixed, 'logs_sd_RE')))
fullRE_params[[ct]]$par.fixed
cov2_full <- L_to_cov(python_like_select_name(fullRE_params[[ct]]$par.fixed, 'cov_par_RE'), d = 13)
diag(cov2_full) <- (exp(python_like_select_name(fullRE_params[[ct]]$par.fixed, 'logs_sd_RE')))
image(cov2_full)
image(cov(matrix(fullRE_params[[ct]]$par.random, ncol=13)), main='cov from U')
plot(as.vector(cov2_full),
as.vector(cov(matrix(fullRE_params[[ct]]$par.random, ncol=13))), col=
  as.factor(as.vector(diag(1, 13))))
abline(coef = c(0,1), lty='dashed') ## this is NOT BAD AT ALL!
lambdas <- python_like_select_name(diagRE_params[[ct]]$par.fixed, 'log_lambda')
# ## lambdas: we don't save them but add them in the Snakemake code
exp(lambdas)
saveRDS(cov2_full,  paste0(flder_save, "multiple_fixed_covmatPCAWG3.RDS"))
saveRDS(betas2[1,], paste0(flder_save, "multiple_fixed_betainterceptPCAWG3.RDS"))
saveRDS(betas2[2,], paste0(flder_save, "multiple_fixed_betaslopePCAWG3.RDS"))
median(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y))
# # [1] 1516
# read_info_list[[ct]]$diagRE_DMDL_SP
# ## log_lambda  -3.38334738   0.17484245
# # log_lambda  -4.73237424   0.17066885
mean(exp(lambdas))*1000 ## for diagRE, where lambdas are scaled

dim(read_info_list[[ct]]$dataset_active_sigs$Y)[1]/2
# [1] 208 ## num samples
## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##
ct <- 'Lung-SCC'

betas2 <- matrix(python_like_select_name(diagRE_params[[ct]]$par.fixed, 'beta'), nrow=2)
cov2 <- diag(exp(python_like_select_name(diagRE_params[[ct]]$par.fixed, 'logs_sd_RE')))
plot(as.vector(cov2),
     as.vector(cov(matrix(diagRE_params[[ct]]$par.random, ncol=ncol(betas2)))), col=
       as.factor(as.vector(diag(1, ncol(betas2)))))
abline(coef = c(0,1), lty='dashed') ## this is NOT BAD AT ALL!
lambdas <- python_like_select_name(diagRE_params[[ct]]$par.fixed, 'log_lambda')
# ## lambdas: we don't save them but add them in the Snakemake code
exp(lambdas)
saveRDS(cov2,  paste0(flder_save, "multiple_fixed_covmatPCAWG4.RDS"))
saveRDS(betas2[1,], paste0(flder_save, "multiple_fixed_betainterceptPCAWG4.RDS"))
saveRDS(betas2[2,], paste0(flder_save, "multiple_fixed_betaslopePCAWG4.RDS"))
saveRDS(rep(0, ncol(betas2)), paste0(flder_save, "multiple_fixed_betaslopezerosPCAWG4.RDS"))
saveRDS(rep(0.1, ncol(betas2)), paste0(flder_save, "multiple_fixed_betaslopeonechangingPCAWG4.RDS")) ## only baseline signature changes
median(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y))
# 14072.5
mean(exp(lambdas))*1000 ## for diagRE, where lambdas are scaled

dim(read_info_list[[ct]]$dataset_active_sigs$Y)[1]/2 ## number of samples

## now with fullRE params
fullDM_extra[[ct]] <- wrapper_run_TMB(read_info_list[[ct]]$dataset_active_sigs,
                                      model = "fullRE_DM", use_nlminb=T, smart_init_vals=T)
saveRDS(object = fullDM_extra[[ct]],
        file = paste0("../../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS"))
cov2_full <- L_to_cov(python_like_select_name(fullDM_extra[[ct]]$par.fixed, 'cov_par_RE'), d = ncol(betas2))
diag(cov2_full) <- (exp(python_like_select_name(fullDM_extra[[ct]]$par.fixed, 'logs_sd_RE')))
cov2_full
saveRDS(cov2_full,  paste0(flder_save, "multiple_fixed_covmatFULLPCAWG4.RDS"))
mean(exp(python_like_select_name(fullDM_extra[[ct]]$par.fixed, 'log_lambda'))) ## for fullRE
## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##
sapply(read_info_list, function(i) ncol(i$dataset_active_sigs$Y))
cbind(sapply(read_info_list, function(i) ncol(i$dataset_active_sigs$Y)),
      sapply(read_info_list, function(i) nrow(i$dataset_active_sigs$Y)/2))
## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##
## PCAWG5.RDS are deprecated: they were the ones below but without lambda having been scaled, i.e. it had not been *1000
## moreover, before the betas had been taken from diagREDM, but now all the parameters are taken from fullREDM, 
## where no scaling is needed
ct = 'Eso-AdenoCA'


# fullDM_extra[[ct]] <- wrapper_run_TMB(read_info_list[[ct]]$dataset_active_sigs,
#                                       model = "fullRE_DM", use_nlminb=T, smart_init_vals=T)
# fullRE_params[[ct]]
# fullDM_extra[[ct]]
# saveRDS(object = fullDM_extra[[ct]],
#         file = paste0("../../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS"))
# fullDM_extra[[ct]] <- readRDS(file = paste0("../../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS"))

exp(python_like_select_name(diagRE_params[[ct]]$par.fixed, 'log_lambda'))*1000; exp(python_like_select_name(fullRE_params[[ct]]$par.fixed, 'log_lambda'))
plot(python_like_select_name(diagRE_params[[ct]]$par.fixed, 'beta'), python_like_select_name(fullRE_params[[ct]]$par.fixed, 'beta')); abline(coef=c(0,1))

betas5 <- matrix(python_like_select_name(fullRE_params[[ct]]$par.fixed, 'beta'), nrow=2)
# cov1 <- diag(exp(python_like_select_name(diagRE_params[[ct]]$par.fixed, 'logs_sd_RE')))
cov5_full <- L_to_cov(python_like_select_name(fullRE_params[[ct]]$par.fixed, 'cov_par_RE'), d = ncol(betas5))
plot(as.vector(cov5_full),
     as.vector(cov(matrix(fullRE_params[[ct]]$par.random, ncol=ncol(betas5))))); abline(coef = c(0,1))
## in diagRE to get lambda one must multiply by 1000, i.e. exp(log_lambda)*1000
lambdas5 <- python_like_select_name(fullRE_params[[ct]]$par.fixed, 'log_lambda')
## lambdas: we don't save them but add them in the Snakemake code
exp(lambdas5)
saveRDS(exp(lambdas5), paste0(flder_save,"multiple_fixed_lambdaPCAWG6.RDS"))
saveRDS(exp(lambdas5)/2, paste0(flder_save,"multiple_fixed_lowlambdaPCAWG6.RDS"))
saveRDS(exp(lambdas5)*2, paste0(flder_save,"multiple_fixed_highlambdaPCAWG6.RDS"))
saveRDS(cov5_full, paste0(flder_save,"multiple_fixed_covmatFULLPCAWG6.RDS"))
# cov5 <- readRDS(paste0(flder_save,"multiple_fixed_covmatPCAWG6.RDS"))
saveRDS(betas5[1,], paste0(flder_save,"multiple_fixed_betainterceptPCAWG6.RDS"))
saveRDS(betas5[2,], paste0(flder_save,"multiple_fixed_betaslopePCAWG6.RDS"))
saveRDS(ncol(betas5)+1, paste0(flder_save,"multiple_fixed_dPCAWG6.RDS"))
saveRDS(nrow(read_info_list[[ct]]$dataset_active_sigs$Y)/2, paste0(flder_save,"multiple_fixed_nPCAWG6.RDS"))
## 9 log-ratios

nlambdas5 = lapply(split_matrix_in_half((read_info_list[[ct]]$dataset_active_sigs$Y)), rowSums)
median(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) ## median num of mutations
saveRDS(nlambdas5, paste0(flder_save,"multiple_fixed_nlambdaPCAWG6.RDS"))
saveRDS(lapply(lapply(nlambdas5, '/', 2), round), paste0(flder_save,"multiple_fixed_lownlambdaPCAWG6.RDS"))
saveRDS(lapply(lapply(nlambdas5, '/', 10), round), paste0(flder_save,"multiple_fixed_low2nlambdaPCAWG6.RDS"))
saveRDS(lapply(lapply(nlambdas5, '/', 100), round), paste0(flder_save,"multiple_fixed_low3nlambdaPCAWG6.RDS"))
## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##

## ----------------------------------------------------------------------------------------- ##
