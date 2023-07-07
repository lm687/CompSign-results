
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

multiple_runs = T
generation = "GenerationJnormTwoLambdasOneChangingBeta"
generation = "GenerationJnormBTwoLambdasOneChangingBeta"

source("../../../2_inference_TMB/helper_TMB.R")
source("../../../1_create_ROO/roo_functions.R")
source("helper_model_assessment.R")

library(grid)
library(gridExtra)
library(reshape2)
library(jcolors)
library(cowplot)
library(ggrepel)
library(dplyr)

if(multiple_runs){
  flder_out <- "../../../../results/results_TMB/simulated_datasets/mixed_effects_models_multiple/"
  flder_in <- "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries_multiple/"
}else{
  flder_out <- "../../../../results/results_TMB/simulated_datasets/mixed_effects_models/"
  flder_in <- "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries/"
}

runs_fullREM0 = readRDS(paste0(flder_in, generation, "_fullREM.RDS"))
runs_fullREDMSL0 = readRDS(paste0(flder_in, generation, "_fullREDMsinglelambda.RDS"))
runs_diagREDMSL0 = readRDS(paste0(flder_in, generation, "_diagREDMsinglelambda.RDS"))
runs_diagREDM0 = readRDS(paste0(flder_in, generation, "_diagREDM.RDS"))

## match them all (wrt fullREM)
runs_fullREDMSL0 <- runs_fullREDMSL0[match(rownames(runs_fullREM0), rownames(runs_fullREDMSL0)),]
runs_diagREDMSL0 <- runs_diagREDMSL0[match(rownames(runs_fullREM0), rownames(runs_diagREDMSL0)),]
runs_diagREDM0 <- runs_diagREDM0[match(rownames(runs_fullREM0), rownames(runs_diagREDM0)),]

## Problem with convergence is acute in DM
table(is.na(runs_fullREM0$beta_est))
table(is.na(runs_fullREDMSL0$beta_est))
table(is.na(runs_diagREDMSL0$beta_est))
table(is.na(runs_diagREDM0$beta_est))

table((runs_fullREM0$converged))
table((runs_fullREDMSL0$converged))
table((runs_diagREDMSL0$converged))
table((runs_diagREDM0$converged))

barplot(c(runs_fullREM0=sum((runs_fullREM0$converged)),
          runs_fullREDMSL0=sum((runs_fullREDMSL0$converged)),
          runs_diagREDMSL0=sum((runs_diagREDMSL0$converged)),
          runs_diagREDM0=sum((runs_diagREDM0$converged))))
image(cbind(runs_fullREM0=((runs_fullREM0$converged)),
            runs_fullREDMSL0=((runs_fullREDMSL0$converged)),
            runs_diagREDMSL0=((runs_diagREDMSL0$converged)),
            runs_diagREDM0=((runs_diagREDM0$converged))))

runs_fullREM <- runs_fullREM0[runs_fullREM0$converged,]
runs_fullREDMSL <- runs_fullREDMSL0[runs_fullREDMSL0$converged,]
runs_diagREDMSL <- runs_diagREDMSL0[runs_diagREDMSL0$converged,]
runs_diagREDM <- runs_diagREDM0[runs_diagREDM0$converged,]

joint_df = cbind.data.frame(fullRE_M=runs_fullREM,
                            fullRE_DMSL=runs_fullREDMSL[match(rownames(runs_fullREM),
                                                              rownames(runs_fullREDMSL)),],
                            diagRE_DMSL=runs_diagREDMSL[match(rownames(runs_fullREM),
                                                              rownames(runs_diagREDMSL)),],
                            diagRE_DM=runs_diagREDM[match(rownames(runs_fullREM),
                                                          rownames(runs_diagREDM)),])
joint_df$names <- sapply(1:nrow(joint_df), function(i) gsub(paste0(joint_df[i,]$fullRE_M.idx_within_dataset, "$"), "", rownames(joint_df)[i]))

sort(unique(gsub("\\..*","",rownames(joint_df))))

## Read in datasets
cat('Reading datasets\n')
datasets_files = list.files("../../../../data/assessing_models_simulation/datasets/", full.names = TRUE)
datasets_files = datasets_files[grep(pattern = paste0('/multiple_', generation, '_'), datasets_files)]
length(datasets_files)

datasets_files[grepl('multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_NA_NA_NA_dataset', datasets_files)]

datasets_files = datasets_files[match(unique(joint_df$names),
                                      gsub(".RDS", "", basename(datasets_files)))]

datasets_files[grepl('multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_NA_NA_NA_dataset', datasets_files)]

length(datasets_files)
length(unique(joint_df$fullRE_M.idx_within_dataset))
datasets = lapply(datasets_files, readRDS)

if((generation %in% c("GenerationMixturePCAWG", "GenerationMixturefewersignaturesPCAWG", "GenerationMixturefewersignaturespairedPCAWG")) | grepl('GenerationMixturefewersignaturespaired', generation)  | grepl('signaturespaired', generation) ){
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
  if(grepl('Mixture', generation)){
    stop('Are you sure you are using probabilities beta_gamma_shape for and not log-ratios?\n')
  }
}
names(datasets) = (gsub("_dataset.RDS", "", basename(datasets_files)))

DA_bool = ( sapply(datasets, function(i) i$beta_gamma_shape) > 0 )

sapply(gsub(".RDS", "", names(DA_bool)[1:3]), function(i) grep(i, rownames(joint_df)))

length(datasets)
dim(joint_df)

runs_fullREM

datasets[[1]]$beta_gamma_shape
rownames(joint_df[1,])

joint_df$diagRE_DM.idx
datasets[[1]]$d

joint_df[match(unique(joint_df$diagRE_DM.idx), joint_df$diagRE_DM.idx),]$diagRE_DM.pvals_adj < 0.05

minimalperts <- lapply(unique(joint_df$fullRE_M.idx), function(i){
  if(is.na(i)){
    list(betas_perturbed=NA)
  }else{
    joint_df_subset <- joint_df[which(joint_df$diagRE_DM.idx == i),]
    if(nrow(joint_df_subset) == 0 ){
      list(betas_perturbed=NA)
    }else{
      joint_df[which(joint_df$diagRE_DM.idx == i), 'diagRE_DM.beta_est']
      # joint_df[which(joint_df$diagRE_DM.idx == i),'fullRE_M.beta_est']
      if(any(is.na(joint_df_subset[,'diagRE_DM.beta_stderr']))){
        list(betas_perturbed=NA)
      }else{
        give_min_pert(idx_sp = NA, list_runs = NA, logR_names_vec = NA,
                      df_betas = joint_df_subset[,c('diagRE_DM.beta_est', 'diagRE_DM.beta_stderr')])
      }
      }
  }
})

## minimal perturbation accuracy
minimalperts_signifbool <- sapply(sapply(minimalperts, `[`, 'betas_perturbed'), function(i) !(all(i == "FALSE")))
minimalperts_signifbool_in_DA <- sapply(sapply(minimalperts, `[`, 'betas_perturbed'), function(i) !(i[1] == "FALSE"))
names(minimalperts) <- basename(datasets_files)
names(minimalperts_signifbool_in_DA) <- basename(datasets_files)
minimalperts_signif <- ifelse(minimalperts_signifbool, 0.001, 1)
minimalperts_signif_in_DA <- ifelse(minimalperts_signifbool_in_DA, 0.001, 1)

length(minimalperts_signif)
length(sapply(datasets, `[`, 'beta_gamma_shape') != 0)

varying_n_all_converged <- give_accuracies_with_varying_var(var = 'n', datasets_arg = datasets, single_pval = T,
                                    pvals_data_frame_arg = data.frame(test=minimalperts_signif,
                                                                      true=sapply(datasets, `[`, 'beta_gamma_shape') != 0))
varying_netagamma_all_converged <- give_accuracies_with_varying_var(var = 'beta_gamma_shape', datasets_arg = datasets, single_pval = T,
                                                            pvals_data_frame_arg = data.frame(test=minimalperts_signif,
                                                                                              true=sapply(datasets, `[`, 'beta_gamma_shape') != 0))
varying_n_all_converged_in_DA <- give_accuracies_with_varying_var(var = 'n', datasets_arg = datasets, single_pval = T,
                                                            pvals_data_frame_arg = data.frame(test=minimalperts_signif_in_DA,
                                                                                              true=sapply(datasets, `[`, 'beta_gamma_shape') != 0))

if(F){
  .xxxidx <- which(grepl('multiple_GenerationJnormBTwoLambdasOneChangingBeta_20_100_80_6_0.01_NA_NA_NA_dataset1', names(datasets)))
  names(datasets)[.xxxidx]
  datasets[[.xxxidx]]$beta
  runs_diagREDM[grepl("multiple_GenerationJnormBTwoLambdasOneChangingBeta_20_100_80_6_0.01_NA_NA_NA_dataset11", rownames(runs_diagREDM)),]
  
  a1 <- readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnormBTwoLambdasOneChangingBeta_20_100_80_6_0.01_diagREDM_NA_NA_NA_dataset1.RDS"))
  plot_betas(a1) ## this is differentially abundant but not because of the actually DA run
  createBarplot(give_dummy_rownames(give_dummy_colnames(normalise_rw(datasets[[.xxxidx]]$W))))
  give_min_pert(idx_sp = 1, list_runs = list(a1), logR_names_vec = list(vector_cats_to_logR(paste0('c', 1:6))))
  wald_TMB_wrapper(a1)
  
  # datasets[[.xxxidx]]$
  dataset_xxx = load_PCAWG(ct = datasets_files[.xxxidx], typedata = 'simulation', simulation = T,
                       path_to_data = NA, read_directly=T)
  resort_columns <- function(i, order){
    i$Y = i$Y[,order]
  }
  dataset_xxx_bs1 <- resort_columns(dataset_xxx, 1)
  
}

ggplot(varying_n_all_converged, aes(y=Accuracy, x=n))+geom_line()+geom_point()+theme_bw()
ggplot(varying_n_all_converged, aes(y=FPR, x=n))+geom_line()
ggplot(varying_netagamma_all_converged, aes(y=Sensitivity, x=factor(beta_gamma_shape)))+geom_line(aes(group=1))+theme_bw()+
  labs(x='Betta gamma')+geom_point()
ggsave(paste0("../../../../results/results_TMB/pcawg/assessing_models/minimal_perturbation/Sensitivity_", generation, '.pdf'),
       width=3, height=3)

ggplot(varying_n_all_converged_in_DA, aes(y=Accuracy, x=n))+geom_line()+geom_point()+theme_bw()+
  geom_line(data = varying_n_all_converged,  aes(y=Accuracy, x=n), col='blue')+
  geom_point(data = varying_n_all_converged,  aes(y=Accuracy, x=n), col='blue')

joint_df$fullRE_M.idx 
duplicated(names(datasets))




length(minimalperts_signif)
length(datasets)
length(datasets_files)
length(unique(joint_df$fullRE_M.idx_within_dataset))
length(unique(unique(joint_df$diagRE_DM.idx)))

minimalperts_signifbool
plot(joint_df$diagRE_DMSL.beta_intercept_true,
     joint_df$diagRE_DMSL.beta_true)
plot(joint_df$diagRE_DMSL.beta_intercept_true,
     scale(joint_df$diagRE_DMSL.beta_true, scale = T, center = T))

minimal_perturbation_df <- cbind.data.frame(diagRE_DMSL.beta_intercept_true=joint_df$diagRE_DMSL.beta_intercept_true,
                 diagRE_DMSL.beta_true=joint_df$diagRE_DMSL.beta_true,
                 diagRE_DMSL.beta_intercept_est=joint_df$diagRE_DMSL.beta_intercept_est,
                 diagRE_DMSL.beta_est=joint_df$diagRE_DMSL.beta_est,
                 minimalperts_signifbool=minimalperts_signifbool[as.numeric(factor(joint_df$fullRE_M.idx))],
                 d=joint_df$fullRE_M.d,
                 n=joint_df$fullRE_M.n,
                 name=joint_df$names)


ggplot(minimal_perturbation_df[minimal_perturbation_df$diagRE_DMSL.beta_true != 0,],
       aes(x=diagRE_DMSL.beta_intercept_true, y=diagRE_DMSL.beta_true, col=minimalperts_signifbool))+
  geom_point()+theme_bw()+theme(legend.position = "bottom")+facet_wrap(.~d)
ggplot(minimal_perturbation_df[minimal_perturbation_df$diagRE_DMSL.beta_true != 0,],
       aes(x=diagRE_DMSL.beta_intercept_true, y=diagRE_DMSL.beta_true, col=minimalperts_signifbool))+
  geom_point()+theme_bw()+theme(legend.position = "bottom")+facet_wrap(.~n)


d4 <- minimal_perturbation_df %>% dplyr::filter(d==4)
(d4[order(d4$diagRE_DMSL.beta_true, decreasing = T)[1],])
(d4[order(d4$diagRE_DMSL.beta_true, decreasing = T)[2],])

datasets['multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_NA_NA_NA_dataset0.RDS']
datasets['multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_NA_NA_NA_dataset1.RDS']
data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnormTwoLambdasOneChangingBeta_100_100_80_5_0_fullREDMsinglelambda_NA_NA_NA_dataset0.RDS


which(basename(datasets_files) == 'multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_NA_NA_NA_dataset0.RDS')
which(basename(datasets_files) == 'multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_NA_NA_NA_dataset1.RDS')
minimalperts$multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_NA_NA_NA_dataset0.RDS
minimalperts$multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_NA_NA_NA_dataset1.RDS

datasets_files[grepl('multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_NA_NA_NA_dataset', datasets_files)]

a1 <- readRDS("../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_diagREDM_NA_NA_NA_dataset0.RDS")
a2 <- readRDS("../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnormTwoLambdasOneChangingBeta_20_100_80_4_4_diagREDM_NA_NA_NA_dataset1.RDS")

grid.arrange(plot_betas((a1)), plot_betas((a2)))

give_min_pert(idx_sp = NA, list_runs = NA, logR_names_vec = NA, df_betas = python_like_select_rownames(summary(a1), 'beta')[c(F,T),])
give_min_pert(idx_sp = NA, list_runs = NA, logR_names_vec = NA, df_betas = python_like_select_rownames(summary(a2), 'beta')[c(F,T),])

ggplot(minimal_perturbation_df[minimal_perturbation_df$diagRE_DMSL.beta_true != 0,],
       aes(x=diagRE_DMSL.beta_intercept_est, y=diagRE_DMSL.beta_est, col=minimalperts_signifbool))+
  geom_point()+theme_bw()+theme(legend.position = "bottom")

only_logR_with_change <- cbind(joint_df[,grepl('diagRE_DM.', colnames(joint_df))], names=joint_df$names)
colnames(only_logR_with_change) <- gsub("diagRE_DM.", "", colnames(only_logR_with_change))
only_logR_with_change <- only_logR_with_change[only_logR_with_change$idx_within_dataset == 1,]

ggplot(data = only_logR_with_change,
       aes(x=factor(beta_true), y=beta_intercept_true, col=pvals_adj))+geom_point()+theme_bw()
