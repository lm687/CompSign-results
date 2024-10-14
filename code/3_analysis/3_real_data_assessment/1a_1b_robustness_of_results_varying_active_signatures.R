##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
source("helper.R")

library(TMB)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(dplyr)
library(jcolors)
library(viridis)
library(mutSigExtractor)
library(reshape2)
library(CompSign)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
##' Including all signatures, even with artefacts (SBS_SIGNATURE_PROFILES_V2 and SBS_SIGNATURE_PROFILES_V3 from
##' mutSigExtractor do not include artefactual signatures)
sigs_cosmic0 <- read.table(paste0( "../../../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
                           stringsAsFactors = FALSE, sep = ',', header = TRUE)
rownames(sigs_cosmic0) <- paste0(substr(sigs_cosmic0$SubType, 1, 1),'[',
                                 sigs_cosmic0$Type, ']', substr(sigs_cosmic0$SubType, 3, 3))
sigs_cosmic0 <- sigs_cosmic0[-c(1,2)];
sigs_cosmic <- colnames(sigs_cosmic0)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../../data/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
source("../../3_analysis/helper/pcawg.colour.palette.R") ## pcawg_palette
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
read_info <- function(ct){
  .x <- list(diagRE_DMDL_SP = try(readRDS(paste0("../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS"))),
             fullRE_DMDL_SP = try(readRDS(paste0("../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS"))),
             dataset_active_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../data/", override_warning_X_Z = T),
             dataset_nucleotidesubstitution3 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution3", path_to_data = "../../../data/", override_warning_X_Z = T))
  .x
}
read_info_list <- lapply(enough_samples, function(ct){
  read_info(ct)
}); names(read_info_list) <- enough_samples
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## From the PCAWG dataset, as shown in the paper
diagRE_DM <- sapply(read_info_list, function(i) i$diagRE_DMDL_SP, simplify = F)
diagRE_DM_tests <-  sapply(diagRE_DM, CompSign::wald_TMB_wrapper, simplify = F)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
pcawg_palette <- pcawg.colour.palette(x = gsub("\\..*", "", names(read_info_list)),  scheme = "tumour.subtype")
pcawg_palette_2 <- c(  "Bone-Osteosarc"="#FFD700", "Breast-AdenoCA"="#CD6090" , "CNS-GBM"="#3D3D3D",
                       "CNS-Medullo"="#D8BFD8", "CNS-PiloAstro"="#B0B0B0", "ColoRect-AdenoCA"="#191970",
                       "Eso-AdenoCA"="#1E90FF",  "Head-SCC" = "#8B2323", "Kidney-ChRCC"= "#B32F0B", 
                       "Kidney-RCC.clearcell"= "#FF4500", "Kidney-RCC.papillary"= "#FF4500",    "Liver-HCC"= "#006400",
                       "Lung-SCC" =  "#FDF5E6", "Lymph-BNHL"= "#698B22", "Lymph-CLL"= "#F4A35D", "Ovary-AdenoCA"= "#008B8B",
                       "Panc-AdenoCA" = "#7A378B", "Panc-Endocrine"= "#E066FF", "Prost-AdenoCA"= "#87CEFA",
                       "Skin-Melanoma.cutaneous"= "#000000",  "Stomach-AdenoCA"= "#BFEFFF" , "Thy-AdenoCA"= "#9370DB",
                       "Uterus-AdenoCA"=   "#FF8C69")
cbind(pcawg_palette, pcawg_palette_2)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Get signatures when active signatures are modified, and run the model to compare results
## leave_one_out_exposures: list of exposures, for each cancer type, where one of each active signature has been left out
## diagDM_leave_one_out_exposures: The model run on these exposures

give_leave_one_out_exposures <- function(signature_fitting_method){
  sapply(enough_samples, function(ct){
    active_sigs_in_ct <- colnames(read_info_list[[ct]]$dataset_active_sigs$Y)
    x <- lapply(active_sigs_in_ct, function(active_sig_it){
      extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                           subset_signatures = active_sigs_in_ct[!(active_sigs_in_ct %in% active_sig_it)],
                           signature_version=NULL, signature_definition = sigs_cosmic0,
                           signature_fitting_method = signature_fitting_method)
    })
    names(x) <- paste0('-', active_sigs_in_ct)
    x
  }, simplify = F)
}

give_add_one_signature_exposures <- function(signature_fitting_method){
  sapply(enough_samples, function(ct){
    active_sigs_in_ct <- colnames(read_info_list[[ct]]$dataset_active_sigs$Y)
    ## sample 4 signatures from sigs_cosmic0
    additional_signatures <- sample(colnames(sigs_cosmic0)[! (colnames(sigs_cosmic0) %in% active_sigs_in_ct)], 4)
    x <- lapply(additional_signatures, function(additional_signature_ct){
      extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                           subset_signatures = c(additional_signature_ct, active_sigs_in_ct),
                           signature_version=NULL, signature_definition = sigs_cosmic0,
                           signature_fitting_method = signature_fitting_method)
    })
    names(x) <- paste0('+', additional_signatures)
    x
  }, simplify = F)
}


##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Run model

wrapper_leave_one_out_and_add_one_exposures <- function(dataset_name, signature_fitting_method){
  if(signature_fitting_method == 'mutSigExtractor'){
    add_to_name = ''
  }else{
    add_to_name = signature_fitting_method
  }
  
  if(dataset_name=='leave_one_out'){
    if(signature_fitting_method == 'mutSigExtractor'){
      df = leave_one_out_exposures
    }else if(signature_fitting_method == 'QP'){
      df = leave_one_out_exposures_QP
    }
  }else if(dataset_name == 'add_one_signature'){
    if(signature_fitting_method == 'mutSigExtractor'){
      df = add_one_signature_exposures
    }else if(signature_fitting_method == 'QP'){
      df = add_one_signature_exposures_QP
    }
  }else{
    stop('Unknown <dataset_name>')
  }

  out <- sapply(enough_samples, function(ct){
    sapply(df[[ct]], function(it_from_dataset){
      try(CompSign::wrapper_run_TMB(object = it_from_dataset,
                                model = "diagRE_DM", use_nlminb=T, smart_init_vals=F))
    }, simplify = F)
  }, simplify = F)
  
  fileout <- paste0("../../../data/assessing_models_real_data/inference_results/TMB/1a/1a_", dataset_name, "_diagREDM",
         add_to_name, ".RDS")
  
  saveRDS(out,
          file = fileout)
  cat('Output file ', fileout, ' saved.\n')
}


# ct='Panc-AdenoCA'
# leave_one_out_it='SBS1'
# diagDM_leave_one_out_exposures[[ct]][[paste0('-', leave_one_out_it)]]
# extra_run = CompSign::wrapper_run_TMB(object = leave_one_out_exposures[[ct]][[paste0('-', leave_one_out_it)]],
#                                       model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
# extra_run
# diagDM_leave_one_out_exposures[[ct]][[paste0('-', leave_one_out_it)]] <- extra_run
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##

# ct='Skin-Melanoma.cutaneous'
# add_one_it='SBS45'
# diagDM_add_one_signature_exposures[[ct]][[paste0('+', add_one_it)]]
# extra_run = CompSign::wrapper_run_TMB(object = add_one_signature_exposures[[ct]][[paste0('+', add_one_it)]],
#                                       model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
# extra_run
# diagDM_add_one_signature_exposures[[ct]][[paste0('+', add_one_it)]] <- extra_run

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
re_run <- FALSE
if(re_run){
  ## Simulate datasets and run model
  leave_one_out_exposures <- give_leave_one_out_exposures(signature_fitting_method = 'mutSigExtractor')
  
  add_one_signature_exposures <- give_add_one_signature_exposures(signature_fitting_method = 'mutSigExtractor')
  wrapper_leave_one_out_and_add_one_exposures(dataset_name = 'leave_one_out', signature_fitting_method = 'mutSigExtractor')
  
  
  ## QP signature extraction
  ## Creating datasets
  leave_one_out_exposures_QP <- give_leave_one_out_exposures(signature_fitting_method = 'QP')
  add_one_signature_exposures_QP <- give_add_one_signature_exposures(signature_fitting_method = 'QP')
  
  ## save datasets
  saveRDS(leave_one_out_exposures, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/1a/1a_leave_one_out_exposures.RDS"))
  saveRDS(add_one_signature_exposures, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/1a/1a_add_one_signature_exposures.RDS"))
  saveRDS(leave_one_out_exposures_QP, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/1a/1a_leave_one_out_exposures_QP.RDS"))
  saveRDS(add_one_signature_exposures_QP, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/1a/1a_add_one_signature_exposures_QP.RDS"))
  
  ## running code
  wrapper_leave_one_out_and_add_one_exposures(dataset_name = 'leave_one_out', signature_fitting_method = 'QP')
  wrapper_leave_one_out_and_add_one_exposures(dataset_name = 'add_one_signature', signature_fitting_method = 'QP')
  
}
##-----------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------##
## re-run if necessary (leave-one-out)
# name_diagRE_temp <- "../../../data/assessing_models_real_data/inference_results/TMB/leave_one_out_diagREDMQP.RDS"
# diagDM_leave_one_out_exposures_QP_to_rerun <- readRDS(file = name_diagRE_temp)
# max(sapply(diagDM_leave_one_out_exposures_QP_to_rerun, function(i) length(table(sapply(i, 'class'))))) ## if 2, there are errors
# # ct <- 'Prost-AdenoCA'

# which(sapply(diagDM_leave_one_out_exposures_QP_to_rerun[[ct]], 'class') == 'try-error')
# j = 5
# diagDM_leave_one_out_exposures_QP_to_rerun[[ct]][[j]]
# tmp_res <- CompSign::wrapper_run_TMB(object = leave_one_out_exposures_QP[[ct]][[j]],
#                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
# tmp_res
# diagDM_leave_one_out_exposures_QP_to_rerun[[ct]][[j]] <- tmp_res
# saveRDS(diagDM_leave_one_out_exposures_QP_to_rerun, file = name_diagRE_temp)

##-----------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## reading the results from code which were saved in the function wrapper functions
diagDM_leave_one_out_exposures_QP <- readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1a/1a_leave_one_out_diagREDMQP.RDS"))
diagDM_add_one_signature_exposures_QP <- readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1a/1a_add_one_signature_diagREDMQP.RDS"))

diagDM_leave_one_out_exposures_mutSigExtractor <- readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1a/1a_leave_one_out_diagREDM.RDS"))
diagDM_add_one_signature_exposures_mutSigExtractor <- readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1a/1a_add_one_signature_diagREDM.RDS"))

##-----------------------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------------------##
## Tests for DA
diagDM_leave_one_out_exposures_tests_mutSigExtractor <- sapply(diagDM_leave_one_out_exposures, function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE)

diagDM_add_one_signature_exposures_tests_mutSigExtractor <- sapply(diagDM_add_one_signature_exposures, function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE)

diagDM_leave_one_out_exposures_tests_QP <- sapply(diagDM_leave_one_out_exposures_QP, function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE)

diagDM_add_one_signature_exposures_tests_QP <- sapply(diagDM_add_one_signature_exposures_QP, function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE)

##-----------------------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------##
## Plots: general correlation between QP and mutSigExtractor
source("../../1_create_ROO/helper_1_create_ROO.R") ## for QPsig
## QPsig is the function used for the PCAWG data in the paper. extract_sigs is a wrapper
## example of comparison of exposures
xxx1 <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                             subset_signatures = active_sigs_in_ct[!(active_sigs_in_ct %in% active_sig_it)],
                             signature_version=NULL, signature_definition = sigs_cosmic0, signature_fitting_method = 'mutSigExtractor')

xxxQP <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                              subset_signatures = active_sigs_in_ct[!(active_sigs_in_ct %in% active_sig_it)],
                              signature_version=NULL, signature_definition = sigs_cosmic0, signature_fitting_method = 'QP')

plot(log(as.vector(xxx1$Y)), log(as.vector(xxxQP$Y)))
##-----------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------------------##
## Plots: correlation of betas (function)
theme_set(theme_bw())
table(is.na(beta_cor_leave_one_out_exposures$value))

plot_beta_cor <- function(title_arg, remove_col_legend=F, additional_title='', ...){
  beta_cor_leave_one_out_exposures <- give_beta_cor(...)
  beta_cor_leave_one_out_exposures <- melt(beta_cor_leave_one_out_exposures)
  ## remove NAs
  beta_cor_leave_one_out_exposures <- beta_cor_leave_one_out_exposures[!is.na(beta_cor_leave_one_out_exposures$value),]
  beta_cor_leave_one_out_exposures <- beta_cor_leave_one_out_exposures %>% group_by(L1) %>% arrange(-value, .by_group = TRUE)
  beta_cor_leave_one_out_exposures$x = paste0(beta_cor_leave_one_out_exposures$L2, beta_cor_leave_one_out_exposures$L1)
  beta_cor_leave_one_out_exposures$x <- factor(beta_cor_leave_one_out_exposures$x, levels=beta_cor_leave_one_out_exposures$x)
  res <- ggplot(beta_cor_leave_one_out_exposures, aes(x=x, y=value, col=(L1)))+geom_point()+geom_line(aes(group=L1))+
    # facet_wrap(.~L1, scales = 'free_x', drop = T, nrow=1)+
    scale_color_manual(values = c(pcawg_palette_2, "white"="white"), na.value="cyan")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    labs(y=paste0('Correlation', additional_title), x=title_arg, col='Cancer type')+
    geom_hline(yintercept = 0.95, lty='dashed')+
    geom_hline(yintercept = 0.90, lty='dashed')+
    geom_hline(yintercept = 0.5, lty='dashed')
  if(remove_col_legend){
    res <- res+guides(col='none')
  }
  return(res)
}
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------##
## Plots: correlation of betas (results)
give_all_plots_leave_one_out_add_one <- function(res_list_with_modifications_LOO, res_list_with_modifications_ADD, res_list){
  cowplot::plot_grid(plot_beta_cor(res_list_with_modifications_LOO, res_list, mode='leave_one_out_all_betas',
                                   title_arg=latex2exp::TeX(r"(One signature removed ($\beta_0, \beta_1$))"), remove_col_legend=T, additional_title = ' in\nleave-one-out (D1A)'),
                     plot_beta_cor(res_list_with_modifications_LOO, res_list, mode='leave_one_out_beta_slope',
                                   title_arg=latex2exp::TeX(r"(One signature removed ($\beta_1$))"), remove_col_legend=T, additional_title = ' in\nleave-one-out (D1A)'),
                     plot_beta_cor(res_list_with_modifications_ADD, res_list, mode='add_one',
                                   title_arg=latex2exp::TeX(r"(One added signature ($\beta_0, \beta_1$))"), remove_col_legend=T, additional_title = 'in\nadd-one (D1B)'),
                     plot_beta_cor(res_list_with_modifications_ADD, res_list, mode='add_one_beta_slope',
                                   title_arg=latex2exp::TeX(r"(One added signature ($\beta_1$))"), remove_col_legend=T, additional_title = ' in\nadd-one\n(D1B)'), nrow=1)
}

## for mutSigExtractor 
pdf("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/1_robustness_of_results_varying_active_signatures_LOO_ADD_mutSigExtractor.pdf", height = 2, width = 12)
give_all_plots_leave_one_out_add_one(diagDM_leave_one_out_exposures, diagDM_add_one_signature_exposures, diagRE_DM)
dev.off()
png("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/1_robustness_of_results_varying_active_signatures_LOO_ADD_mutSigExtractor.png", height = 2, width = 12, units = 'in', res=300)
give_all_plots_leave_one_out_add_one(diagDM_leave_one_out_exposures, diagDM_add_one_signature_exposures, diagRE_DM)
dev.off()

## to plot legend
# plot_beta_cor(diagDM_leave_one_out_exposures, diagRE_DM, mode='leave_one_out_all_betas',
#               title_arg=latex2exp::TeX(r"(One signature removed ($\beta_0, \beta_1$))"), remove_col_legend=F)+
#   theme(legend.position = 'bottom', legend.spacing.y = unit(-0.3, 'cm'), legend.title = element_blank())+
#   guides(color=guide_legend(ncol=8))

## for QP
pdf("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/1_robustness_of_results_varying_active_signatures_LOO_ADD_QP.pdf", height = 2, width = 12)
give_all_plots_leave_one_out_add_one(diagDM_leave_one_out_exposures_QP, diagDM_add_one_signature_exposures_QP, diagRE_DM)
dev.off()
png("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/1_robustness_of_results_varying_active_signatures_LOO_ADD_QP.png", height = 2, width = 12, units = 'in', res=300)
give_all_plots_leave_one_out_add_one(diagDM_leave_one_out_exposures_QP, diagDM_add_one_signature_exposures_QP, diagRE_DM)
dev.off()


table(unlist(sapply(diagDM_leave_one_out_exposures_QP, function(i) sapply(i, 'class'))))
table(unlist(sapply(diagDM_add_one_signature_exposures_QP, function(i) sapply(i, 'class'))))

##-----------------------------------------------------------------------------------------------------##
## Plots: agreement of DA results between the original (non-modified) dataset and the modified datasets

plt1_mutSigExtractor <- give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_leave_one_out_exposures_tests_mutSigExtractor, ylabel='Number of datasets (leave-one-out)', colour_version = 'colourversion2')+guides(fill='none')
plt2_mutSigExtractor <- give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_add_one_signature_exposures_tests_mutSigExtractor, ylabel='Number of datasets (add-one)', colour_version = 'colourversion2')+guides(fill='none')

plt1_QP <- give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_leave_one_out_exposures_tests_QP, ylabel='Number of datasets (leave-one-out)', colour_version = 'colourversion2')+guides(fill='none')
plt2_QP <- give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_add_one_signature_exposures_tests_QP, ylabel='Number of datasets (add-one)', colour_version = 'colourversion2')+guides(fill='none')

# plt1_legend <- ggpubr::get_legend(give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_leave_one_out_exposures_tests, ylabel='Number of datasets (add-one)', colour_version = 'colourversion2'))
# plt2_legend <- ggpubr::get_legend(give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_add_one_signature_exposures_tests, ylabel='Number of datasets (add-one)', colour_version = 'colourversion2'))
# plot(plt1_legend)
# plot(plt2_legend)

## mutSigExtractor
pdf("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/check1_DAagrement_mutSigExtractor.pdf",
    height = 3, width = 6.2)
cowplot::plot_grid(cowplot::plot_grid(plt1_mutSigExtractor, plt2_mutSigExtractor, nrow=1),
                   ggdraw()+draw_image("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/DA_agreement_legend.png"),
                   nrow=2, rel_heights = c(6,1))
dev.off()
png("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/check1_DAagrement_mutSigExtractor.png",
    height = 3, width = 6.2, units = 'in', res=300)
cowplot::plot_grid(cowplot::plot_grid(plt1_mutSigExtractor, plt2_mutSigExtractor, nrow=1),
                   ggdraw()+draw_image("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/DA_agreement_legend.png"),
                   nrow=2, rel_heights = c(6,1))
dev.off()

## QP
pdf("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/check1_DAagrement_QP.pdf",
    height = 3, width = 6.2)
cowplot::plot_grid(cowplot::plot_grid(plt1_QP, plt2_QP, nrow=1),
                   ggdraw()+draw_image("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/DA_agreement_legend.png"),
                   nrow=2, rel_heights = c(6,1))
dev.off()
png("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/check1_DAagrement_QP.png",
    height = 3, width = 6.2, units = 'in', res=300)
cowplot::plot_grid(cowplot::plot_grid(plt1_QP, plt2_QP, nrow=1),
                   ggdraw()+draw_image("../../../results/results_TMB/3_real_data_assessment/1_robustness_of_results_varying_active_signatures/DA_agreement_legend.png"),
                   nrow=2, rel_heights = c(6,1))
dev.off()
