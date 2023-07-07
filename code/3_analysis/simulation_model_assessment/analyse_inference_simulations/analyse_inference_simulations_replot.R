## replot the output of analyse_inference_simulations_integrate.R, if more plots are needed
## this way the p-values for competing methods do not need to be re-computed

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(ggrepel)
library(dplyr)
library(latex2exp)
library(reshape2)

source("helper_model_assessment.R")
source("../../../2_inference_TMB/helper_TMB.R")
source("../../../1_create_ROO/roo_functions.R")

multiple_runs <- T
if(multiple_runs){
  flder_out <- "../../../../results/results_TMB/simulated_datasets/mixed_effects_models_multiple/"
  flder_in <- "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries_multiple/"
}else{
  flder_out <- "../../../../results/results_TMB/simulated_datasets/mixed_effects_models/"
  flder_in <- "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries/"
}

# datasetgeneration=['GenerationMixturePCAWG', #'GenerationMixturefewersignaturesPCAWG', \
#                    'GenerationMixturefewersignaturespairedPCAWG', \
#                    'GenerationMixturefewersignaturespairedstomachPCAWG',\
#                    'GenerationMixturefewersignaturespairedKidneyRCCPCAWG',\

# generation <- c('GenerationMixturefewersignaturespairedKidneyRCCPCAWG')
# generation <- c('GenerationMixturefewersignaturespairedPCAWG')
# generation <- c('GenerationMixturefewersignaturespairedstomachPCAWG')
# generation <- c('GenerationMixturefewersignaturesPCAWG')
# generation <- c('GenerationMixturefewersmallsignaturespairedKidneyRCCPCAWG')

for(generation in c(
  # 'GenerationHnormtwolambdas', 'GenerationJnorm', 'GenerationK',
  # 'GenerationMixturefewersignaturespairedKidneyRCCPCAWG', 'GenerationMixturefewersignaturespairedPCAWG',
  #                   'GenerationMixturefewersignaturespairedstomachPCAWG', 'GenerationMixturefewersignaturesPCAWG',
  #                   'GenerationMixturefewersmallsignaturespairedKidneyRCCPCAWG',
  #                   'GenerationMixturefewersignaturespairedProstAdenoCAPCAWG',
  #                   'GenerationMixturefewersignaturespairedCNSGBMPCAWG',
  #                   'GenerationMixturefewersignaturespairedObsNmPancEndocrinePCAWG',
  #                   'GenerationMixturefewersignaturespairedObsNmUterusAdenoCAPCAWG',
  # 'GenerationJnormBTwoLambdasOneChangingBeta',
  # 'GenerationPois',
  'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGProstAdenoCAPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMOvaryAdenoCAPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMLungSCCPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMKidneyRCCpapillaryPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMPancEndocrinePCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMColoRectAdenoCAPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMLymphBNHLPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMHeadSCCPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMEsoAdenoCAPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMLymphCLLPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMPancAdenoCAPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMKidneyChRCCPCAWG',
  'GenerationMixturefewersignaturespairedObsNmObsDMCNSGBMPCAWG',
  'GenerationMixtureallsignaturespairedObsNmCNSGBMPCAWG',
  'GenerationMixturefewersignaturespairedObsNmPoissonSigPCAWGColoRectAdenoCAPCAWG',
  'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGEsoAdenoCAPCAWG',
  'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGHeadSCCPCAWG',
  'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGLungSCCPCAWG',
  'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGKidneyChRCCPCAWG'
                    )){
  
  try(rm(a))
  a <- readRDS(paste0("../../../../data/assessing_models_simulation/summaries_synthetic_DA/", generation, ".RDS"))
  
  all_converged <- !apply(a$pvals_data_frame, 1, function(i) any(is.na(i)))
  a$pvals_data_frame_adj <- a$pvals_data_frame
  
  a$pvals_data_frame_adj <- adjust_all(a$pvals_data_frame_adj, new_version = T, method = "bonferroni")
  varying_betashape <-give_accuracies_with_varying_var(var = 'beta_gamma_shape',
                                                       datasets_arg = a$datasets,
                                                       pvals_data_frame_arg = a$pvals_data_frame)
  varying_betashape_adj <-give_accuracies_with_varying_var(var = 'beta_gamma_shape',
                                                       datasets_arg = a$datasets,
                                                       pvals_data_frame_arg = a$pvals_data_frame_adj)
  varying_betashapeAC <-give_accuracies_with_varying_var(var = 'beta_gamma_shape',
                                                       datasets_arg = a$datasets[all_converged],
                                                       pvals_data_frame_arg = a$pvals_data_frame[all_converged,])
  
  varying_betashape_n_adj <-give_accuracies_with_varying_var(var = c('beta_gamma_shape', 'n'), two_var = T,
                                                       datasets_arg = a$datasets,
                                                       pvals_data_frame_arg = a$pvals_data_frame_adj)
  
  varying_betashape_n <-give_accuracies_with_varying_var(var = c('beta_gamma_shape', 'n'), two_var = T,
                                                             datasets_arg = a$datasets,
                                                             pvals_data_frame_arg = a$pvals_data_frame)
  varying_betashape_d <-give_accuracies_with_varying_var(var = c('beta_gamma_shape', 'd'), two_var = T,
                                                         datasets_arg = a$datasets,
                                                         pvals_data_frame_arg = a$pvals_data_frame)
  
  if((generation %in% c("GenerationMixturePCAWG", "GenerationMixturefewersignaturesPCAWG", "GenerationMixturefewersignaturespairedPCAWG")) | grepl('GenerationMixturefewersignaturespaired', generation) ){
    varying_betashape$beta_gamma_shape <- signif(varying_betashape$beta_gamma_shape, 2)
    varying_betashape_adj$beta_gamma_shape <- signif(varying_betashape_adj$beta_gamma_shape, 2)
  }
  
  varying_betashape$model <- gsub("pvals_", "", varying_betashape$model)
  varying_betashape_adj$model <- gsub("pvals_", "", varying_betashape_adj$model)
  
  .xx <- varying_betashape[which(varying_betashape$beta_gamma_shape == max(varying_betashape$beta_gamma_shape)),]
  ## sometimes there are problems with numeric precision. if that is the case, select as the last value only one of the selected rows
  .xx <- .xx[!duplicated(.xx$model),]
  
  remove_duplicated_rows <- function(.xx){
    .xx <- .xx[!duplicated(.xx$model),]
  }
  
  
  if(generation == 'GenerationMixturefewersignaturespairedKidneyRCCPCAWG'){
    title = 'Kidney-RCC (paired, larger sigs)'
  }else if( generation == 'GenerationMixturefewersignaturespairedPCAWG'){
    title = 'Liver-HCC (paired, larger sigs)'
  }else if ( generation  == 'GenerationMixturefewersignaturespairedstomachPCAWG'){
    title = 'Stomach-AdenoCa (paired, larger sigs)'
  }else if( generation == 'GenerationMixturefewersignaturesPCAWG'){
    title = 'Liver-HCC (not paired, larger sigs)'
  }else if( generation == 'GenerationMixturefewersmallsignaturespairedKidneyRCCPCAWG'){
    title = 'Kidney-RCC (paired, small sigs)'
  }else if( generation == 'GenerationMixturefewersignaturespairedProstAdenoCAPCAWG'){
    title = 'Prost-AdenoCA (paired, larger sigs)'
  }else if( generation == 'GenerationMixturefewersignaturespairedCNSGBMPCAWG'){
    title = 'CNS-GBM (paired, larger sigs)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmPancEndocrinePCAWG'){
    title = 'Panc-Endocrine (paired, larger sigs)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmUterusAdenoCAPCAWG'){
    title = 'Uterus-AdenoCA (paired, larger sigs)'
  }else if( generation ==  'GenerationMixtureallsignaturespairedObsNmCNSGBMPCAWG'){
    title = 'CNS-GBM (paired, all sigs)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMOvaryAdenoCAPCAWG'){
    title = 'Ovary-Adeno (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMLungSCCPCAWG'){
    title = 'Lung-SCC (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMKidneyRCCpapillaryPCAWG'){
    title = 'Kidney-RCCp (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMPancEndocrinePCAWG'){
    title = 'Panc-Endocrine (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMColoRectAdenoCAPCAWG'){
    title = 'ColoRect-AdenoCA (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMLymphBNHLPCAWG'){
    title = 'Lymph-BNHL (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMHeadSCCPCAWG'){
    title = 'Head-SCC (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMEsoAdenoCAPCAWG'){
    title = 'Eso-Adeno (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMLymphCLLPCAWG'){
    title = 'Lymph-CLL (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMPancAdenoCAPCAWG'){
    title = 'Panc-AdenoCA (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMKidneyChRCCPCAWG'){
    title = 'Kidney-ChRCC (paired, larger sigs, DM)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMCNSGBMPCAWG'){
    title = 'CNS-GBM (paired, larger sigs, DM)'
  }else if(generation == 'GenerationJnormBTwoLambdasOneChangingBeta'){
    title = "Single category change"
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmPoissonSigPCAWGColoRectAdenoCAPCAWG'){
    title = 'Colo-RectAdenoCA (paired, larger sigs, Poisson)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGEsoAdenoCAPCAWG'){
    title = 'Eso-AdenoCA (paired, larger sigs, Norm)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGHeadSCCPCAWG'){
    title = 'Head-SCC (paired, larger sigs, Gaussian)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGLungSCCPCAWG'){
    title = 'Lung-SCC (paired, larger sigs, Gaussian)'
  }else if( generation == 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGKidneyChRCCPCAWG'){
    title = 'Kidney-ChRCC (paired, larger sigs, Norm)'
  }else{
    title = generation
  }

  varying_betashape$beta_gamma_shape <- signif(varying_betashape$beta_gamma_shape, 2)
  varying_betashape_adj$beta_gamma_shape <- signif(varying_betashape_adj$beta_gamma_shape, 2)
  
  title_x <- 'Percentage of mixture'
  if(generation %in% c("GenerationJnormBTwoLambdasOneChangingBeta", "GenerationPois", "GenerationJnorm", "GenerationHnormtwolambdas", "GenerationK")){
    title_x <- TeX("$\\gamma$")
  }
  
  ggplot(varying_betashape, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
                                lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
    ggtitle(title)+
    geom_label_repel(data = .xx,
                     max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
                     nudge_x=max(varying_betashape$beta_gamma_shape)+2,
                     size=2.5)+
    coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.3))+ ## used to be *1.5
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_palette2_factorv2.pdf"),
         height = 3.0, width = 4.0)
  
  ggplot(varying_betashape_adj, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
                                lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
    ggtitle(title)+
    geom_label_repel(data = remove_duplicated_rows(varying_betashape_adj[which(varying_betashape_adj$beta_gamma_shape == max(varying_betashape_adj$beta_gamma_shape)),]),
                     max.overlaps = Inf, aes(x=factor(max(varying_betashape_adj$beta_gamma_shape)),
                                             lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')), direction = "y",
                     nudge_x=max(varying_betashape_adj$beta_gamma_shape)+2,
                     size=2.5)+
    coord_cartesian(xlim = c(0, length(unique(varying_betashape_adj$beta_gamma_shape))*1.5))+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_adjbonferroni_with_betagammashape_palette2_factorv2.pdf"),
         height = 3.0, width = 4.0)
  
  system(paste0("open ", flder_out, generation, "/summaries/"))
  
  ggplot(varying_betashapeAC, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
                                lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
    ggtitle(title)+
    geom_label_repel(data = .xx,
                     max.overlaps = Inf, aes(x=factor(max(varying_betashapeAC$beta_gamma_shape))), direction = "y",
                     nudge_x=max(varying_betashapeAC$beta_gamma_shape)+2,
                     size=2.5)+
    coord_cartesian(xlim = c(0, length(unique(varying_betashapeAC$beta_gamma_shape))*1.5))+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
  
  
  varying_betashape_n_adj$beta_gamma_shape <- signif(varying_betashape_n_adj$beta_gamma_shape, 2)
  varying_betashape_n$beta_gamma_shape <- signif(varying_betashape_n$beta_gamma_shape, 2)
  varying_betashape_d$beta_gamma_shape <- signif(varying_betashape_d$beta_gamma_shape, 2)
  ggplot(varying_betashape_n_adj, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
                                lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
    ggtitle(title)+
    geom_label_repel(data = .xx,
                     max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
                     nudge_x=max(varying_betashape$beta_gamma_shape)+2,
                     size=2.5)+
    coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.5))+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
    facet_wrap(.~n)
  
  ggplot(varying_betashape_n_adj[varying_betashape_n_adj$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
         aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
                                      lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
    ggtitle(title)+
    geom_label_repel(data = .xx[varying_betashape_n_adj$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
                     max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
                     nudge_x=max(varying_betashape$beta_gamma_shape)+2,
                     size=2.5)+
    coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.5))+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
    facet_wrap(.~n)
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_adjbonferroni_with_betagammashape_palette2_factorv2_subset.pdf"),
         height = 3.0, width = 6.6)
  
  ggplot(varying_betashape_n[varying_betashape_n$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
         aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
             lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
    ggtitle(title)+
    geom_label_repel(data = .xx[varying_betashape_n$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
                     max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
                     nudge_x=max(varying_betashape$beta_gamma_shape)+2,
                     size=2.5)+
    coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.5))+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
    facet_wrap(.~n)
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_palette2_factorv2_subset.pdf"),
         height = 3.0, width = 6.6)
  
  ggplot(varying_betashape_d[varying_betashape_d$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
         aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
             lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
    ggtitle(title)+
    geom_label_repel(data = .xx[varying_betashape_d$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
                     max.overlaps = Inf, aes(x=factor(max(varying_betashape_d$beta_gamma_shape))), direction = "y",
                     nudge_x=max(varying_betashape_d$beta_gamma_shape)+2,
                     size=2.5)+
    coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.5))+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
    facet_wrap(.~d)
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_facetD_palette2_factorv2_subset.pdf"),
         height = 3.0, width = 6.6)
  
  ggplot(varying_betashape_n_adj[varying_betashape_n_adj$beta_gamma_shape == 0,],
         aes(x=n, y = Accuracy, col=model, group=model, label=model,
                                  lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models2)+labs(col='', x='Number of samples')+
    guides(lty='none')+
    ggtitle(title)+
    geom_label_repel(data = remove_duplicated_rows(varying_betashape_n_adj[which(varying_betashape_n_adj$beta_gamma_shape == 0),]),
                     max.overlaps = Inf, aes(x=(max(varying_betashape_n_adj$n))), direction = "y",
                     size=2.5)+
    coord_cartesian(xlim = c(0, max(varying_betashape_n_adj$n)*1.5))+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
  
  pvals_noDA <- a$pvals_data_frame[as.vector(sapply(a$datasets, function(i) i$beta_gamma_shape)) == 0,]
  pvals_DA <- a$pvals_data_frame[as.vector(sapply(a$datasets, function(i) i$beta_gamma_shape)) != 0,]
  plot(pvals_noDA$pvals_diagREDM, pvals_noDA$pvals_fullREM); abline(coef = c(0,1))
  
  table(pvals_noDA$pvals_diagREDM < 0.05)
  table(pvals_noDA$pvals_fullREM < 0.05)
  
  hist(pvals_noDA$pvals_diagREDM, breaks = 20)
  hist(pvals_noDA$pvals_fullREM, breaks = 20)
  
  pvals_noDA_adj <- a$pvals_data_frame_adj[as.vector(sapply(a$datasets, function(i) i$beta_gamma_shape)) == 0,]
  plot(pvals_noDA_adj$pvals_diagREDM, pvals_noDA_adj$pvals_fullREM); abline(coef = c(0,1), main=generation)
  table(pvals_noDA_adj$pvals_diagREDM < 0.05)
  table(pvals_noDA_adj$pvals_fullREDMSL < 0.05)
  table(pvals_noDA_adj$pvals_fullREM < 0.05)
  

  length_out_rocauc <- 50
  rocauc <- lapply(seq(0, 1, length.out = length_out_rocauc), function(thresholdit) (as(t(sapply(colnames(a$pvals_data_frame), function(col_it){
    summarise_DA_detection(true = a$pvals_data_frame$true, predicted = a$pvals_data_frame[,col_it] < thresholdit, verbose = F)[c('Specificity', 'Sensitivity')]
  })), 'matrix')))
  names(rocauc) <- seq(0, 1, length.out = length_out_rocauc)
  # rocauc <- lapply(rocauc, as.matrix)
  rocauc <- (reshape2::melt(rocauc, id.vars=c('Specificity', 'Sensitivity')))
  rocauc <- dcast(rocauc, Var1+L1~Var2, value.var = "value")
  rocauc$model <- gsub("pvals_", "", rocauc$Var1)
  rocauc <- rocauc[rocauc$model != "true",]
  unique(rocauc$model)
  rocauc$model[rocauc$model == "ttest_ilr_adj"] = "ILR"
  rocauc$model[rocauc$model == "ttest_props"] = "ttest"
  
  ggplot(rocauc, aes(x=1-Specificity, y= Sensitivity, col=model,
                     lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_line()+theme_bw()+theme(legend.position = "bottom")+
    scale_color_manual(values=colours_models2)+labs(col='', x='1- Specificity', y='Sensitivity')+guides(col='none', lty='none')
  ggsave(paste0(flder_out, generation, "/summaries/ROCAUCcurve_palette.pdf"),
         height = 3.0, width = 3)
  
  ggplot(rocauc, aes(x=1-Specificity, y= Sensitivity, col=model,
                     lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_line()+theme_bw()+theme(legend.position = "bottom")+
    scale_color_manual(values=colours_models2)+labs(col='', x='1- Specificity', y='Sensitivity')+guides(col='none', lty='none')+facet_wrap(.~model)
  
  ## for fixed cut-offs, changing effect sizes
  ## information about effect sizes
  beta_gamma_shapes <- as.vector(unlist(sapply(a$datasets, `[`, 'beta_gamma_shape')))
  ns <- as.vector(unlist(sapply(a$datasets, `[`, 'n')))
  thresholdit = 0.05
  rocauc_2 <- lapply(sort(unique(ns)), function(betagammashapeit){
    (as(t(sapply(colnames(a$pvals_data_frame), function(col_it){
    summarise_DA_detection(true = a$pvals_data_frame$true[ns == betagammashapeit],
                           predicted = a$pvals_data_frame[ns == betagammashapeit,col_it] < thresholdit, verbose = F)[c('Specificity', 'Sensitivity')]
    })), 'matrix'))
  })
  names(rocauc_2) <- sort(unique(ns))
  rocauc_2 <- (reshape2::melt(rocauc_2, id.vars=c('Specificity', 'Sensitivity')))
  rocauc_2 <- dcast(rocauc_2, Var1+L1~Var2, value.var = "value")
  rocauc_2$model <- gsub("pvals_", "", rocauc_2$Var1)
  rocauc_2$model[rocauc_2$model == "ttest_ilr_adj"] = "ILR"
  rocauc_2$model[rocauc_2$model == "ttest_props"] = "ttest"
  ggplot(rocauc_2, aes(x=1-Specificity, y= Sensitivity, col=model,
                     lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_line()+theme_bw()+theme(legend.position = "bottom")+
    scale_color_manual(values=colours_models2)+labs(col='', x='1- Specificity', y='Sensitivity')+guides(col='none', lty='none')

  
}

a$pvals_data_frame$pvals_diagREDM[log(as.vector(sapply(a$datasets, function(i) i$beta_gamma_shape))) == -1.2]
