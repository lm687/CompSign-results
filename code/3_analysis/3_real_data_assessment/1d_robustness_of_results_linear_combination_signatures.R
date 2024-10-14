## Select a signature at random, then divide it into two such that the first one is a linear combination of the other two.
## Vary the weight of the two in the combination

##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
source("../../1_create_ROO/helper_1_create_ROO.R") ## for QPsig
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

theme_set(theme_bw())

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
weight_first_sig_in_linear_comb = c(0.75, 0.5, 0.25) ## signature under analysis is a combination of 75% of the first and 25% of the second artificial signature,
##.                                                                                                  50-50%
##.                                                                                                  25-70%
nreplicates = 5
##-----------------------------------------------------------------------------------------------------##

active_sigs_in_ct_list <- sapply(enough_samples, function(ct){
  active_sigs_in_ct <- colnames(read_info_list[[ct]]$dataset_active_sigs$Y)
  active_sigs_in_ct
})

additional_signatures_list <- sapply(enough_samples, function(ct){
  active_sigs_in_ct <- active_sigs_in_ct_list[[ct]]
  sample(active_sigs_in_ct, 1)
})


## signature exposures with the signatures in the order above
give_original_signature_exposures <- function(signature_fitting_method){
  sapply(enough_samples, function(ct){
    ## sample 1 signature from active_sigs_in_ct. These will be replaced by two signatures such that the original signature is a linear combination of the other two
    additional_signature <- additional_signatures_list[[ct]]
    x <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                           subset_signatures = c(additional_signature, active_sigs_in_ct_list[[ct]][!(active_sigs_in_ct_list[[ct]] == additional_signature)]),
                           signature_version=NULL, signature_definition = sigs_cosmic0,
                           signature_fitting_method = signature_fitting_method)
    x
  }, simplify = F)
}

## signature under consideration is a linear combination of the two but there is no correlation between these two, or little correlation
give_linearlycombined_nocor_signature_exposures <- function(signature_fitting_method){
  sapply(enough_samples, function(ct){
    ## sample 1 signature from active_sigs_in_ct. These will be replaced by two signatures such that the original signature is a linear combination of the other two
    additional_signature <- additional_signatures_list[[ct]]
    ## get two signatures which are a linear combination of the above, in each of the three weights

    x <- lapply(weight_first_sig_in_linear_comb, function(weight_first_sig_in_linear_comb_it){
      
      ## first bit of the signature
      first_sig <- sigs_cosmic0[,additional_signature][cumsum(sigs_cosmic0[,additional_signature]) < weight_first_sig_in_linear_comb_it]
      first_sig <- c(first_sig, weight_first_sig_in_linear_comb_it-sum(first_sig), rep(0, nrow(sigs_cosmic0)-1-length(first_sig)))
      second_sig <- sigs_cosmic0[,additional_signature]-first_sig
      ## scale them so that they add up to 1
      first_sig <- first_sig/sum(first_sig)
      second_sig <- second_sig/sum(second_sig)
      ## chech that they are indeed the linear combination
      stopifnot( (sum(first_sig*weight_first_sig_in_linear_comb_it+second_sig*(1-weight_first_sig_in_linear_comb_it)) - 1) < 0.02) ## allowing some error
      ## add these artificial signatures as the first two columns
      sigs_cosmic0 <- cbind(first_sig, second_sig, sigs_cosmic0)
      
      colnames(sigs_cosmic0)[1:2] <- c(paste0(additional_signature, '_part1'), paste0(additional_signature, '_part2'))
      
      extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                           subset_signatures = c(paste0(additional_signature, '_part1'),
                                                 paste0(additional_signature, '_part2'),
                                                 active_sigs_in_ct_list[[ct]][!(active_sigs_in_ct_list[[ct]] == additional_signature)]),
                           signature_version=NULL, signature_definition = sigs_cosmic0,
                           signature_fitting_method = signature_fitting_method)
    })
    names(x) <- paste0('_weightsig1incomb', weight_first_sig_in_linear_comb)
    x
  }, simplify = F)
}


re_run=FALSE
if(re_run){
  for(repl in 1:nreplicates){
    ## Run diagREDM
    dataset_original_signature_exposuresQP <- give_original_signature_exposures(signature_fitting_method = 'QP')
    dataset_linearlycombined_nocor_signature_exposuresQP <- give_linearlycombined_nocor_signature_exposures(signature_fitting_method = 'QP')
    saveRDS(dataset_original_signature_exposuresQP, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/1d_dataset_original_signature_exposuresQP_repl", repl, ".RDS"))
    saveRDS(dataset_linearlycombined_nocor_signature_exposuresQP, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/1d_dataset_linearlycombined_nocor_signature_exposuresQP_repl", repl, ".RDS"))
    
    diagREDM_original_signature_exposuresQP <- sapply(enough_samples, function(ct){
        try(CompSign::wrapper_run_TMB(object = dataset_original_signature_exposuresQP[[ct]],
                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F))
    }, simplify = F)
    
    diagREDM_linearlycombined_nocor_signature_exposuresQP <- sapply(enough_samples, function(ct){
      sapply(dataset_linearlycombined_nocor_signature_exposuresQP[[ct]], function(it_from_dataset){
        try(CompSign::wrapper_run_TMB(object = it_from_dataset,
                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F))
      }, simplify = F)
    }, simplify = F)
    
    saveRDS(diagREDM_original_signature_exposuresQP, file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1d/1d_dataset_original_signature_exposuresQP_repl", repl, ".RDS"))
    saveRDS(diagREDM_linearlycombined_nocor_signature_exposuresQP, file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1d/1d_dataset_linearlycombined_nocor_signature_exposuresQP_repl", repl, ".RDS"))
  }
}

diagREDM_original_signature_exposuresQP_all <- sapply(1:nreplicates, function(repl) readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1d/1d_dataset_original_signature_exposuresQP_repl", repl, ".RDS")), simplify=F)
diagREDM_linearlycombined_nocor_signature_exposuresQP_all <- sapply(1:nreplicates, function(repl) readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1d/1d_dataset_linearlycombined_nocor_signature_exposuresQP_repl", repl, ".RDS")), simplify=F)
for(i in 1:length(diagREDM_linearlycombined_nocor_signature_exposuresQP_all)){
  names(diagREDM_linearlycombined_nocor_signature_exposuresQP_all[[i]]) <- enough_samples
}

plotbettas_wrapper <- function(df, name_cats, nsigs_to_select=1){
  name_cats_noSBS = gsub("SBS", "", name_cats)
  plot_betas(df, return_plot = F, return_ggplot = T, plot = F, names_cats = name_cats, remove_SBS = T)+
    geom_hline(yintercept = c(-1, -2), lty='dashed')+geom_point(aes(x=factor(LogR, levels=name_cats_noSBS),
                                                                    y=`Estimate`,
                                                                    col=ifelse(LogR %in% name_cats_noSBS[1:nsigs_to_select], 'Split signature', 'Other signature'),
                                                                    shape=ifelse(LogR %in% name_cats_noSBS[1:nsigs_to_select], 'Split signature', 'Other signature')),
                                                                size=3)+
    labs(col='Signature in numerator', shape='Signature in numerator')+scale_color_manual(values = c('#0b30c9', '#d1b631'))+
    guides(col='none', shape='none')
}


re_run_plots=FALSE

if(re_run_plots){
  
  diagREDM_original_signature_exposuresQP <- readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1d/1d_dataset_original_signature_exposuresQP_repl", 1, ".RDS"))
  diagREDM_linearlycombined_nocor_signature_exposuresQP <- readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1d/1d_dataset_linearlycombined_nocor_signature_exposuresQP_repl", 1, ".RDS"))

  # ct='Bone-Osteosarc'
  for(ct in enough_samples){
    names_cats_true <- colnames(dataset_original_signature_exposuresQP[[ct]]$Y)
    names_cats_splitsig = CompSign::vector_cats_to_logR(c(paste0(names_cats_true[1], '*'), paste0(names_cats_true[1], '#'), names_cats_true[-1]))
    names_cats_true = CompSign::vector_cats_to_logR(names_cats_true)
    
    plt <- cowplot::plot_grid(plotbettas_wrapper(diagREDM_original_signature_exposuresQP[[ct]], names_cats_true)+ggtitle(ct),
                        plotbettas_wrapper(diagREDM_linearlycombined_nocor_signature_exposuresQP[[ct]]$`_weightsig1incomb0.75`, names_cats_splitsig, nsigs_to_select=2),
                        plotbettas_wrapper(diagREDM_linearlycombined_nocor_signature_exposuresQP[[ct]]$`_weightsig1incomb0.5`, names_cats_splitsig, nsigs_to_select=2),
                        plotbettas_wrapper(diagREDM_linearlycombined_nocor_signature_exposuresQP[[ct]]$`_weightsig1incomb0.25`, names_cats_splitsig, nsigs_to_select=2),
                        ncol = 1)
      
    pdf(paste0("../../../results/results_TMB/3_real_data_assessment/1d_robustness_of_results_linear_combination_signatures/betas/",
               "betas_", ct, ".pdf"), height = 10, width = 0.7*length(names_cats_true))
    print(plt)
    dev.off()
    png(paste0("../../../results/results_TMB/3_real_data_assessment/1d_robustness_of_results_linear_combination_signatures/betas/",
               "betas_", ct, ".png"), height = 10, width = 0.7*length(names_cats_true), units = 'in', res=300)
    print(plt)
    dev.off()
  }
}

## the structure we want is
## - ct
##   - parameter in 1,...length_out_iterations
##     - replicate
diagREDM_linearlycombined_nocor_signature_exposuresQP_all_new = list()
for(ct in enough_samples){
  diagREDM_linearlycombined_nocor_signature_exposuresQP_all_new[[ct]] <- sapply(weight_first_sig_in_linear_comb, function(j) sapply(1:nreplicates, function(repl) diagREDM_linearlycombined_nocor_signature_exposuresQP_all[[repl]][[ct]][[paste0('_weightsig1incomb', j)]], simplify = F), simplify = F)
  names(diagREDM_linearlycombined_nocor_signature_exposuresQP_all_new[[ct]]) <- paste0('_weightsig1incomb', weight_first_sig_in_linear_comb)
}
names(diagREDM_linearlycombined_nocor_signature_exposuresQP_all_new) <- enough_samples
diagREDM_linearlycombined_nocor_signature_exposuresQP_all = diagREDM_linearlycombined_nocor_signature_exposuresQP_all_new
diagREDM_linearlycombined_nocor_signature_exposuresQP_all_new = NULL

diagREDM_original_signature_exposuresQP_all <- sapply(enough_samples, function(ct) sapply(diagREDM_original_signature_exposuresQP_all, `[[`, ct, simplify = F), simplify = F, USE.NAMES = T)

plot_betacors_all_betas <- give_plot_betacors_nested(diagRE_DM = diagREDM_original_signature_exposuresQP_all,
                                                     fold_decrease_nmuts = 1:length(weight_first_sig_in_linear_comb),
                                                     diagREDM_varying_nmuts_with_replicates_signature_exposuresQP = diagREDM_linearlycombined_nocor_signature_exposuresQP_all,
                                                     mode='splitting_first_signature_all_betas', title_y = latex2exp::TeX(r"(Correlation of $\beta_0$ and $\beta_1$ with original dataset)"),
                                                     diagRE_DM_is_nested=T, level2_is_fold_decrease=F, gsub_from_x='_weightsig1incomb',
                                                     ncol_arg = 7)

plot_betacors_betaslope <- give_plot_betacors_nested(diagRE_DM = diagREDM_original_signature_exposuresQP_all,
                                                     fold_decrease_nmuts = 1:length(weight_first_sig_in_linear_comb),
                                                     diagREDM_varying_nmuts_with_replicates_signature_exposuresQP = diagREDM_linearlycombined_nocor_signature_exposuresQP_all,
                                                     mode='splitting_first_signature_beta_slope', title_y = latex2exp::TeX(r"(Correlation of $\beta_1$ with original dataset)"),
                                                     diagRE_DM_is_nested=T, level2_is_fold_decrease=F, gsub_from_x='_weightsig1incomb',
                                                     ncol_arg = 7)

plot_betacors_all_betas

# diagREDM_original_signature_exposuresQP$`Bone-Osteosarc`
# diagREDM_linearlycombined_nocor_signature_exposuresQP$`Bone-Osteosarc`$`_weightsig1incomb0.75`
# 
# cors <- give_beta_cor(df = diagREDM_linearlycombined_nocor_signature_exposuresQP,
#               diagRE_DM = diagREDM_original_signature_exposuresQP, mode = 'splitting_first_signature_all_betas')
# 
# plot_betacors_all_betas <- give_plot_betacors_notnested(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP = diagREDM_linearlycombined_nocor_signature_exposuresQP,
#                              diagRE_DM = diagREDM_original_signature_exposuresQP,
#                              mode = 'splitting_first_signature_all_betas',
#                              fold_decrease_nmuts = 1:length(weight_first_sig_in_linear_comb),
#                              level2_is_fold_decrease = F, gsub_from_x = "_weightsig1incomb",
#                              title_y = latex2exp::TeX(r"(Correlation of $\beta_0$ and $\beta_1$ with original dataset)"),
#                              ncol_arg = 7)
# plot_betacors_all_betas
# 
# plot_betacors_betaslope <- give_plot_betacors_notnested(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP = diagREDM_linearlycombined_nocor_signature_exposuresQP,
#                                                         diagRE_DM = diagREDM_original_signature_exposuresQP,
#                                                         mode = 'splitting_first_signature_beta_slope',
#                                                         fold_decrease_nmuts = 1:length(weight_first_sig_in_linear_comb),
#                                                         level2_is_fold_decrease = F, gsub_from_x = "_weightsig1incomb",
#                                                         title_y = latex2exp::TeX(r"(Correlation of $\beta_1$ with original dataset)"),
#                                                         ncol_arg = 7)

# 'intact_beta_slope', title_y = latex2exp::TeX(r"(Correlation of $\beta_1$ with original dataset)"))

plot_betacors_all_betas
ggsave("../../../results/results_TMB/3_real_data_assessment/1d_robustness_of_results_linear_combination_signatures/1d_robustness_split_linearcorrsigs_plot_betacors_allbetas_QP.pdf", height = 3.6, width = 6)
ggsave("../../../results/results_TMB/3_real_data_assessment/1d_robustness_of_results_linear_combination_signatures/1d_robustness_split_linearcorrsigs_plot_betacors_allbetas_QP.png", height = 3.6, width = 6)

plot_betacors_betaslope
ggsave("../../../results/results_TMB/3_real_data_assessment/1d_robustness_of_results_linear_combination_signatures/1d_robustness_split_linearcorrsigs_plot_betacors_betaslope_QP.pdf", height = 3.6, width = 6)
ggsave("../../../results/results_TMB/3_real_data_assessment/1d_robustness_of_results_linear_combination_signatures/1d_robustness_split_linearcorrsigs_plot_betacors_betaslope_QP.png", height = 3.6, width = 6)


