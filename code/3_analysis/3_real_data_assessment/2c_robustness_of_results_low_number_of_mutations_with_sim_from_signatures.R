##' Simulation to show how inference changes as fewer mutations are included in the PCAWG datasets.
##' The number of mutations is determined by drawing from a Poisson distribution. The same
##' per-subsample exposures, together with this scaling factor for the number of mutations, are
##' used to simulate trinucleotide abundances, and from those signature exposures are re-extracted
##' and the model is run.


##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
source("../../1_create_ROO/helper_1_create_ROO.R") ## for QPsig
source("../../3_analysis/2_simulation_model_assessment/1_generate_datasets_simulations/helper.R")
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
fold_decrease_nmuts <- c(1, 1.2, 1.4, 1.8, 2.5, 4)
nreplicates <- 5
signature_definitions = sigs_cosmic0
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##

give_signature_exposures_given_nmuts_with_replicates <- function(signature_fitting_method){
  sapply(enough_samples, function(ct){
    
    Npois = mean(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y))/fold_decrease_nmuts
    
    x <- lapply(Npois, function(Npois_it){
      
      lapply(1:nreplicates, function(repl_it){
        ## get ground truth signature proportions by normalising the observed data
        sigs_in_simplex <- normalise_rw(read_info_list[[ct]]$dataset_active_sigs$Y)
        
        ## get a number of mutations to simulate from (Npois), and from that a number of mutations per patient
        Npois_it
        Npois_it_vec <- rpois(n = nrow(read_info_list[[ct]]$dataset_active_sigs$Y), lambda = Npois_it)
        
        ## get ground truth signature exposures (in count)
        sigs_in_count <- round(sigs_in_simplex*Npois_it_vec)
        
        ## create trinucleotide count data from the ground truth signature exposures (in count)
        trinucleotides_matrix <- simulate_mutations_from_signatures(exposure_mat = sigs_in_count)
        
        ## re-extract signatures using QP or mutSigExtractor
        
        stopifnot(colnames(trinucleotides_matrix) == rownames(sigs_cosmic0))
        
        reextracted_sigs_in_count <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=list(Y=trinucleotides_matrix,
                                                                                         x=read_info_list[[ct]]$dataset_nucleotidesubstitution3$x,
                                                                                         z=read_info_list[[ct]]$dataset_nucleotidesubstitution3$z),
                                                          subset_signatures = colnames(read_info_list[[ct]]$dataset_active_sigs$Y),
                                                          signature_version=NULL, signature_definition = sigs_cosmic0,
                                                          signature_fitting_method = signature_fitting_method)
      })
    })
    names(x) <- paste0('_fold_decrease_nmuts', fold_decrease_nmuts)
    x
  }, simplify = F)
}

##-----------------------------------------------------------------------------------------------------##
## Same as above but with replicates (to account for e.g. how Npois varies across patients)

re_run <- FALSE

if(re_run){
  dataset_varying_nmuts_with_replicates_signature_exposuresQP <- give_signature_exposures_given_nmuts_with_replicates(signature_fitting_method = 'QP')
  saveRDS(dataset_varying_nmuts_with_replicates_signature_exposuresQP, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/2c_robustness_of_results_low_number_of_mutations_with_sim_from_signatures_with_replicates_QP.RDS"))
  
  diagREDM_varying_nmuts_with_replicates_signature_exposuresQP <- sapply(enough_samples, function(ct){
    sapply(dataset_varying_nmuts_with_replicates_signature_exposuresQP[[ct]], function(it_from_dataset){
      sapply(it_from_dataset, function(replicate_it){
        try(CompSign::wrapper_run_TMB(object = replicate_it,
                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F))
        }, simplify = F)
    }, simplify = F)
  }, simplify = F)
  saveRDS(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP, file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/2c_robustness_of_results_low_number_of_mutations_with_sim_from_signatures_with_replicates_diagREDM_QP.RDS"))
  
  fullREDM_varying_nmuts_with_replicates_signature_exposuresQP <- sapply(enough_samples, function(ct){
    sapply(dataset_varying_nmuts_with_replicates_signature_exposuresQP[[ct]], function(it_from_dataset){
      sapply(it_from_dataset, function(replicate_it){
        try(CompSign::wrapper_run_TMB(object = replicate_it,
                                      model = "fullRE_DM", use_nlminb=T, smart_init_vals=F))
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)
  saveRDS(fullREDM_varying_nmuts_with_replicates_signature_exposuresQP, file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/2c_robustness_of_results_low_number_of_mutations_with_sim_from_signatures_with_replicates_fullREDM_QP.RDS"))

}

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## see how the results from diagREDM deteriorate as the number of mutations used (Npois) decreases
diagREDM_varying_nmuts_with_replicates_signature_exposuresQP$`Bone-Osteosarc`$`_fold_decrease_nmuts1`
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
plot_betas(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP$`Bone-Osteosarc`$`_fold_decrease_nmuts1`[[1]], return_df = T)
plot_betas(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP$`Bone-Osteosarc`$`_fold_decrease_nmuts1`[[4]], return_df = T)
##-----------------------------------------------------------------------------------------------------##

plot_betacors_nested_allbetas_QP <- give_plot_betacors_nested(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP,
                                      diagRE_DM, mode='intact', title_y = latex2exp::TeX(r"(Correlation of $\beta_0$ and $\beta_1$ with original dataset)"))
plot_betacors_nested_beta1_QP <- give_plot_betacors_nested(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP,
                                  diagRE_DM, mode='intact_beta_slope', title_y = latex2exp::TeX(r"(Correlation of $\beta_1$ with original dataset)"))

plot_betacors_nested_allbetas_QP
ggsave("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2c_plot_betacors_nested_allbetas_QP.pdf", height = 3.6, width = 8)
ggsave("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2c_plot_betacors_nested_allbetas_QP.png", height = 3.6, width = 8)

plot_betacors_nested_beta1_QP
ggsave("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2c_plot_betacors_nested_beta1_QP.pdf", height = 3.6, width = 8)
ggsave("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2c_plot_betacors_nested_beta1_QP.png", height = 3.6, width = 8)


