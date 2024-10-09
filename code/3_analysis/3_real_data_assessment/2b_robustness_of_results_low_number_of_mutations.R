## Is the overdispersion in the subclonal group merely a reflection of the lower number of mutations in this group?

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
##-----------------------------------------------------------------------------------------------------##

sort(sapply(enough_samples, function(ct) sum(read_info_list[[ct]]$dataset_nucleotidesubstitution3$Y)))

# ct <- 'CNS-PiloAstro'

lapply(enough_samples, function(ct, n = 6){
  ## progressively remove mutations from the sample
  
  dataset_it = read_info_list[[ct]]$dataset_nucleotidesubstitution3
  samples_to_remove_per_it <- round(c(sum(dataset_it$Y)*0.1))
  progressive_removal_mutations <- list()
  
  for(i in 1:n){
    counts_to_remove = sample(1:sum(dataset_it$Y), size = samples_to_remove_per_it, replace = F)
    
    dataset_it_prev = dataset_it
    ## remove these counts
    dataset_it = give_missclassified_matrix(original_matrix = dataset_it,
                                                  counts_to_misclassify = counts_to_remove, missclassify = F)
    stopifnot( (sum(dataset_it_prev$Y) == (sum(dataset_it$Y)+samples_to_remove_per_it)))
    progressive_removal_mutations[[i]] = dataset_it
  }
  
  active_sigs_in_ct <- colnames(read_info_list[[ct]]$dataset_active_sigs$Y)
  
  ## extract signatures
  signatures <- lapply(progressive_removal_mutations, function(i){
    try({extract_sigs_TMB_obj(dataset_obj_trinucleotide=i,
                         subset_signatures = active_sigs_in_ct,
                         signature_version=NULL, signature_definition = sigs_cosmic0,
                         signature_fitting_method = 'QP')
      })
  })
  
  signatures_DA <- lapply(signatures, function(i){
    if(class(i) == "try-error"){
      NA
    }else{
      CompSign::wrapper_run_TMB(object = i, model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
    }
    })
  
  
  saveRDS(signatures, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/2_robustness_of_results_varying_number_of_mutations/2b_lowering_number_of_mutations/signatures_", ct, ".RDS"))
  saveRDS(signatures_DA, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/2_robustness_of_results_varying_number_of_mutations/2b_lowering_number_of_mutations/TMBres_", ct, ".RDS"))

})