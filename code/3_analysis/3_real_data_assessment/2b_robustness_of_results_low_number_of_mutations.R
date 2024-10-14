##' Simulation to show how inference changes as fewer mutations are included in the PCAWG datasets.
##' A percentage of mutations is removed iteratively.

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
source("../../2_inference_TMB/helper_TMB.R")
source("../../3_analysis/helper/pcawg.colour.palette.R") ## pcawg_palette
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../../data/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
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
## From the PCAWG dataset, as shown in the paper
diagRE_DM <- sapply(read_info_list, function(i) i$diagRE_DMDL_SP, simplify = F)
diagRE_DM_tests <-  sapply(diagRE_DM, CompSign::wald_TMB_wrapper, simplify = F)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
length_out_iterations <- 6 ## 90% mutations, 80%, ...., 40%
nreplicates <- 5
##-----------------------------------------------------------------------------------------------------##

re_run <- FALSE

if(re_run){
  sort(sapply(enough_samples, function(ct) sum(read_info_list[[ct]]$dataset_nucleotidesubstitution3$Y)))
  

  for(repl in 1:nreplicates){
    lapply(enough_samples, function(ct, n = length_out_iterations){
      ## progressively remove mutations from the sample
      
      outTMB <- paste0("../../../data/assessing_models_real_data/simulated_datasets/2_robustness_of_results_varying_number_of_mutations/2b_lowering_number_of_mutations/TMBres_", ct, "_QP_repl", repl, ".RDS")
      
      if(file.exists(outTMB)){
        print('File exists\n')
      }else{
        
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
            try({CompSign::wrapper_run_TMB(object = i, model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)})
          }
          })
        
        
        saveRDS(signatures, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/2_robustness_of_results_varying_number_of_mutations/2b_lowering_number_of_mutations/signatures_", ct, "_QP_repl", repl, ".RDS"))
        saveRDS(signatures_DA, file = outTMB)
      
      }
    })
  }
}

# signatures <- readRDS(file = paste0("../../../data/assessing_models_real_data/simulated_datasets/2_robustness_of_results_varying_number_of_mutations/2b_lowering_number_of_mutations/signatures_", ct, ".RDS"))
signatures_DA <- sapply(enough_samples, function(ct) sapply(1:nreplicates, function(repl){
    readRDS(file = paste0("../../../data/assessing_models_real_data/simulated_datasets/2_robustness_of_results_varying_number_of_mutations/2b_lowering_number_of_mutations/TMBres_", ct, "_QP_repl", repl, ".RDS"))
  }, simplify=F), simplify=F)

signatures_DA

## so far the structure is
## - ct
##   - replicate
##      - parameter in 1,...length_out_iterations
## whereas what we want is
## - ct
##   - parameter in 1,...length_out_iterations
##     - replicate
## Make this change:
for(ct in 1:length(signatures_DA)){
  signatures_DA[[ct]] <- sapply(1:length_out_iterations, function(i) sapply(1:length(signatures_DA[[ct]]), function(j) signatures_DA[[ct]][[j]][[i]], simplify = F), simplify = F)
}

plot_betacors_nested_allbetas_QP <- give_plot_betacors_nested(signatures_DA, fold_decrease_nmuts = 1:length_out_iterations,
                                                              diagRE_DM, mode='intact', title_y = latex2exp::TeX(r"(Correlation of $\beta_0$ and $\beta_1$ with original dataset)"))

plot_betacors_nested_beta1_QP <- give_plot_betacors_nested(signatures_DA, fold_decrease_nmuts = 1:length_out_iterations,
                                                                 diagRE_DM, mode='intact_beta_slope', title_y = latex2exp::TeX(r"(Correlation of $\beta_1$ with original dataset)"))


# plot_betacors_notnested_allbetas_QP <- give_plot_betacors_notnested(signatures_DA, fold_decrease_nmuts = 1:length_out_iterations,
#                                                                     diagRE_DM, mode='intact', title_y = latex2exp::TeX(r"(Correlation of $\beta_0$ and $\beta_1$ with original dataset)"))
# 
# plot_betacors_notnested_beta1_QP <- give_plot_betacors_notnested(signatures_DA, fold_decrease_nmuts = 1:length_out_iterations,
#                                                                  diagRE_DM, mode='intact_beta_slope', title_y = latex2exp::TeX(r"(Correlation of $\beta_1$ with original dataset)"))


plot_betacors_nested_allbetas_QP
ggsave("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2b_plot_betacors_allbetas_QP.pdf", height = 3.6, width = 8)
ggsave("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2b_plot_betacors_allbetas_QP.png", height = 3.6, width = 8)

plot_betacors_nested_beta1_QP
ggsave("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2b_plot_betacors_beta1_QP.pdf", height = 3.6, width = 8)
ggsave("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2b_plot_betacors_beta1_QP.png", height = 3.6, width = 8)


