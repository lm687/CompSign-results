## In each cancer type, add signatures SBS1, SBS5, and SBS40

## for each cancer type:
# SBS1
# SBS5
# SBS40
# SBS1 and SBS5
# SBS5 and SBS40
# SBS1 and SBS40
# SBS1, SBS5 and SBS40

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

##-----------------------------------------------------------------------------------------------------##
signature_vec <- c('SBS1', 'SBS5', 'SBS40')
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
rm_empty <- function(i) i[i!=""]
## combination of presence and absence of each of these three signatures
combinations_flat_and_clock_signature_exposures <- apply(expand.grid(c('', 'SBS1'), c('', 'SBS5'), c('', 'SBS40')), 1, paste0, collapse='_')
combinations_flat_and_clock_signature_exposures <- sapply(combinations_flat_and_clock_signature_exposures,
                                                          function(i) (strsplit(i, '_')[[1]]))
combinations_flat_and_clock_signature_exposures
names_removal_flat_clock_signatures <- sapply(combinations_flat_and_clock_signature_exposures, function(i) paste0(paste0('-', paste0(rm_empty(i))), collapse = ''))
##-----------------------------------------------------------------------------------------------------##


re_run <- FALSE
if(re_run){
  
  read_info_list_cp = read_info_list

  ##-----------------------------------------------------------------------------------------------------##
  ## and note which of those correspond to the actual set of active signatures
  table(sapply(enough_samples, function(ct) rev(colnames(read_info_list_cp[[ct]]$dataset_active_sigs$Y))[1] %in% signature_vec))
  
  ## first of all, check if there is any ct for which SSB1, SBS5 or SBS40 are the baseline category
  sapply(enough_samples, function(ct) colnames(read_info_list_cp[[ct]]$dataset_active_sigs$Y))
  ## for several cancer types, one of these 3 signatures SBS1', 'SBS5', 'SBS40' is the last one.
  
  ## Re-arrange signatures so that neither of the three are used as baseline. Put them all in the beginning
  
  for(ct in enough_samples){
    idx_flat_and_clock_sigs <- sapply(signature_vec, function(i) which(colnames(read_info_list_cp[[ct]]$dataset_active_sigs$Y) == i))
    if(any(sapply(idx_flat_and_clock_sigs, length) == 0)){
      ## some of these signatures is not found
      idx_flat_and_clock_sigs <- unlist(idx_flat_and_clock_sigs[sapply(idx_flat_and_clock_sigs, length) == 1])
    }
    read_info_list_cp[[ct]]$dataset_active_sigs$Y <- cbind(read_info_list_cp[[ct]]$dataset_active_sigs$Y[,idx_flat_and_clock_sigs],
                                                           read_info_list_cp[[ct]]$dataset_active_sigs$Y[,-idx_flat_and_clock_sigs])
  }
  
  sapply(enough_samples, function(ct) colnames(read_info_list_cp[[ct]]$dataset_active_sigs$Y))
  ##-----------------------------------------------------------------------------------------------------##
  
  ##-----------------------------------------------------------------------------------------------------##
  give_remove_flat_and_clock_signature_exposures <- function(signature_fitting_method){
    sapply(enough_samples, function(ct){
      active_sigs_in_ct <- colnames(read_info_list_cp[[ct]]$dataset_active_sigs$Y)
      x <- lapply(combinations_flat_and_clock_signature_exposures, function(flat_clock_signatures_to_remove){
        try(extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list_cp[[ct]]$dataset_nucleotidesubstitution3,
                             subset_signatures = active_sigs_in_ct[!(active_sigs_in_ct %in% flat_clock_signatures_to_remove)],
                             signature_version=NULL, signature_definition = sigs_cosmic0,
                             signature_fitting_method = signature_fitting_method))
      })
      names(x) <- names_removal_flat_clock_signatures
      x
    }, simplify = F)
  }
  
  datasets_remove_flat_and_clock_signature_exposures <- give_remove_flat_and_clock_signature_exposures(signature_fitting_method = 'QP')
  saveRDS(datasets_remove_flat_and_clock_signature_exposures, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/1c_datasets_remove_flat_and_clock_signature_exposuresQP.RDS"))
  
  df = datasets_remove_flat_and_clock_signature_exposures
  diagREDM_remove_flat_and_clock_signature_exposures <- sapply(enough_samples, function(ct){
    sapply(df[[ct]], function(it_from_dataset){
      try(CompSign::wrapper_run_TMB(object = it_from_dataset,
                                    model = "diagRE_DM", use_nlminb=T, smart_init_vals=F))
    }, simplify = F)
  }, simplify = F)
  saveRDS(diagREDM_remove_flat_and_clock_signature_exposures, file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1c_datasets_remove_flat_and_clock_signature_exposuresQP.RDS"))
}

datasets_remove_flat_and_clock_signature_exposures <- readRDS(file = paste0("../../../data/assessing_models_real_data/simulated_datasets/1c_datasets_remove_flat_and_clock_signature_exposuresQP.RDS"))
diagREDM_remove_flat_and_clock_signature_exposures <- readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/1c_datasets_remove_flat_and_clock_signature_exposuresQP.RDS"))


## 1. compare DA
diagREDM_remove_flat_and_clock_signature_exposures_tests <- sapply(diagREDM_remove_flat_and_clock_signature_exposures,
                                                                   function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE)

give_opposite_sigs <- function(i){
  i2 <- strsplit(i, '-')[[1]]
  paste0(signature_vec[!(signature_vec %in% i2)], collapse = '+')
}


diagREDM_remove_flat_and_clock_signature_exposures_tests_df <- melt(sapply(diagREDM_remove_flat_and_clock_signature_exposures_tests, t, USE.NAMES = T, simplify = F))
diagREDM_remove_flat_and_clock_signature_exposures_tests_df$Var2 <- as.character(diagREDM_remove_flat_and_clock_signature_exposures_tests_df$Var2)
diagREDM_remove_flat_and_clock_signature_exposures_tests_df$Var2b <- sapply(diagREDM_remove_flat_and_clock_signature_exposures_tests_df$Var2, give_opposite_sigs)
diagREDM_remove_flat_and_clock_signature_exposures_tests_df$Var2b[diagREDM_remove_flat_and_clock_signature_exposures_tests_df$Var2b == ""] = 'None'
diagREDM_remove_flat_and_clock_signature_exposures_tests_df$Var2b <- factor(diagREDM_remove_flat_and_clock_signature_exposures_tests_df$Var2b,
                                                                            levels=c("None", "SBS1",  "SBS5", "SBS40",  "SBS1+SBS5",  "SBS1+SBS40",  "SBS5+SBS40", "SBS1+SBS5+SBS40"))
diagREDM_remove_flat_and_clock_signature_exposures_tests_df$pval_in_original_dataset <- diagRE_DM_tests[match(diagREDM_remove_flat_and_clock_signature_exposures_tests_df$L1, names(diagRE_DM_tests))]
diagREDM_remove_flat_and_clock_signature_exposures_tests_df <- diagREDM_remove_flat_and_clock_signature_exposures_tests_df %>%
  dplyr::mutate(agreement =  ((pval_in_original_dataset<0.05) & (value<0.05)) | ((pval_in_original_dataset>0.05) & (value>0.05)))
diagREDM_remove_flat_and_clock_signature_exposures_tests_df$agreement <- ifelse(diagREDM_remove_flat_and_clock_signature_exposures_tests_df$agreement,
                                                                                'Agreement', 'Disagreement')

df_true_sigs_present <- melt(sapply(read_info_list, function(x) paste0(signature_vec[(signature_vec %in% colnames(x$dataset_active_sigs$Y))], collapse = '+'), simplify = F)) 
colnames(df_true_sigs_present) <- c('Var2b', 'L1')

ggplot(diagREDM_remove_flat_and_clock_signature_exposures_tests_df,
       aes(x=Var2b, y=factor(L1, levels=rev(unique(L1)))))+
  geom_tile(aes(fill=agreement))+
  scale_fill_manual(values = c("Disagreement"="red", "Agreement"="#eaf7e9"), na.value="#fccfd1")+
  geom_tile(data = df_true_sigs_present, fill=NA, col='black')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'bottom',
        legend.title=element_blank())+
  labs(x='Signatures removed from the dataset', fill='Agreement in DA with original dataset', y='Cancer types')
ggsave("../../../results/results_TMB/3_real_data_assessment/1c_flat_sigs/1c_DAagreement.pdf", height = 5, width = 5)
ggsave("../../../results/results_TMB/3_real_data_assessment/1c_flat_sigs/1c_DAagreement.png", height = 5, width = 5)

give_barplot_agreement_in_DA(diagRE_DM_tests = diagRE_DM_tests,
                             diagDM_leave_one_out_exposures_tests = diagREDM_remove_flat_and_clock_signature_exposures_tests)

diagRE_DM_tests$`Lymph-CLL`
diagREDM_remove_flat_and_clock_signature_exposures_tests$`Lymph-CLL`

## 2. compare betas. To compare betas, we must know which signatures have been added (or removed) wrt the set of active signatures

do.call('grid.arrange',
        lapply(diagREDM_remove_flat_and_clock_signature_exposures$`Bone-Osteosarc`,
               plot_betas, return_ggplot = T))

