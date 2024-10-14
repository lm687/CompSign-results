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

# > num_muts_clonal_subclonal
# Bone-Osteosarc Breast-AdenoCA CNS-GBM CNS-Medullo CNS-PiloAstro ColoRect-AdenoCA Eso-AdenoCA Head-SCC Kidney-ChRCC Kidney-RCC.clearcell
# [1,]          91626         591762  174061       94652          5280          2735108     1052644   355753        37560               432350
# [2,]          37003         362322  237861       53933          6756          2965529      564454   120722        34088               146562
# Kidney-RCC.papillary Liver-HCC Lung-SCC Lymph-BNHL Lymph-CLL Ovary-AdenoCA Panc-AdenoCA Panc-Endocrine Prost-AdenoCA Skin-Melanoma.cutaneous
# [1,]               133055   2028626  1000076     360857     84613        628571       844225         118502        519160                 2256247
# [2,]                31400    596333   247050      84572     46831        236328       458903         114581        317790                  383657
# Stomach-AdenoCA Thy-AdenoCA Uterus-AdenoCA
# [1,]          449787       33974         484817
# [2,]          410947       26275         208795

nmuts_to_remove_per_iteration_list <- list()
nmuts_to_remove_per_iteration_list[['Breast-AdenoCA']] <- 50000
nmuts_to_remove_per_iteration_list[['Thy-AdenoCA']] <- 2000

ct <- 'Thy-AdenoCA' # 'Breast-AdenoCA'

## progressively remove mutations from the clonal group, so that the number of mutations in the clonal group ends by being lower than in the subclonal


clonal_submatrix0 = split_matrix_in_half(read_info_list[[ct]]$dataset_nucleotidesubstitution3$Y)[[1]]
clonal_submatrix = clonal_submatrix0
subclonal_submatrix = split_matrix_in_half(read_info_list[[ct]]$dataset_nucleotidesubstitution3$Y)[[2]]
sum(clonal_submatrix)
sum(subclonal_submatrix)

nmuts_to_remove_per_iteration <- nmuts_to_remove_per_iteration_list[[ct]]
n=6

sum(clonal_submatrix0)-(nmuts_to_remove_per_iteration*1:n)
sum(subclonal_submatrix)

progressive_removal_mutations <- list()
for(i in 1:n){
  ## subclonal_submatrix is left untouched: subclonal_submatrix
  
  ## in clonal_submatrix we remove mutations, one step at a time. We select indices at random, which is equivalent to
  ## removing a percentage of mutations from each patient 
  ## (so that we preserve the variability in the number of mutations in each patient)
  patient_clonal_nmuts = rowSums(clonal_submatrix)
  clonal_submatrix
  
  counts_to_remove = sample(1:sum(clonal_submatrix), size = nmuts_to_remove_per_iteration, replace = F)
  
  clonal_submatrix_prev = clonal_submatrix
  clonal_submatrix = give_missclassified_matrix(original_matrix = list(Y=clonal_submatrix_prev),
                                                counts_to_misclassify = counts_to_remove, missclassify = F)$Y
  stopifnot( (sum(clonal_submatrix_prev)-sum(clonal_submatrix)) == nmuts_to_remove_per_iteration)
  cat('Number of mutations in clonal group: ', sum(clonal_submatrix), '\n')
  progressive_removal_mutations[[i]] = rbind(clonal_submatrix, subclonal_submatrix)
}

active_sigs_in_ct <- colnames(read_info_list[[ct]]$dataset_active_sigs$Y)

signatures <- lapply(progressive_removal_mutations, function(i){
  extract_sigs_TMB_obj(dataset_obj_trinucleotide=list(x=read_info_list[[ct]]$dataset_nucleotidesubstitution3$x,
                                                      z=read_info_list[[ct]]$dataset_nucleotidesubstitution3$z,
                                                      Y=i),
                     subset_signatures = active_sigs_in_ct,
                     signature_version=NULL, signature_definition = sigs_cosmic0,
                     signature_fitting_method = 'QP')
})
sapply(signatures, function(i) sum(i$Y)) ## decreasing

saveRDS(signatures, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/2_robustness_of_results_varying_number_of_mutations/overdispersion/signatures_", ct, ".RDS"))


signatures_DA <- lapply(signatures, function(i) CompSign::wrapper_run_TMB(object = i,
                          model = "diagRE_DM", use_nlminb=T, smart_init_vals=F))

saveRDS(signatures_DA, file = paste0("../../../data/assessing_models_real_data/simulated_datasets/2_robustness_of_results_varying_number_of_mutations/overdispersion/TMBres_", ct, ".RDS"))

betas_cor <- outer(signatures_DA, signatures_DA, Vectorize(function(i,j){
  cor(plot_betas(i, return_df = T)$Estimate,
      plot_betas(j, return_df = T)$Estimate)
}))
stopifnot(min(betas_cor)> 0.95)
min(betas_cor) ## all betas are identical

signatures_DA[[1]]
signatures_DA[[2]]

lambda = sapply(signatures_DA, function(i) python_like_select_name(i$par.fixed, 'log_lambda'))
colnames(lambda) = paste0(round((sum(clonal_submatrix0)-(nmuts_to_remove_per_iteration*1:n))/sum(subclonal_submatrix), digits = 2), 'X')
rownames(lambda) = c('Clonal', 'Subclonal')

head(melt((lambda)))

ggplot(melt(lambda),
       aes(x=Var2, y=value, col=Var1, group=Var1))+geom_line()+
  labs(x=paste0('Size of clonal group with respect to subclonal'), y=latex2exp::TeX(r"($\ln(\hat{\lambda})$)"), col='Group')+
  scale_color_manual(values=c('#3b4d61', '#ef9d10'))+
  theme(legend.position = "bottom")+ggtitle(paste0('Overdispersion as clonal group is\nsynthetically reduced: ', ct))
ggsave(paste0("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2a_robustness_of_results_varying_number_of_mutations_", ct, ".pdf"),
       height = 4, width = 4)
ggsave(paste0("../../../results/results_TMB/3_real_data_assessment/2_robustness_of_results_low_number_of_mutations/2a_robustness_of_results_varying_number_of_mutations_", ct, ".png"),
       height = 4, width = 4)

