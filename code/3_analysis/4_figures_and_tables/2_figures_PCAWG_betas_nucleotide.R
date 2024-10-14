## Comments on the cancer-specific results of PCAWG

##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../2_inference_TMB/")
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../3_analysis/4_figures_and_tables/helper_figures_manuscript.R")
source("../2_inference_TMB/helper_TMB.R")
# source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
# source("../3_analysis/recovery_COSMIC_signatures/recover_COSMIC_signatures.R")

library(TMB)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(dplyr)
library(jcolors)
library(viridis)
library(reshape2)
library(mutSigExtractor)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../data/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
enough_samples
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
read_info <- function(ct){
  .x <- list(#fullRE_M_SP = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREM_", ct, "_signaturesPCAWG.RDS"))),
             #fullRE_DMSL_SP = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDMsinglelambda_", ct, "_signaturesPCAWG.RDS"))),
             #fullRE_M_nonexo_SP = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullRE_DMSL_nonexo_SP = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDMsinglelambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             diagRE_DMDL_SP = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS"))),
             #diagRE_DMDL_nonexo_SP =  try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             #diagRE_DMDL_wSBS1SBS5nonexo_SP = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDMwSBS1SBS5nonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMnoscaling_SP_nonexo =  try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDMnoscalingnonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMnoscaling_SP_nonexo_subsets_and_amalgamations <- try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMonefixedlambdanonexo_SP = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMonefixedlambda2nonexo_SP = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDMonefixedlambda2nonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMonefixedlambdanonexo_SPSaA = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWGSaA.RDS"))),
             #fullREM_MSE = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREM_", ct, "_signaturesMSE.RDS"))),
             #fullREDM_MSE = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_signaturesMSE.RDS"))),
             fullREDM_nucleotide1 = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_nucleotidesubstitution1.RDS"))),
             #diagREDM_MSE = try(readRDS(paste0("../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDM_", ct, "_signaturesMSE.RDS"))),
             #dataset_all_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../data/", load_all_sigs = T, override_warning_X_Z = T),
             dataset_active_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../data/", override_warning_X_Z = T),
             dataset_nucleotidesubstitution1 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution1", path_to_data = "../../data/", override_warning_X_Z = T)
             #dataset_nucleotidesubstitution3 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution3", path_to_data = "../../data/", override_warning_X_Z = T),
             #dataset_nucleotidesubstitution3MSE = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution3MSE", path_to_data = "../../data/", override_warning_X_Z = T),
             #dataset_active_sigs_MSE = load_PCAWG(ct = ct, typedata = "signaturesMSE", path_to_data = "../../data/", load_all_sigs = F, override_warning_X_Z = T),
             #DMM = list(z_DMM=lapply(1:8, function(k) try(read.table(paste0("../../data/roo_for_DMM_SPpcawg/DMM_output/", ct, "_signaturesPCAWG_all", k, "_dmm.z"), sep = ',', skip = 1))),
            #            fit_DMM = lapply(1:8, function(k) try(read.table(paste0("../../data/roo_for_DMM_SPpcawg/DMM_output/", ct, "_signaturesPCAWG_all", k, "_dmm.fit"), sep = ' '))))
  )
  return(.x)
}
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
read_info_list <- lapply(enough_samples, read_info)
names(read_info_list) <- enough_samples
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../3_analysis/helper/pcawg.colour.palette.R")
pcawg_palette <- pcawg.colour.palette(x = gsub("\\..*", "", names(read_info_list)),  scheme = "tumour.subtype")
names(pcawg_palette) <- names(read_info_list)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
diagRE_DMDL <- lapply(read_info_list, function(i) i$diagRE_DMDL_SP)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
percentage_loss <- sapply(enough_samples, function(ct){
  x <- normalise_cl(sapply(split_matrix_in_half(read_info_list[[ct]]$dataset_nucleotidesubstitution1$Y), colSums))
  (x[2,]-x[1,])
})


sort(rowMeans(apply(percentage_loss, 2, function(i) i>0)))
sort(rowSums(apply(percentage_loss, 2, function(i) i>0)))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##

nucleotide1 <- sapply(read_info_list, `[`, 'fullREDM_nucleotide1')
names(nucleotide1) <- names(read_info_list)
names_trinucleotide <- vector_cats_to_logR(colnames(read_info_list[[1]]$dataset_nucleotidesubstitution1$Y))

give_plot_betas_all_ct(nucleotide1, names_trinucleotide)

# tikzDevice::tikz("../../results/results_TMB/pcawg/reports_per_cancer_type/cors_trinucleotide3sorted_v3.tex", height = 4, width = 5)
tikzDevice::tikz("../../results/results_TMB/pcawg/reports_per_cancer_type/cors_trinucleotide3sorted_v3.tex", height = 4, width = 5.5)
give_plot_betas_all_ct(nucleotide1, names_trinucleotide)
dev.off()

betas_nucleotides_slopes_softmax <- apply(betas_nucleotides_slopes, 2, function(i) softmax(c(i,0)))
rownames(betas_nucleotides_slopes_softmax) <- names(nucleotide_colours_dollar)
tikzDevice::tikz(file ="../../results/results_TMB/pcawg/betas_nucleotides_slopes_softmax_cor.tex",
                 width=3.2, height = 3)
pheatmap::pheatmap(cor(t(betas_nucleotides_slopes_softmax)))
dev.off()


## Correlations of coefficients

betas_nucleotides_slopes_mat = apply(betas_nucleotides_slopes, 2, as.vector)
rownames(betas_nucleotides_slopes_mat) <- gsub("/", "wrt", gsub("[$]|>", "", (rownames(betas_nucleotides_slopes))))

## Correlation of C>T with all others
(apply(betas_nucleotides_slopes_mat, 1, function(i) (cor(i, betas_nucleotides_slopes_mat['CTwrtTG',], method = 'pearson'))))

cor(betas_nucleotides_slopes_mat['TAwrtTG',],
     betas_nucleotides_slopes_mat['TCwrtTG',])
cor.test(betas_nucleotides_slopes_mat['TAwrtTG',],
    betas_nucleotides_slopes_mat['TCwrtTG',])

plot(betas_nucleotides_slopes_mat['TAwrtTG',],
     betas_nucleotides_slopes_mat['TCwrtTG',])

plot(betas_nucleotides_slopes_mat['TAwrtTG',],
     betas_nucleotides_slopes_mat['TCwrtTG',])

plot(betas_nucleotides_slopes_mat[4,],
     betas_nucleotides_slopes_mat[6,])


outer(1:nrow(betas_nucleotides_slopes_softmax), 1:nrow(betas_nucleotides_slopes_softmax), Vectorize(function(i,j){
  cor(betas_nucleotides_slopes_softmax[i,],
      betas_nucleotides_slopes_softmax[j,], method = 'pearson')
}))

cor(betas_nucleotides_slopes_softmax['T$>$A',],
    betas_nucleotides_slopes_softmax['T$>$C',], method = 'pearson')
        
cor.test(betas_nucleotides_slopes_softmax['T$>$A',],
         betas_nucleotides_slopes_softmax['T$>$C',], method = 'kendall')
