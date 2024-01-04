## Comments on the cancer-specific results of PCAWG

##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../2_inference_TMB/")
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../2_inference_TMB/helper_TMB.R")
source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
source("../3_analysis/recovery_COSMIC_signatures/recover_COSMIC_signatures.R")

library(TMB)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(dplyr)
library(jcolors)
library(viridis)
library(mutSigExtractor)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
enough_samples

nucleotide_colours_logR <- c('C$>$A/T$>$G'= '#3cb371', 'C$>$G/T$>$G'= '#90ee90', 'C$>$T/T$>$G'= '#66cdaa',
                             'T$>$A/T$>$G'= '#cd5c5c', 'T$>$C/T$>$G'= '#f4a460')
nucleotide_colours <- c('C>A' = '#3cb371', 'C>G'= '#90ee90', 'C>T'= '#66cdaa',
                        'T>A'= '#cd5c5c', 'T>C'= '#f4a460', 'T>G'='red')
nucleotide_colours_dollar <- c('C$>$A' = '#3cb371', 'C$>$G'= '#90ee90', 'C$>$T'= '#66cdaa',
                               'T$>$A'= '#cd5c5c', 'T$>$C'= '#f4a460', 'T$>$G'='red')

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

diagRE_DMDL <- lapply(read_info_list, function(i) i$diagRE_DMDL_SP)

##-----------------------------------------------------------------------------------------------------##

nucleotide1 <- sapply(read_info_list, `[`, 'fullREDM_nucleotide1')
names(nucleotide1) <- names(read_info_list)
names_trinucleotide <- vector_cats_to_logR(colnames(read_info_list[[1]]$dataset_nucleotidesubstitution1$Y))

betas_nucleotides <- lapply(nucleotide1, function(i) plot_betas(i, return_df = T))
betas_nucleotides <- lapply(betas_nucleotides, function(i){
  i$LogR <- names_trinucleotide[i$LogR]
  # rownames(i) <- make.names(i$LogR, unique = T)
  i
})

betas_nucleotides_slopes <- do.call('cbind', lapply(betas_nucleotides, function(i) i%>% filter(type_beta == 'Slope' ) %>% select(Estimate)))
colnames(betas_nucleotides_slopes) <- names(nucleotide1)
rownames(betas_nucleotides_slopes) <- names_trinucleotide
rownames(betas_nucleotides_slopes) <- gsub(">", "$>$", rownames(betas_nucleotides_slopes))

tikzDevice::tikz("../../results/results_TMB/pcawg/reports_per_cancer_type/cors_trinucleotide3sorted_v3.tex", height = 4, width = 5)
ggplot(melt(as(betas_nucleotides_slopes, 'matrix')),
       aes(x=factor(Var2,levels=names(sort(colMeans(betas_nucleotides_slopes)))),
           col=Var1, y=value))+geom_point()+
  geom_hline(yintercept = 0, lty='dashed')+theme_bw()+geom_line(aes(group=Var1))+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank())+
  # labs(y=("$\\widehat{\\betab}_1$"))+
  labs(y=("$\\hat{\\beta}_1$"))+
  guides(col=guide_legend(nrow=1,byrow=TRUE))+
  scale_color_manual(values = nucleotide_colours_logR)
dev.off()

betas_nucleotides_slopes_softmax <- apply(betas_nucleotides_slopes, 2, function(i) softmax(c(i,0)))
rownames(betas_nucleotides_slopes_softmax) <- names(nucleotide_colours_dollar)
tikzDevice::tikz(file ="../../results/results_TMB/pcawg/betas_nucleotides_slopes_softmax_cor.tex",
                 width=3.2, height = 3)
pheatmap::pheatmap(cor(t(betas_nucleotides_slopes_softmax)))
dev.off()

