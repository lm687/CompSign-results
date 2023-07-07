## Comments on the cancer-specific results of PCAWG

##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
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
  .x <- list(#fullRE_M_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct, "_signaturesPCAWG.RDS"))),
             #fullRE_DMSL_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambda_", ct, "_signaturesPCAWG.RDS"))),
             #fullRE_M_nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullRE_DMSL_nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             diagRE_DMDL_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS"))),
             #diagRE_DMDL_nonexo_SP =  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             #diagRE_DMDL_wSBS1SBS5nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMwSBS1SBS5nonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMnoscaling_SP_nonexo =  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMnoscaling_SP_nonexo_subsets_and_amalgamations <- try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMonefixedlambdanonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMonefixedlambda2nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambda2nonexo_", ct, "_signaturesPCAWG.RDS"))),
             #fullREDMonefixedlambdanonexo_SPSaA = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWGSaA.RDS"))),
             #fullREM_MSE = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct, "_signaturesMSE.RDS"))),
             #fullREDM_MSE = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", ct, "_signaturesMSE.RDS"))),
             fullREDM_nucleotide1 = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", ct, "_nucleotidesubstitution1.RDS"))),
             #diagREDM_MSE = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_", ct, "_signaturesMSE.RDS"))),
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

##-----------------------------------------------------------------------------------------------------##

ovrdisp <- do.call('rbind.data.frame', lapply(1:length(diagRE_DMDL), try(function(idx){
  if(diagRE_DMDL[[idx]]$pdHess){
    cbind.data.frame( plot_lambdas(diagRE_DMDL[[idx]], return_df=T, plot=F), ct=names(diagRE_DMDL)[idx])
  }else{
    c(NA, NA)
  }
})))
ovrdisp[ovrdisp$name == 'Lambda 1','name'] = 'Clonal'
ovrdisp[ovrdisp$name == 'Lambda 2','name'] = 'Subclonal'
# ovrdisp$differentially_abundant = ifelse(ovrdisp$ct %in% names(differential_precision[(differential_precision <= 0.05)]), yes = '*', no = '')
differential_precision_2 <- p.adjust(sapply(diagRE_DMDL, ttest_TMB_wrapper_overdisp), method = 'fdr')
names(differential_precision_2) <- names(diagRE_DMDL)
which(differential_precision_2 <= 0.05) ## which are differentially overdispersed
length(which(differential_precision_2 <= 0.05))
cat(names(which(differential_precision_2 <= 0.05)), sep=', ')
head(dcast(ovrdisp, formula = ct~name, value.var=c('Estimate'))) ## in whih cases subclonal is more dispersed
table(apply(dcast(ovrdisp, formula = ct~name, value.var=c('Estimate')), 1, function(i) i[2] > i[3]))


ovrdisp$differential_precision_2 = ifelse(ovrdisp$ct %in% names(differential_precision_2[(differential_precision_2 <= 0.05)]), yes = '*', no = '')
# tikzDevice::tikz(file = "../../results/figures_paper/overdispersion_params_groups_test_2.tex", height = 3.5, width = 6.5)
ggplot(ovrdisp, aes(x=ct,  y=`Estimate`, group=name, col=name))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=`Estimate`-`Std..Error`,
                    ymax=`Estimate`+`Std..Error`), width=.1, position=position_dodge(width=0.5))+
  theme_bw()+
  geom_text(aes(y=Inf, label=differential_precision_2, vjust=1.8), col='black')+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='', y='Log lambda', col='Group')+theme(legend.position = "bottom")+
  theme(        legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(-10,-10,-10,-10),
                plot.margin = unit(c(1,1,1,1), "cm"))+
  scale_color_manual(values=c('#3b4d61', '#ef9d10'))
# dev.off()
## add manually in file: $\log(\widehat{\lambda})$


##-----------------------------------------------------------------------------------------------------##
all_diagREDMDL_betas <- lapply(enough_samples, function(ct){
  plot_betas(TMB_obj = read_info_list[[ct]]$diagRE_DMDL_SP,
             names_cats= vector_cats_to_logR(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)),
             return_df=T, plot=F, only_slope = T, line_zero=F, add_confint = T)
})
all_diagREDMDL_betas <- lapply(all_diagREDMDL_betas, function(i) cbind(i, numerator_LogR = gsub("/.*", "", i$LogR)))


all_diagREDMDL_betas_softmax <- lapply(all_diagREDMDL_betas,
                                       function(i) cbind.data.frame(Estimate=softmax(c(i$Estimate[i$type_beta == 'Slope'], 0)),
                                                                    sig=c(i$numerator_LogR[i$type_beta == 'Slope'],
                                                                          strsplit(i$LogR[1], '/')[[1]][2])))
names_sigs_unique <- gtools::mixedsort(unique(do.call('rbind', all_diagREDMDL_betas_softmax)$sig))

all_diagREDMDL_betas_softmax_1_5_40 <- t(sapply(all_diagREDMDL_betas_softmax, function(i) i[match(names_sigs_unique, i$sig,),'Estimate']))
colnames(all_diagREDMDL_betas_softmax_1_5_40) <- paste0('SBS', names_sigs_unique)
all_diagREDMDL_betas_softmax_1_5_40 <- data.frame(all_diagREDMDL_betas_softmax_1_5_40, ct=enough_samples)
pairs(all_diagREDMDL_betas_softmax_1_5_40)

give_scatter_two_betas <- function(sig1, sig2){
  ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes_string(x=sig1, y=sig2, col='ct', group=1))+
    geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
    theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
    geom_smooth(method = "lm", size=0.2)+
    ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.38, vjust=1, hjust=-.02, size=3.5, label.sep='\n')+
    labs(x=paste0('Beta slope of ', sig1), y=paste0('Beta slope of ', sig2))+guides(col='none')
}

do.call('grid.arrange', c(list(
give_scatter_two_betas('SBS1', 'SBS5'),
give_scatter_two_betas('SBS5', 'SBS40'),
give_scatter_two_betas('SBS2', 'SBS13'),
give_scatter_two_betas('SBS2', 'SBS8'),
give_scatter_two_betas('SBS2', 'SBS18'),
give_scatter_two_betas('SBS3', 'SBS18')), nrow=2))

##-----------------------------------------------------------------------------------------------------##

increases_matrix <- outer(names_sigs_unique, names_sigs_unique, Vectorize(function(sig_it1, sig_it2){
  mean(as.numeric(sapply(all_diagREDMDL_betas_softmax, function(j) j[j$sig == sig_it1,'Estimate'] > j[j$sig == sig_it2,'Estimate'] )), 
       na.rm = T)
}))
num_ct_in_common_increases_matrix <- outer(names_sigs_unique, names_sigs_unique, Vectorize(function(sig_it1, sig_it2){
  sum(!is.na(as.numeric(sapply(all_diagREDMDL_betas_softmax, function(j) j[j$sig == sig_it1,'Estimate'] > j[j$sig == sig_it2,'Estimate'] ))))
}))
colnames(increases_matrix) <- rownames(increases_matrix) <- names_sigs_unique

pheatmap_increases_matrix <- pheatmap::pheatmap(increases_matrix)

increases_matrix_melt <- melt(increases_matrix)
increases_matrix_melt$size= melt(num_ct_in_common_increases_matrix)$value
increases_matrix_melt$Var1 = factor(increases_matrix_melt$Var1, levels=pheatmap_increases_matrix$tree_row$labels[pheatmap_increases_matrix$tree_row$order])
increases_matrix_melt$Var2 = factor(increases_matrix_melt$Var2, levels=pheatmap_increases_matrix$tree_col$labels[pheatmap_increases_matrix$tree_row$order])
increases_matrix_melt <- increases_matrix_melt[!is.na(increases_matrix_melt$value),]

increases_matrix_melt$Var1 <- factor(paste0('SBS', increases_matrix_melt$Var1), levels=paste0('SBS', levels(increases_matrix_melt$Var1)))
increases_matrix_melt$Var2 <- factor(paste0('SBS', increases_matrix_melt$Var2), levels=paste0('SBS', levels(increases_matrix_melt$Var2)))

cowplot::plot_grid(cowplot::plot_grid(plot.new(), ggplot(increases_matrix_melt,
                                                         aes(x=Var1, y=Var2, col=value,size=size))+geom_point()+theme_bw()+ # shape=(value>0.5), 
                                        scale_color_viridis() + #theme(legend.position = "left", legend.box="vertical")+
                                        scale_y_discrete(position = "right")+
                                        theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1), legend.position = "bottom")+
                                        labs(x='',y='', size='Number of signatures in common', col='Fraction of higher coefficients'),
                                      nrow=2,
                                      rel_heights = c(0.03, 0.9)),
                   plot_grid(plot.new(),ggdendro::ggdendrogram(pheatmap_increases_matrix$tree_row, rotate=T), plot.new(),
                             rel_heights = c(0.005, 0.9, 0.1), ncol=1), rel_widths = c(3,1))
# ggsave("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/comparison_order_coefficients_heatmap_with_tree_bottomlegend.pdf",
#        height = 8, width = 12)
