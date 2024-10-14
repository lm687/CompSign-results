rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(CompSign)
library(ggplot2)
library(reshape2)
library(pheatmap)
# data(package='CompSign')
##------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------##
# source("../../../../../../CompSign-do-not-modify_do_everything_online/R/DA_functions.R")
## Looking at the differential abundance of mutational signatures in chromosomes
data(ProstAdenoCA_chrom, package='CompSign')

matrices_exposures_active_sigs <- sapply(ProstAdenoCA_chrom$all_exposures_ct_active, t)

chroms <- rownames(matrices_exposures_active_sigs[[1]])

## we will unify chromosomes: if some is absent, we will add exposures of z
matrices_exposures_active_sigs <- sapply(matrices_exposures_active_sigs, function(i){
  rwni <- rownames(i)
  n_missing <- chroms %in% rwni
  if(any(!n_missing)){
    i <- rbind(i, t(sapply(1:sum(!n_missing), function(j) rep(0, ncol(i)))))
    rownames(i) <- c(rwni, chroms[!(chroms %in% rwni)])
  }
  i[chroms,]
}, simplify = F)
##------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------##
## random effects:
## there are 24 observations per patient (one for each chromosome). We use patient-specific intercepts
object_for_DA <- list(x=cbind(1, do.call('rbind', lapply(matrices_exposures_active_sigs, function(x) diag(nrow(x))))[,-1]), ## removing the first column: the first chromosome will become the baseline
                      z=do.call('rbind', lapply(1:length(chroms), function(unused) diag(length(matrices_exposures_active_sigs)))),
                      Y=do.call('rbind', matrices_exposures_active_sigs))
sapply(object_for_DA, dim)
##------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------##
re_run <- F
if(re_run){
  ## with random effects
  # chrom_diagREDM <- wrapper_run_TMB(model = 'diagRE_DM', object = object_for_DA, smart_init_vals = F)
  # saveRDS(chrom_diagREDM, "../../../../data/additional_use_cases/chromosomes.RDS")

  ## with a chromosome-specific lambda
  # chrom_FEREDM <- wrapper_run_TMB(model = 'FE_DM', object = object_for_DA, smart_init_vals = F)
  # saveRDS(chrom_diagREDM, "../../../../data/additional_use_cases/chromosomes_FE.RDS")
  
  ## with one single lambda shared across chromosomes
  chrom_FEREDMsinglelambda <- wrapper_run_TMB(model = 'FE_DM_singlelambda', object = object_for_DA, smart_init_vals = F)
  chrom_FEREDMsinglelambda
  saveRDS(chrom_FEREDMsinglelambda, "../../../../data/additional_use_cases/chromosomes_FE_singlelambda.RDS")
}

##------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------##

logR_coef_names = c('chr1', vector_cats_to_logR(c(unique(rownames(object_for_DA$Y))[-1], 'chr1'), sep = '-'))
logR_sigs = vector_cats_to_logR(colnames(object_for_DA$Y))

plt_chrom <- plot_betas(chrom_FEREDMsinglelambda, num_covariates = ncol(object_for_DA$x),
           names_cats = logR_sigs, labels_betas=logR_coef_names, keep_order_signatures = T, 
           keep_order_covariates = T, return_ggplot = T, return_plot = F, plot = F)
plt_chrom+labs(x='Signature log-ratios', y=latex2exp::TeX(r"($\hat{\beta}$)"))
ggsave("../../../../results/additional_use_cases/chromosomes_FE_singlelambda_betas.pdf", height = 8, width = 8.5)
ggsave("../../../../results/additional_use_cases/chromosomes_FE_singlelambda_betas.png", height = 8, width = 8.5)
##------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------##

betas_df <- plot_betas(chrom_FEREDMsinglelambda, num_covariates = ncol(object_for_DA$x),
           names_cats = logR_sigs, labels_betas=logR_coef_names, keep_order_signatures = T, 
           keep_order_covariates = T, return_df = T)
betas_mat = dcast(betas_df[,-2], value.var='Estimate', type_beta~LogR)
rownames(betas_mat) <- betas_mat[,1]; betas_mat[,1] <- NULL

betas_cor <- outer(1:nrow(betas_mat), 1:nrow(betas_mat), Vectorize(function(x,y){
  cor(t(betas_mat[x,]), t(betas_mat[y,]))
}))
rownames(betas_cor) <- colnames(betas_cor) <- rownames(betas_mat)

cor_plot <- pheatmap::pheatmap(betas_cor,   color = colorRampPalette(c("#5FFBF1", "#86A8E7", "#D16BA5"))(100))
cor_plot
##------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------##
colnames(betas_cor) <- rownames(betas_cor) <- rownames(betas_mat)
pdf("../../../../results/additional_use_cases/chromosomes_FE_singlelambda_betas_cor.pdf", 
    height = 5, width = 5.5)
print(cor_plot)
dev.off()

png("../../../../results/additional_use_cases/chromosomes_FE_singlelambda_betas_cor.png", 
    height = 5, width = 5.5, units = 'in', res=300)
print(cor_plot)
dev.off()


