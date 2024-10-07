##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## plotting betas from all the cancer types
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
library(TMB)
library(ggplot2)
library(gtools)
library(scico)
library(latex2exp)
library(cowplot)
library(ggpubr)
library(reshape2)
library(dplyr)
library(gridExtra)
library(ggdendro)
library(writexl)
library(pheatmap)
theme_set(theme_bw())
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
             dataset_active_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../data/", override_warning_X_Z = T))
  
  .x
}
read_info_list <- lapply(enough_samples, function(ct){
  read_info(ct)
}); names(read_info_list) <- enough_samples
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
pcawg_palette <- pcawg.colour.palette(x = gsub("\\..*", "", names(read_info_list)),  scheme = "tumour.subtype")
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
# results_TMB_pcawg <- lapply(enough_samples, read_info)
# names(results_TMB_pcawg) <- enough_samples

give_betas_df <- function( ...){
  
  all_diagREDMDL_betas <- lapply(enough_samples, function(ct){
    tmb_res <- read_info(ct)
    plot_betas(TMB_obj = tmb_res$diagRE_DMDL_SP, names_cats= vector_cats_to_logR(colnames(tmb_res$dataset_active_sigs$Y)), ...)
  }); names(all_diagREDMDL_betas) <- enough_samples
  all_diagREDMDL_betas <- lapply(all_diagREDMDL_betas, function(i) cbind(i, numerator_LogR = gsub("/.*", "", i$LogR)))
  # all_diagREDMDL_betas_softmax <- lapply(all_diagREDMDL_betas,
  #                                        function(i) cbind.data.frame(Estimate=softmax(c(i$Estimate[i$type_beta == 'Slope'], 0)),
  #                                                                     sig=c(i$numerator_LogR[i$type_beta == 'Slope'],
  #                                                                           strsplit(i$LogR[1], '/')[[1]][2])))
  ##------------------------------- ##
  
  ##-----------------------------------------------------------------------------------------------------##
  all_diagREDMDL_betas_df <- melt(all_diagREDMDL_betas)
  all_diagREDMDL_betas_df$sig <- all_diagREDMDL_betas_df$numerator_LogR
  all_diagREDMDL_betas_df$denominator <- sapply(strsplit(all_diagREDMDL_betas_df$LogR, '/'), `[[`, 2)
  all_diagREDMDL_betas_df$sig <- paste0('SBS', all_diagREDMDL_betas_df$sig)
  all_diagREDMDL_betas_df$denominator <- paste0('SBS', all_diagREDMDL_betas_df$denominator)
  all_diagREDMDL_betas_df$sig <- factor(all_diagREDMDL_betas_df$sig, 
                                        levels=gtools::mixedsort(unique(all_diagREDMDL_betas_df$sig)))
  
  all_diagREDMDL_betas_df_with_estimate = all_diagREDMDL_betas_df
  all_diagREDMDL_betas_df_with_estimate = dcast(all_diagREDMDL_betas_df_with_estimate, 
                                                type_beta+LogR+numerator_LogR+L1+sig+denominator~variable) %>%
    dplyr::arrange(L1, (Estimate)) %>% group_by(L1) %>% mutate(idx_in_ct = 1:n())
  
  all_diagREDMDL_betas_df = all_diagREDMDL_betas_df %>% filter(variable == 'Estimate')
  all_diagREDMDL_betas_df = all_diagREDMDL_betas_df %>% 
    group_by(L1) %>% 
    mutate(value_minSBS1 = value-value[which(sig == 'SBS1')],
           value_minSBS5 = value-value[which(sig == 'SBS5')])
  all_diagREDMDL_betas_df = all_diagREDMDL_betas_df %>% 
    group_by(L1) %>% mutate(idx_in_ct = 1:n())
  return(list(all_diagREDMDL_betas_df = all_diagREDMDL_betas_df, 
              all_diagREDMDL_betas_df_with_estimate = all_diagREDMDL_betas_df_with_estimate))
}

df_slopes <- give_betas_df(return_df=T, plot=F, only_slope = T, line_zero=F, add_confint = T)
all_diagREDMDL_betas_df = df_slopes$all_diagREDMDL_betas_df
all_diagREDMDL_betas_df_with_estimate = df_slopes$all_diagREDMDL_betas_df_with_estimate

df_interceots <- give_betas_df(return_df=T, plot=F, only_slope = F, line_zero=F, add_confint = T)
all_diagREDMDL_betas_df_intercept = df_interceots$all_diagREDMDL_betas_df %>% filter(type_beta == 'Intercept')
all_diagREDMDL_betas_df_with_estimate_intercept = df_interceots$all_diagREDMDL_betas_df_with_estimate %>% filter(type_beta == 'Intercept')
all_diagREDMDL_betas_df_with_estimate_intercept = all_diagREDMDL_betas_df_with_estimate_intercept %>%
  group_by(L1) %>% mutate(idx_in_ct = 1:n()) %>% ungroup()

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
palette_name <- 'cork'

ggplot(all_diagREDMDL_betas_df, aes(y=L1, x=sig, col=value))+
  geom_tile(aes(fill=value_minSBS1))+
  geom_point(aes(col=value_minSBS5), size=1.9)+
  # geom_point()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='Signatures', y='Cancer types', 
       fill= TeX(r"($\beta_1-\beta_1^{SBS1}$)"),
       col=TeX(r"($\beta_1-\beta_1^{SBS5}$)"))+
  scale_color_scico(palette = palette_name)+
  scale_fill_scico(palette = palette_name)
ggsave("../../../results/results_TMB/pcawg/betas_comparing_with_SBS1_SBS5.pdf", height = 5, width = 10)
ggsave("../../../results/results_TMB/pcawg/betas_comparing_with_SBS1_SBS5.png", height = 5, width = 10)

ggplot(all_diagREDMDL_betas_df, aes(y=L1, x=idx_in_ct, col=value))+
  geom_tile(aes(fill=value_minSBS1))+
  geom_point(aes(col=value_minSBS5), size=4.5)+
  # geom_point()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='Active signatures', y='Cancer types', fill= 'Softmax', col='das')+
  scale_color_scico(palette = palette_name)+
  scale_fill_scico(palette = palette_name)+
  geom_text(aes(label=gsub('SBS', '', sig)), col='white', size=1.9, hjust=0.4)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
signature_to_col <- function(i){
  sapply(i, function(j){
    if(j == 1){
      'SBS1'
    }else if(j == 5){
      'SBS5'
    }else{
      'Other'
    }
  })
}
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
ggplot(all_diagREDMDL_betas_df_with_estimate, 
       aes(x=idx_in_ct, y=Estimate, col=signature_to_col(numerator_LogR)))+
  geom_line(aes(group=L1))+
  geom_point(aes(shape=signature_to_col(numerator_LogR)))+
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`))+
  scale_color_manual(values=c('black', '#ffba81', '#a3f6ab'))
ggsave("../../../results/results_TMB/pcawg/betas_comparing_with_SBS1_SBS5_line_stderror.png", height = 5, width = 10)


plot_row_sorted_betas <- function(df, return_legend=F, nudge_y_arg=0.3, isintercept=FALSE, ...){
  a <- ggplot(df, 
              aes(x=idx_in_ct, y=Estimate, col=signature_to_col(numerator_LogR),
                  shape=signature_to_col(numerator_LogR)))+
    geom_line(aes(group=L1), ...)+
    geom_point()+geom_label(aes(y =Estimate+`Std. Error`, label=numerator_LogR),
                            # direction="y", segment.color = 'transparent',
                                  nudge_y=nudge_y_arg,
                                  label.size = NA,
                            fill = alpha(c("white"),0.5))+
    geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`))+
    scale_color_manual(values=c('black', '#FFC107', '#C86518'))+
    facet_wrap(.~L1_rows, nrow=2)+
    facet_grid(.~L1, scales='free', space='free_x')+
    theme(legend.position = 'bottom')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  if(isintercept){
    a <- a + labs(y=TeX(r"($\hat{\beta_0}$)"), col='Signature', shape='Signature')
    
  }else{
    a <- a + labs(y=TeX(r"($\hat{\beta_1}$)"), col='Signature', shape='Signature')
  }
  if(return_legend){
    a
  }else{
    a+guides(col='none', shape='none')
  }
}

# all_diagREDMDL_betas_df_with_estimate$L1_rows = all_diagREDMDL_betas_df_with_estimate$L1 %in% unique(all_diagREDMDL_betas_df_with_estimate$L1)[1:12]
# png("../../../results/results_TMB/pcawg/betas_comparing_with_SBS1_SBS5_line_stderror_facets.png", height = 12, width = 20, units = 'cm', res = 200)
# plot_grid(plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate %>% filter(L1_rows)),
#           plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate %>% filter(!L1_rows)),
#           cowplot::get_legend(plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate %>% filter(L1_rows), return_legend = T)), rel_heights = c(1, 1, 0.2), nrow=3)
# dev.off()
.tab_rows <- table(all_diagREDMDL_betas_df_with_estimate$L1, c(rep(1:4, each=rep(round(nrow(all_diagREDMDL_betas_df_with_estimate)/4)))[-1]))
all_diagREDMDL_betas_df_with_estimate$L1_rows = apply(.tab_rows, 1, which.max)[match(all_diagREDMDL_betas_df_with_estimate$L1, names(apply(.tab_rows, 1, which.max)))]
all_diagREDMDL_betas_df_with_estimate_intercept$L1_rows = all_diagREDMDL_betas_df_with_estimate$L1_rows 
png("../../../results/results_TMB/pcawg/betas_comparing_with_SBS1_SBS5_line_stderror_facetsv2.png",
    height = 25, width = 32, units = 'cm', res = 200)
# pdf("../../../results/results_TMB/pcawg/betas_comparing_with_SBS1_SBS5_line_stderror_facetsv2.pdf",
#     height = 10.9375, width = 14)
plot_grid(plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate %>% filter(L1_rows == 1)),
          plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate %>% filter(L1_rows == 2)),
          plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate %>% filter(L1_rows == 3)),
          plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate %>% filter(L1_rows == 4)),
          cowplot::get_legend(plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate %>% filter(L1_rows == 1), return_legend = T)),
          rel_heights = c(1, 1, 1, 1, 0.2), nrow=5)
dev.off()

## same plot but for the intercept
# png("../../../results/results_TMB/pcawg/betas_comparing_with_SBS1_SBS5_line_stderror_facets_v2_intercept.png",
#     height = 25, width = 32, units = 'cm', res = 200)
pdf("../../../results/results_TMB/pcawg/betas_comparing_with_SBS1_SBS5_line_stderror_facetsv2_intercept.pdf",
    height = 10.9375, width = 14)
plot_grid(plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate_intercept %>% filter(L1_rows == 1), lty='dashed', col='blue', nudge_y_arg = 0.8, isintercept = TRUE),
          plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate_intercept %>% filter(L1_rows == 2), lty='dashed', col='blue', nudge_y_arg = 0.8, isintercept = TRUE),
          plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate_intercept %>% filter(L1_rows == 3), lty='dashed', col='blue', nudge_y_arg = 0.8, isintercept = TRUE),
          plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate_intercept %>% filter(L1_rows == 4), lty='dashed', col='blue', nudge_y_arg = 0.8, isintercept = TRUE),
          cowplot::get_legend(plot_row_sorted_betas(all_diagREDMDL_betas_df_with_estimate_intercept %>% filter(L1_rows == 1), return_legend = T)),
          rel_heights = c(1, 1, 1, 1, 0.2), nrow=5)
dev.off()
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
all_fullREDMDL <- lapply(enough_samples, function(ct){
  tmb_res <- read_info(ct)
  tmb_res$fullRE_DMDL_SP
}); names(all_fullREDMDL) <- enough_samples

all_diagREDMDL <- lapply(enough_samples, function(ct){
  tmb_res <- read_info(ct)
  tmb_res$diagRE_DMDL_SP
}); names(all_diagREDMDL) <- enough_samples

all_datasets_active <- lapply(enough_samples, function(ct){
  tmb_res <- read_info(ct)
  tmb_res$dataset_active_sigs
}); names(all_datasets_active) <- enough_samples


a <- cbind.data.frame(num_sigs=sapply(all_diagREDMDL, function(j) try(length(python_like_select_name(j$par.fixed, 'beta'))/2)),
                      pval=sapply(all_diagREDMDL, wald_TMB_wrapper),
                      num_samples=sapply(all_diagREDMDL, function(j) try(length(j$par.random)/(length(python_like_select_name(j$par.fixed, 'beta'))/2))),
                      num_sigs_full=sapply(all_fullREDMDL, function(j) try(length(python_like_select_name(j$par.fixed, 'beta'))/2)),
                      pval_full=sapply(all_fullREDMDL, wald_TMB_wrapper))
# View(a)

# matrix(all_fullREDMDL_betas[[1]]$diag.cov.random, ncol=length(python_like_select_name(all_fullREDMDL_betas[[1]]$par.fixed, 'beta'))/2)

read_info_all <- lapply(enough_samples, read_info); names(read_info_all) <- enough_samples
nrow(read_info_all$`CNS-GBM`$dataset_active_sigs$Y)/2
sapply(read_info_all, function(i) nrow(i$dataset_active_sigs$Y)/2)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
df_num_mutations_ct_group <- melt(lapply(read_info_all, function(i) lapply(split_matrix_in_half(i$dataset_active_sigs$Y), rowSums))) %>%
  mutate(Group=ifelse(L2 == 1, 'Clonal', 'Subclonal'))
table(df_num_mutations_ct_group$L2, df_num_mutations_ct_group$Group) ## should be diagonal

df_num_mutations_ct_group_table <- df_num_mutations_ct_group[,-2] %>% group_by(L1, Group) %>%
                                           summarise(min=min(value), max=max(value),
                                                     mean=mean(value), median=median(value))
head(df_num_mutations_ct_group_table)
colnames(df_num_mutations_ct_group_table) <- c('Cancer type', 'Group', 'Min',  'Max',  'Mean', "Median")
writexl::write_xlsx(df_num_mutations_ct_group_table, "../../../results/supplementary_data/number_of_mutations_per_ct_group.xlsx")
print(xtable::xtable(df_num_mutations_ct_group_table), include.rownames=FALSE)

ggplot(df_num_mutations_ct_group,
       aes(x=value, col=Group))+geom_density()+
  facet_wrap(.~L1, scales = "free", ncol=4)+
  scale_x_log10()+
  scale_color_manual(values=c('#3b4d61', '#ef9d10'))+
  geom_vline(xintercept = 180, lty='dashed')+labs(x='Number of mutations', y='Density')+
  theme(legend.position = 'bottom')
ggsave("../../../results/results_TMB/pcawg/nlambda_in_samples_density.png", 
       height = 7, width = 8)
ggsave("../../../results/results_TMB/pcawg/nlambda_in_samples_density.pdf", 
       height = 7, width = 8)

##-----------------------------------------------------------------------------------------------------##


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
# pairs(all_diagREDMDL_betas_softmax_1_5_40)

cbind(names(pcawg_palette), enough_samples)
stopifnot(length(names(pcawg_palette)) == length(enough_samples))
names(pcawg_palette) <- enough_samples

give_scatter_two_betas <- function(sig1, sig2){
  ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes_string(x=sig1, y=sig2, col='ct', group=1))+
    geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
    theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
    geom_smooth(method = "lm", size=0.2)+
    ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.38, vjust=1, hjust=-.02, size=3.5, label.sep='\n')+
    labs(x=paste0('Beta slope of ', sig1), y=paste0('Beta slope of ', sig2))+guides(col='none')
}


give_scatter_two_betas_MSE <- function(sig1, sig2){
  subsetct = !(is.na(all_diagREDMDL_betas_softmax_1_5_40[,sig1])) & !(is.na(all_diagREDMDL_betas_softmax_1_5_40[,sig2]))
  mse <- mean((all_diagREDMDL_betas_softmax_1_5_40[subsetct,sig1]-all_diagREDMDL_betas_softmax_1_5_40[subsetct,sig2])**2)
  ggplot(all_diagREDMDL_betas_softmax_1_5_40[subsetct,], aes_string(x=sig1, y=sig2, col='ct', group=1))+
    geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
    theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
    # geom_smooth(method = "lm", size=0.2)+
    annotate("text", 
             x = min(all_diagREDMDL_betas_softmax_1_5_40[subsetct,sig1], na.rm = T)+0.3*(max(all_diagREDMDL_betas_softmax_1_5_40[subsetct,sig1], na.rm = T)-min(all_diagREDMDL_betas_softmax_1_5_40[subsetct,sig1], na.rm = T)),
             y = max(all_diagREDMDL_betas_softmax_1_5_40[subsetct,sig2], na.rm = T)-0.15*(max(all_diagREDMDL_betas_softmax_1_5_40[subsetct,sig2], na.rm = T)-min(all_diagREDMDL_betas_softmax_1_5_40[subsetct,sig2], na.rm = T)),
             label = paste0("MSE=", signif(mse, 2)), size=3.5, )+
    labs(x=paste0('Beta slope of ', sig1), y=paste0('Beta slope of ', sig2))+guides(col='none')
}

do.call('grid.arrange', c(list(
  give_scatter_two_betas('SBS1', 'SBS5'),
  give_scatter_two_betas('SBS5', 'SBS40'),
  give_scatter_two_betas('SBS2', 'SBS13'),
  give_scatter_two_betas('SBS2', 'SBS8'),
  give_scatter_two_betas('SBS2', 'SBS18'),
  give_scatter_two_betas('SBS3', 'SBS18')), nrow=2))

tikzDevice::tikz("../../../results/figures_paper/betas_MSE_scatter.tex", height = 5)
do.call('grid.arrange', c(list(
  give_scatter_two_betas_MSE('SBS1', 'SBS5'),
  give_scatter_two_betas_MSE('SBS5', 'SBS40'),
  give_scatter_two_betas_MSE('SBS2', 'SBS13'),
  give_scatter_two_betas_MSE('SBS2', 'SBS8'),
  give_scatter_two_betas_MSE('SBS2', 'SBS18'),
  give_scatter_two_betas_MSE('SBS3', 'SBS18')), nrow=2))
dev.off()

tikzDevice::tikz("../../../results/figures_paper/betas_MSE_scatter_additional.tex", height = 2.5)
do.call('grid.arrange', c(list(
  give_scatter_two_betas_MSE('SBS3', 'SBS8'),
  give_scatter_two_betas_MSE('SBS5', 'SBS8'),
  give_scatter_two_betas_MSE('SBS5', 'SBS18')#,
), nrow=1))
dev.off()


##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
names_sigs_unique <- gtools::mixedsort(unique(do.call('rbind', all_diagREDMDL_betas_softmax)$sig))
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
                                                         aes(x=Var1, y=Var2, col=value, size=size, shape=value>0.5))+geom_point()+theme_bw()+ # shape=(value>0.5), 
                                        scale_color_scico(palette = 'broc') +
                                        scale_y_discrete(position = "right")+
                                        theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1), legend.position = "bottom")+
                                        guides(shape='none')+
                                        labs(x='',y='', size='Number of cancer types in common', col='Fraction of higher coefficients'),
                                      nrow=2,
                                      rel_heights = c(0.03, 0.9)),
                   plot_grid(plot.new(),ggdendro::ggdendrogram(pheatmap_increases_matrix$tree_row, rotate=T), plot.new(),
                             rel_heights = c(0.005, 0.9, 0.1), ncol=1), rel_widths = c(3,1))
ggsave("../../../results/results_TMB/pcawg/betas_comparison_order_coefficients_heatmap_with_tree_bottomlegend.pdf",
       height = 8, width = 12)
ggsave("../../../results/results_TMB/pcawg/betas_comparison_order_coefficients_heatmap_with_tree_bottomlegend.png", 
       height = 8, width = 12)

num_ct_with_comparison <- dcast(increases_matrix_melt[,c(1:2, 4)], Var1~Var2)
num_ct_with_comparison[is.na(num_ct_with_comparison)] <- 0
writexl::write_xlsx(num_ct_with_comparison, "../../../results/supplementary_data/number_of_cancer_types_signature_combination.xlsx")
write.table(num_ct_with_comparison, "../../../results/supplementary_data/number_of_cancer_types_signature_combination.tsv",
            sep = '\t', quote = F, row.names = F)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
all_diagREDMDL_betas0_and_beta1 <- lapply(lapply(enough_samples, function(ct){
  plot_betas(TMB_obj = read_info_list[[ct]]$diagRE_DMDL_SP,
             names_cats= vector_cats_to_logR(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)),
             return_df=T, plot=F, only_slope = F, line_zero=F, add_confint = T)
}), function(i) cbind(i, numerator_LogR = gsub("/.*", "", i$LogR)))
names(all_diagREDMDL_betas0_and_beta1) <- enough_samples

all_diagREDMDL_betas0_and_beta1 <- melt(all_diagREDMDL_betas0_and_beta1, id.vars=colnames(all_diagREDMDL_betas0_and_beta1[[1]]))

all_diagREDMDL_betas0_and_beta1 %>% filter(numerator_LogR %in% c(1,5), L1=='Lung-SCC')
ggplot(dcast(all_diagREDMDL_betas0_and_beta1 %>% filter(numerator_LogR %in% c(1,5)), type_beta+L1~numerator_LogR, value.var = 'Estimate'),
       aes(x=`1`, y=`5`, col=L1))+geom_point()+facet_wrap(.~type_beta)+geom_abline(slope = 1, intercept = 0)+
  scale_color_manual(values = pcawg_palette)

ggplot(dcast(all_diagREDMDL_betas0_and_beta1 %>% filter(numerator_LogR %in% c(2,13)), type_beta+L1~numerator_LogR, value.var = 'Estimate'),
       aes(x=`2`, y=`13`, col=L1))+geom_point()+facet_wrap(.~type_beta)+geom_abline(slope = 1, intercept = 0)+
  scale_color_manual(values = pcawg_palette)

##-----------------------------------------------------------------------------------------------------##
all_fullREDMDL_converged <- sapply(all_fullREDMDL, wald_TMB_wrapper)

cors_betas_fullREDM_diagREDM <- t(sapply(names(all_fullREDMDL_converged[!is.na(all_fullREDMDL_converged)]), function(ct){
  x1= give_betas(all_fullREDMDL[[ct]])
  x2 = give_betas(all_diagREDMDL[[ct]])
  c(cor_beta=cor(as.vector(x1), as.vector(x2)),
    cor_beta0 = cor(x1[1,] ,x2[1,]),
    cor_beta1 = cor(x1[2,] ,x2[2,]))
}))
cors_betas_fullREDM_diagREDM

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
all_fullREDMDL_converged
##-----------------------------------------------------------------------------------------------------##


# idx <- 2

plot_intercepts_from_mat <- function(RI_mat, order_x, colour_scheme_bool, order_both){
  require(scico)
  if(colour_scheme_bool){
    palette_name <- 'vik'
  }else{
    palette_name <- 'broc'
  }
  
  RI_mat <- #data.frame(RI_mat)#, idx=1:nrow(RI_mat))
    RI_mat <- (melt(as.matrix(RI_mat)))
  RI_mat$Var1 <- factor(RI_mat$Var1, levels=unique(RI_mat$Var1)[order_x])
  if(order_both){
    RI_mat$Var2 <- factor(RI_mat$Var2, levels=unique(RI_mat$Var2)[order_x])
  }
  ggplot(RI_mat, aes(x=Var2, y=Var1, fill=value))+
    geom_tile()+
    scale_fill_scico(palette = palette_name)+
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.y=element_blank())
}

get_intercepts_and_plot <- function(idx, plot_correlation_mat=F, order_both=F, return_df=F, use_empirical_cor=FALSE){
  if(!is.na(all_fullREDMDL_converged[idx])){
    cat('Converged run\n')
    converged <- T
  }else{
    cat('Not converged run\n')
    converged <- F
  }
  dmin1 <- length(python_like_select_name(all_fullREDMDL[[idx]]$par.fixed, 'beta'))/2
  
  random_intercepts_diag <- matrix(all_diagREDMDL[[idx]]$par.random, ncol=dmin1)
  
  if(converged & !use_empirical_cor){
    random_intercepts <- matrix(all_fullREDMDL[[idx]]$par.random, ncol=dmin1)
    cormat <- L_to_cov(python_like_select_name(all_fullREDMDL[[idx]]$par.fixed, 'cov_par_RE'),
             d = dmin1)
  }else{
    ## Use the empirical correlation, either because it has not converged or because we are forcing it with use_empirical_cor
    random_intercepts <- matrix(all_diagREDMDL[[idx]]$par.random, ncol=dmin1)
    cormat = cor(random_intercepts)
  }
  
  active_sigs <- colnames(all_datasets_active[[idx]]$Y)
  
  if(plot_correlation_mat){
    # plt1 <- plot(as.vector(cov(random_intercepts)),
    # as.vector(xx))
    # abline(coef = c(0,1))
    x <- pheatmap::pheatmap(cormat)
    rownames(cormat) <- colnames(cormat) <- rev(rev(active_sigs)[-1])#, '/', rev(active_sigs)[1])
    
    if(return_df){
      return(cormat)
    }else{
      plot_intercepts_from_mat(RI_mat = cormat, order_x = x$tree_row$order, colour_scheme_bool=converged, order_both=T)+
        ggtitle(paste0(enough_samples[idx]))+labs(fill= paste0('/', rev(active_sigs)[1]))
    }
  }else{
    
    x <- pheatmap::pheatmap(random_intercepts)
    # plt2 <- pheatmap::pheatmap(random_intercepts[x$tree_row$order,], cluster_cols = F, cluster_rows = F)
    # plt3 <- pheatmap::pheatmap(random_intercepts_diag[x$tree_row$order,], cluster_cols = F, cluster_rows = F)
    
    # colnames(random_intercepts) <- paste0(rev(rev(active_sigs)[-1]), '/', rev(active_sigs)[1])
    colnames(random_intercepts) <- rev(rev(active_sigs)[-1])#, '/', rev(active_sigs)[1])
    
    if(return_df){
      return(random_intercepts)
    }else{
      # RI_mat = random_intercepts; order_x = x$tree_row$order
      plot_intercepts_from_mat(RI_mat = random_intercepts, order_x = x$tree_row$order, colour_scheme_bool=converged, order_both=F)+
        ggtitle(paste0(enough_samples[idx]))+labs(fill= paste0('/', rev(active_sigs)[1]))
    }
  }
}

get_intercepts_and_plot(1)
get_intercepts_and_plot(2)
get_intercepts_and_plot(2, plot_correlation_mat = T, order_both=T)

all_random_intercepts_plot <- lapply(1:length(enough_samples), get_intercepts_and_plot, plot_correlation_mat=F); names(all_random_intercepts_df) <- enough_samples
all_random_intercepts_df <- lapply(1:length(enough_samples), get_intercepts_and_plot, plot_correlation_mat=F, return_df=T); names(all_random_intercepts_df) <- enough_samples
all_random_intercepts_plot_cor <- lapply(1:length(enough_samples), get_intercepts_and_plot, plot_correlation_mat = T, order_both=T)
all_random_intercepts_plot_cor_df <- lapply(1:length(enough_samples), get_intercepts_and_plot, plot_correlation_mat = T, order_both=T, return_df=T)

pdf("../../../results/results_TMB/pcawg/randomintercepts_mat.pdf", height = 13, width = 11)
do.call('grid.arrange', all_random_intercepts_plot)
dev.off()
pdf("../../../results/results_TMB/pcawg/randomintercepts_cormat.pdf", height = 13, width = 11)
do.call('grid.arrange', all_random_intercepts_plot_cor)
dev.off()

## In cases where both converged, is the empirical correlation similar to the estimated correlation?
ct_converged_fullREDM = names(all_fullREDMDL_converged[!is.na(all_fullREDMDL_converged)])
all_random_intercepts_df[ct_converged_fullREDM] ## this contains the estimated covariances, or the 
all_random_intercepts_df_empiricalcor <- lapply(1:length(enough_samples), get_intercepts_and_plot, plot_correlation_mat=F, return_df=T, use_empirical_cor=T); names(all_random_intercepts_df_empiricalcor) <- enough_samples

ct_conv = ct_converged_fullREDM[1]
cor_covariances = sapply(ct_converged_fullREDM, function(ct_conv){
  cor(as.vector(all_random_intercepts_df[[ct_conv]]),
      as.vector(all_random_intercepts_df_empiricalcor[[ct_conv]]))
})
sort(cor_covariances)
max(cor_covariances); min(cor_covariances); median(cor_covariances)

cor_covariances = sapply(ct_converged_fullREDM, function(ct_conv){
  cor(as.vector(all_random_intercepts_df[[ct_conv]]),
      as.vector(all_random_intercepts_df_empiricalcor[[ct_conv]]))
})
######################################################################

plot(all_random_intercepts_df$`Eso-AdenoCA`[,c('SBS2', 'SBS13')])
APOBEC <- lapply(all_random_intercepts_df, function(i) try(i[,c('SBS2', 'SBS13')]))
APOBEC <- APOBEC[!(sapply(APOBEC, typeof) == 'character')]
APOBEC <- dcast(melt(APOBEC), Var1+L1~Var2, value.var = 'value')
ggplot(APOBEC, aes(x=SBS2, y=SBS13))+geom_point()+facet_wrap(.~L1, scales = 'free')


compute_cor_random_intercepts_pairs <- function(sigs){
  #df <- lapply(all_random_intercepts_df, function(i) try(i[,c(sigs[1], sigs[2])]))
  #df = sapply(df, function(i) try(cor(i[,1], i[,2])))
  df = sapply(all_random_intercepts_plot_cor_df, function(i) try(i[sigs[1], sigs[2]]))
  names(df) <- enough_samples
  df
}

list_pairs_sigs <- list(c('SBS1', 'SBS5'),
                        c('SBS1', 'SBS3'),
                        c('SBS1', 'SBS4'),
                        c('SBS1', 'SBS40'),
                        c('SBS2', 'SBS13'),
                        c('SBS2', 'SBS8'),
                        c('SBS2', 'SBS18'))
cors_RI <- lapply(list_pairs_sigs,
                       function(pair_sig) data.frame(compute_cor_random_intercepts_pairs(pair_sig)))
cors_RI <- lapply(cors_RI, function(i) subset(i, !grepl('Error', i[,1])))
names(cors_RI) <- sapply(list_pairs_sigs, paste0, collapse='-')
cors_RI <- melt(lapply(cors_RI, as.matrix))
cors_RI$Var2 %>% unique
ggplot(cors_RI, aes(y=factor(Var1, levels=rev(unique(Var1))), 
                    x=factor(L1, levels = sapply(list_pairs_sigs, paste0, collapse='-')), 
                    fill=as.numeric(value)))+geom_tile()+
  scale_fill_scico(palette = 'vik')+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  geom_label(aes(label=round(as.numeric(value), 2)), fill=alpha('black', 0.2), col='white', label.size = NA)+
  labs(x='Signature pair', y='Cancer type', fill='Pearson cor.')
ggsave("../../../results/results_TMB/pcawg/randomintercepts_pairssigscors.pdf", height = 7, width = 7)

all_random_intercepts_plot[[which(enough_samples == 'Liver-HCC')]]
all_random_intercepts_plot[[which(enough_samples == 'Panc-AdenoCA')]]
all_random_intercepts_plot[[which(enough_samples == 'Stomach-AdenoCA')]]

