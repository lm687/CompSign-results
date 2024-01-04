rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../../2_inference_TMB/helper_TMB.R")
library(ggplot2)
library(dplyr)
library(reshape2)
library(latex2exp)
theme_set(theme_bw())

it=200

# basename <- "../../../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm_200_14072_80_6_0_diagREDM_betainterceptPCAWG4_betaslopezerosPCAWG4_covmatPCAWG4/multiple_GenerationJnorm_200_14072_80_6_0_diagREDM_betainterceptPCAWG4_betaslopezerosPCAWG4_covmatPCAWG4_dataset"
basename <- "../../../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm_200_14072_80_6_0_diagREDM_betainterceptPCAWG4_betaslopeonechangingPCAWG4_covmatPCAWG4/multiple_GenerationJnorm_200_14072_80_6_0_diagREDM_betainterceptPCAWG4_betaslopeonechangingPCAWG4_covmatPCAWG4_dataset"


res_inference_defaultbaseline <- lapply(0:(it-1), function(i) readRDS(paste0(basename, i, ".RDS")))
res_inference_maxbaseline <- lapply(0:(it-1), function(i) readRDS(paste0(gsub("GenerationJnorm", "GenerationJnormMax", basename), i, ".RDS")))
res_inference_maxbaseline

pvals_default <- sapply(res_inference_maxbaseline, wald_TMB_wrapper)
pvals_maxbaseline <- sapply(res_inference_maxbaseline, wald_TMB_wrapper)

plot(pvals_default, pvals_maxbaseline); abline(coef = c(0,1), col='blue')

pvals_maxbaseline == pvals_default

res_inference_defaultbaseline[[1]]
res_inference_maxbaseline[[1]]

get_softmax_df <- function(i){
  if(typeof(i) == 'character'){
    return(NA)
  }else{
    i = plot_betas(i, return_df = T) %>% select(-c('Std. Error')) %>% 
      dcast(LogR~type_beta, value.var = 'Estimate') %>% select(-'LogR')
    
    i = apply(i, 2, function(j) softmax(c(j,0)))
    ## change order of categories
    
    return(i)
  }
}

softmaxed_res_defaultbaseline <- lapply(res_inference_defaultbaseline, get_softmax_df)
softmaxed_res_maxbaseline <- lapply(res_inference_maxbaseline, get_softmax_df)
names(softmaxed_res_defaultbaseline) <- names(softmaxed_res_maxbaseline) <- paste0('It', 1:length(softmaxed_res_maxbaseline))

keep_runs <- (!is.na(softmaxed_res_defaultbaseline)) & (!is.na(softmaxed_res_maxbaseline))

plt1 <- ggplot(data.frame(pvals_default, pvals_maxbaseline), aes(x=pvals_default, y=pvals_maxbaseline))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+geom_point()+
  labs(x='p-values using default baseline', y='p-values using most abundant category as baseline')

# From invert_columns_dataset.R
# order_cols
# [1] 1 2 3 5 6 4
order_cols <- c(1, 2, 3, 5, 6, 4)

# softmaxed_res_defaultbaseline[keep_runs][[1]]
softmaxed_res_maxbaseline[keep_runs][[1]]
# softmaxed_res_maxbaseline[keep_runs][[1]][order_cols,]
softmaxed_res_defaultbaseline[keep_runs][[1]][order_cols,]

df_softmaxed <- list(default_baseline=do.call('rbind', lapply(softmaxed_res_defaultbaseline[keep_runs], function(i) cbind.data.frame(i[order_cols,], cat_idx=1:length(order_cols)))),
           max_baseline=do.call('rbind', lapply(softmaxed_res_maxbaseline[keep_runs], function(i) cbind.data.frame(i, cat_idx=1:length(order_cols)))))
# df_softmaxed <- list(default_baseline=do.call('rbind', lapply(softmaxed_res_defaultbaseline[keep_runs], function(i) cbind.data.frame(i, cat_idx=order_cols))),
#                      max_baseline=do.call('rbind', lapply(softmaxed_res_maxbaseline[keep_runs], function(i) cbind.data.frame(i, cat_idx=1:length(order_cols)))))

df_softmaxed <- melt(df_softmaxed, id.vars=c('cat_idx'))
df_softmaxed$Dataset <- rep(rep(paste0('Dataset', 1:sum(keep_runs)), each=length(order_cols)), 2*2)
df_softmaxed <- (dcast(df_softmaxed, cat_idx+variable+Dataset~L1, value.var = 'value'))

ggplot(df_softmaxed, aes(x=default_baseline, y=max_baseline, shape=factor(cat_idx), col=factor(cat_idx)))+geom_point()+
  facet_wrap(.~variable, scales = 'free')+
  geom_abline(slope = 1, intercept = 0)+
  scale_color_manual(values=c("#BDD9BF", "#929084", "#FFC857", "#A997DF", "#E5323B", "#2E4052"))

plt2 <- ggplot(df_softmaxed %>% filter( variable == 'Slope'),
               aes(x=default_baseline, y=max_baseline, shape=factor(cat_idx), col=factor(cat_idx)))+geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  scale_color_manual(values=c("#BDD9BF", "#929084", "#FFC857", "#A997DF", "#E5323B", "#2E4052"))+
  labs(x=TeX('Softmaxed \\beta_1 using default baseline'),
       y=TeX('Softmaxed \\beta_1 using most abundant category as baseline'), shape='Category', color='Category')
  # theme(legend.position = 'bottom')
plt2

cowplot::plot_grid(plt1, plt2, rel_widths = c(0.8, 1.2))
ggsave("../../../results/figures_paper/baseline_comparison.pdf", height = 4.5, width = 8)
ggsave("../../../results/figures_paper/baseline_comparison.png", height = 4.5, width = 8)
