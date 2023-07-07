rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../../../2_inference_TMB/helper_TMB.R")

library(gridExtra)
library(ggrepel)
require(jcolors)
require(reshape2)
require(dplyr)

generation = c("GenerationCnorm",  "generationFnorm", "generationGnorm", "generationMGnorm")

all_runs = lapply(generation, function(i){
  ## these files are not there anymore
  # runs_fullREM = readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/summaries/", i, "_fullREM.RDS"))
  # runs_fullREDM = readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/summaries/", i, "_fullREDM.RDS"))
  return(list(fullREM=runs_fullREM, fullREDM=runs_fullREDM))
})
names(all_runs) = generation
all_runs_m = melt(all_runs, id.vars = c('beta_true', 'idx', 'd', 'n', 'beta_gamma_shape', 'beta_est', 'beta_stderr',
                                        'pvals_adj', 'DA_bool', 'idx_within_dataset', 'bool_zero_true_beta', 'converged'))
all_runs_m = all_runs_m[!(duplicated(all_runs_m[,!(colnames(all_runs_m) %in% c('idx', 'beta_true', 'beta_est', 'beta_stderr',
                                                                               'type', 'idx_within_dataset'))])),]
all_runs_m$type = paste0(all_runs_m$DA_bool, all_runs_m$pvals_adj < 0.05)
# joint_df_grouping_by_n = all_runs_m[sapply(unique(all_runs_m$idx), function(i) which(all_runs_m$idx == i)[1]),] %>%
all_runs_m_grouped_n = all_runs_m %>%
  group_by(n, L1, L2) %>%
  mutate(sensitivity=sum(type == 'TRUETRUE')/sum(type %in% c('TRUETRUE', 'TRUEFALSE')),
         specificity=sum(type == 'FALSEFALSE')/sum(type %in% c('FALSETRUE', 'FALSEFALSE')))
all_runs_m_grouped_d = all_runs_m %>%
  group_by(d, L1, L2) %>%
  mutate(sensitivity=sum(type == 'TRUETRUE')/sum(type %in% c('TRUETRUE', 'TRUEFALSE')),
         specificity=sum(type == 'FALSEFALSE')/sum(type %in% c('FALSETRUE', 'FALSEFALSE')))

cat('Line 35\n')

ggplot(all_runs_m_grouped_n, aes(x=1-specificity, y=sensitivity, col=n, shape=(L1), label=interaction(L1, L2)))+
  geom_point()+facet_wrap(.~L2)+geom_label_repel()+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Number of samples", shape="Simulation type")+
  theme_bw()
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/summary_meta/ROC_meta_analysis_simulations.pdf"),
       width=7, height = 5)

ggplot(all_runs_m_grouped_d, aes(x=1-specificity, y=sensitivity, col=d, shape=(L1), label=interaction(L1, L2)))+
  geom_point()+facet_wrap(.~L2)+geom_label_repel()+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Number of categories", shape="Simulation type")+
  theme_bw()
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/summary_meta/ROC_meta_analysis_simulations_d.pdf"),
       width=7, height = 5)

## changing level of significance
give_ROC_alpha = function(significance_levels){
  all_runs_m$type = paste0(all_runs_m$DA_bool, all_runs_m$pvals_adj < significance_levels)
  # joint_df_grouping_by_n = all_runs_m[sapply(unique(all_runs_m$idx), function(i) which(all_runs_m$idx == i)[1]),] %>%
  all_runs_m %>%
    group_by(L1, L2) %>%
    mutate(sensitivity=sum(type == 'TRUETRUE')/sum(type %in% c('TRUETRUE', 'TRUEFALSE')),
           specificity=sum(type == 'FALSEFALSE')/sum(type %in% c('FALSETRUE', 'FALSEFALSE')),
           significance=significance_levels)
}

all_runs_m_significances = do.call('rbind.data.frame', lapply(seq(from = 0, to = 0.5, length.out = 20), give_ROC_alpha))
table(all_runs_m_significances$significance)

ggplot(droplevels(all_runs_m_significances), aes(x=1-specificity, y=sensitivity, col=L1, shape=interaction(L1)))+
  # geom_point()+
  facet_wrap(.~L2, nrow=1, scales = "free")+
  # geom_line(aes(group=interaction(n, L1, L2, d)))+
  geom_point()+
  geom_line(aes(group=interaction(L1, L2)))+
  scale_color_jcolors("pal5")+labs(col = "Number of samples", shape="Simulation type")+
  theme_bw()
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/summary_meta/ROC_meta_analysis_simulations_varyingalpha.pdf"),
       width=7, height = 4)

