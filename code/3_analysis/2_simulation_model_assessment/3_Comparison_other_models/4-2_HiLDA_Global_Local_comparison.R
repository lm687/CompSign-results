##------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(CompSign)
library(HiLDA)
library(ggplot2)
library(lsa) ## cosine similarity
library(reshape2)
library(dplyr)
library(GGally)
library(gridExtra)
theme_set(theme_bw())
source("helper.R")
source("../../../2_inference_TMB/helper_TMB.R")


all_out_HiLDA <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/HiLDA/",
                                     remove_HiLDAGlobal = F)

names(all_out_HiLDA) <- sapply(names(all_out_HiLDA), rename_datasets_fun)
## some are HiLDA and some are HiLDAGlobal. Make sure they are the same and then only keep the HiLDA files to avoid
## duplication.
all_out_HiLDA_global <- grepl('HiLDAGlobal', names(all_out_HiLDA))
all_out_HiLDA_global <- all_out_HiLDA[all_out_HiLDA_global]
all_out_HiLDA_matchglobal <- all_out_HiLDA[gsub("HiLDAGlobal", "HiLDA", names(all_out_HiLDA_global))]
all_out_HiLDA_global <- all_out_HiLDA_global[!(sapply(all_out_HiLDA_matchglobal, length) == 0)]
all_out_HiLDA_matchglobal <- all_out_HiLDA_matchglobal[!(sapply(all_out_HiLDA_matchglobal, length) == 0)]
comparison_global_HiLDA = do.call('rbind', lapply(1:length(all_out_HiLDA_global), function(idx) data.frame(run=idx, HiLDAGlobal=head(all_out_HiLDA_global[[idx]]$BUGSoutput$summary, n=18)[,1], HiLDA=head(all_out_HiLDA_matchglobal[[idx]]$BUGSoutput$summary, n=18)[,1])))
comparison_global_HiLDA$param_general = gsub("\\[.*", "", rownames(comparison_global_HiLDA))
ggplot(comparison_global_HiLDA, aes(x=HiLDAGlobal, y=HiLDA))+geom_point()+facet_wrap(.~param_general, scales = 'free')

##------------------------------------------------------------------------##

subset_10kit <- grep('10kit', names(all_out_HiLDA), value = T)
length(subset_10kit)

length(subset_10kit)

subset_10kit_list <- list(HiLDA_10kit = all_out_HiLDA[subset_10kit],
                          HiLDA = all_out_HiLDA[gsub("HiLDA10kit", "HiLDA", subset_10kit)],
                          HiLDA_Global = all_out_HiLDA[gsub("HiLDA10kit", "HiLDAGlobal", subset_10kit)])
sapply(subset_10kit_list, length)

HiLDA::hildaLocalResult(subset_10kit_list$HiLDA_10kit$`TwoCT_100_50_NA_NA_-1_HiLDA10kit_NA_NA_NA_dataset0`)
HiLDA::hildaLocalResult(subset_10kit_list$HiLDA_10kit$`TwoCT_100_50_NA_NA_1_HiLDA10kit_NA_NA_NA_dataset0`)

sapply(subset_10kit_list$HiLDA_10kit, hildaGlobal_from_alpha)

ex1 <- list(python_like_select_rownames(subset_10kit_list$HiLDA_10kit[[1]]$BUGSoutput$summary, 'alpha*'),
            python_like_select_rownames(subset_10kit_list$HiLDA[[1]]$BUGSoutput$summary, 'alpha*'),
            python_like_select_rownames(subset_10kit_list$HiLDA_Global[[1]]$BUGSoutput$summary, 'alpha*'))

pairs(sapply(ex1, function(i) i[,1]))
