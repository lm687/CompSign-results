setwd("~/Documents/other/PhD/DM/CompSign-results/code/3_analysis/2_simulation_model_assessment/3_Comparison_other_models/")

library(CompSign)
library(ggplot2)

in_files <- unlist(sapply(grep('GenerationMixtureSimulation', 
                   list.files("../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/", full.names = T),
                   value = T),
              list.files, full.names = T))

## there are the results from inference and the time it took for these models to run - select only the results
in_files <- grep('.RDS', in_files, value = T)

all_out_TMB <- lapply(in_files, readRDS)
names(all_out_TMB) <- paste0(gsub(".RDS", "", gsub("diagREDM_NA_NA_NA_dataset", "", 
                                                            gsub("multiple_GenerationMixtureSimulation_", "", basename(in_files)))))


pvals <- sapply(all_out_TMB, wald_TMB_wrapper)
max(pvals, na.rm = T)
min(pvals, na.rm = T)

df_pvals <- data.frame(DA=pvals <= 0.5, 
                       beta=sapply(as.numeric(sapply(names(pvals), function(i) strsplit(i, '_')[[1]][5])), function(i) softmax(c(i, 0))[1]),
                       betaALR=as.numeric(sapply(names(pvals), function(i) strsplit(i, '_')[[1]][5])),
                       numsamples=sapply(names(pvals), function(i) strsplit(i, '_')[[1]][1]),
                       numsmuts=sapply(names(pvals), function(i) strsplit(i, '_')[[1]][2]))
ggplot(df_pvals, aes(x=factor(betaALR, levels=sort(unique(df_pvals$betaALR))), fill=factor(DA)))+geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(.~interaction(numsamples, numsmuts))

ggplot(df_pvals, aes(x=factor(round(beta, 3), levels=sort(unique(round(beta, 3)))), fill=factor(DA)))+geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(.~interaction(numsamples, numsmuts))

