##------------------------------------------------------------------------##
rm(list = ls())
setwd("~/Documents/other/PhD/DM/CompSign-results/code/3_analysis/2_simulation_model_assessment/3_Comparison_other_models/")

library(CompSign)
library(HiLDA)
library(ggplot2)
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
get_inference_files <- function(folder_in){
  in_files <- unlist(sapply(grep('GenerationMixtureSimulation', 
                     list.files(folder_in, full.names = T),
                     value = T),
                list.files, full.names = T))
  
  ## there are the results from inference and the time it took for these models to run - select only the results
  in_files <- grep('.RDS', in_files, value = T)
  
  all_out_TMB <- lapply(in_files, readRDS)
  names(all_out_TMB) <- paste0(gsub(".RDS", "", 
                                    gsub("multiple_GenerationMixtureSimulation_", "", basename(in_files))))
  names(all_out_TMB) <-  gsub("diagREDM_NA_NA_NA_dataset", "", names(all_out_TMB))
  names(all_out_TMB) <-  gsub("TCSM_NA_NA_NA_dataset", "", names(all_out_TMB))
  names(all_out_TMB) <-  gsub("HilDA_NA_NA_NA_dataset", "", names(all_out_TMB))
  return(all_out_TMB)
}

get_dataset_files <- function(folder_in){
  in_files <- unlist(sapply(grep('GenerationMixtureSimulation', 
                                 list.files(folder_in, full.names = T),
                                 value = T),
                            list.files, full.names = T))
  
  ## there are the results from inference and the time it took for these models to run - select only the results
  in_files <- grep('.RDS', in_files, value = T)
  
  all_out_TMB <- lapply(in_files, readRDS)
  names(all_out_TMB) <- paste0(gsub(".RDS", "",  gsub("multiple_GenerationMixtureSimulation_", "",
                                                      basename(in_files))))
  names(all_out_TMB) <-  gsub("diagREDM_NA_NA_NA_dataset", "", names(all_out_TMB))
  names(all_out_TMB) <-  gsub("TCSM_NA_NA_NA_dataset", "", names(all_out_TMB))
  names(all_out_TMB) <-  gsub("HilDA_NA_NA_NA_dataset", "", names(all_out_TMB))
  return(all_out_TMB)
}
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
# all_in_datasets <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/datasets/")
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
## TMB
all_out_TMB <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/")
names(all_out_TMB)
pvals <- sapply(all_out_TMB, wald_TMB_wrapper)
max(pvals, na.rm = T)
min(pvals, na.rm = T)

df_pvals <- data.frame(DA_TMB=pvals <= 0.05, 
                       beta_TMB=sapply(as.numeric(sapply(names(pvals), function(i) strsplit(i, '_')[[1]][5])), function(i) softmax(c(i, 0))[1]),
                       betaALR_TMB=as.numeric(sapply(names(pvals), function(i) strsplit(i, '_')[[1]][5])),
                       numsamples=sapply(names(pvals), function(i) strsplit(i, '_')[[1]][1]),
                       numsmuts=sapply(names(pvals), function(i) strsplit(i, '_')[[1]][2]))
ggplot(df_pvals, aes(x=factor(betaALR_TMB, levels=sort(unique(df_pvals$betaALR_TMB))), fill=factor(DA_TMB)))+geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(.~interaction(numsamples, numsmuts))

ggplot(df_pvals, aes(x=factor(round(beta_TMB, 3), levels=sort(unique(round(beta_TMB, 3)))), fill=factor(DA_TMB)))+geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(.~interaction(numsamples, numsmuts))


##------------------------------------------------------------------------##
## HiLDA
all_out_HiLDA <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/HiLDA//")

# all_out_HiLDA[[1]]$
# hildaDiffPlots <- hildaDiffPlot(all_out_HiLDA[[1]], hildaLocal, sigOrder=c(1,3,2))
HiLDA::hildaGlobalResult(all_out_HiLDA[[1]])
# HiLDA::

hildaGlobalResult_v2 <- function (jagsOutput, pM1 = 0.5) {
  ## replacing pM2 by p2, as pM2 was not found in jagsOutput$BUGSoutput$sims.list
  if (is(jagsOutput) != "rjags") {
    stop("Not an output object from running the HiLDA tests.")
  }
  if (pM1 <= 0 | pM1 >= 1) {
    stop("Please input a fraction between 0 and 1.")
  }
  freq <- table(jagsOutput$BUGSoutput$sims.list$p2)
  if (length(freq) == 1) {
    stop("It got stuck in the model ", as.numeric(names(freq)) + 
           1)
  }
  return(as.vector(freq)[2]/as.vector(freq)[1] * (pM1)/(1 - 
                                                          pM1))
}

df_pvals_HiLDA <- sapply(all_out_HiLDA, hildaGlobalResult_v2)
table(df_pvals_HiLDA)
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
## TCSM
all_out_TCSM <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/TCSM/")
all_out_TCSM[[1]]$effect
all_out_TCSM[[1]]$sigma
## one value for each signature
as.numeric(all_out_TCSM[[1]]$EstSignificance$Mann.Whitney.U.p.value)
pvals_TCSM <- sapply(all_out_TCSM, function(j) min(as.numeric(j$EstSignificance$Mann.Whitney.U.p.value)))
df_pvals$DA_TCSM <- (pvals_TCSM <= 0.05)[match(names(pvals), names(pvals_TCSM))]
ggplot(df_pvals, aes(x=factor(round(beta_TMB, 3), levels=sort(unique(round(beta_TMB, 3)))), 
                     fill=factor(DA_TCSM)))+geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  facet_wrap(.~interaction(numsamples, numsmuts))

##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
plot(pvals[match(names(pvals_TCSM), names(pvals))], pvals_TCSM)
plot(pvals[match(names(pvals_TCSM), names(pvals))], pvals_TCSM)
table(TMB=df_pvals$DA_TMB, TCSM=df_pvals$DA_TCSM)

##------------------------------------------------------------------------##
