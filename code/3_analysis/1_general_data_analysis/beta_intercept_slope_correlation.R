rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(gridExtra)

source("../../2_inference_TMB/helper_TMB.R")

### correlation between beta shape and intercept, in simulated datasets
flder_in <- "../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries_multiple/"
all_files_TMB <- list.files("../../../data/assessing_models_simulation/inference_results/TMB/nlminb/", full.names = T)

plot_smooth <- function(i, title=''){
  ggplot(i, aes(x=beta_est, y=beta_intercept_est, group=idx))+geom_point()+
    # stat_smooth (geom="line", alpha=0.3, size=3)+
    geom_smooth(method = "lm", se = F)+
    # geom_line(stat="smooth",method = "lm", alpha=0.8, col='blue')+
    theme_bw()+ggtitle(title)
}


generations <- c("GenerationJnorm", ## no intercept
                 "GenerationK",
                 "GenerationK2",
                 "generationFnorm",
                 "GenerationCnorm",
                 "GenerationJnorm2",
                 "GenerationJnorm3",
                 "GenerationJnormTwoLambdas",
                 "GenerationInoREtwolambdas",
                 "generationHnormtwolambdas",
                 "GenerationMixture1",
                 "GenerationMixturePCAWG",
                 "GenerationMixturefewersignaturesPCAWG",
                 "GenerationJnormTwoLambdasOneChangingBeta",
                 "GenerationJnormBTwoLambdasOneChangingBeta",
                 'GenerationMixturefewersignaturespairedPCAWG', "GenerationMixturefewersignaturespairedKidneyRCCPCAWG",
                 'GenerationMixturefewersignaturespairedstomachPCAWG', 'GenerationMixturefewersmallsignaturespairedKidneyRCCPCAWG',
                 'GenerationMixturefewersignaturespairedProstAdenoCAPCAWG', 'GenerationMixturefewersignaturespairedCNSGBMPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmPancEndocrinePCAWG', 'GenerationMixturefewersignaturespairedObsNmUterusAdenoCAPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGCNSGBMPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmPoissonSigPCAWGColoRectAdenoCAPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmInvPCAWGColoRectAdenoCAPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmInvPCAWGCNSGBMPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmObsDMColoRectAdenoCAPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmObsDMCNSGBMPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmObsDMEsoAdenoCAPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmObsDMHeadSCCPCAWG',
                 'GenerationMixturefewersignaturespairedObsNmObsDMKidneyChRCCPCAWG',
                 'GenerationMixtureallsignaturespairedObsNmCNSGBMPCAWG')

for(generation in generations){
  
  # runs_fullREM0 = readRDS(paste0(flder_in, generation, "_fullREM.RDS"))
  # runs_fullREDMSL0 = readRDS(paste0(flder_in, generation, "_fullREDMsinglelambda.RDS"))
  # runs_diagREDMSL0 = readRDS(paste0(flder_in, generation, "_diagREDMsinglelambda.RDS"))
  # runs_diagREDM0 = readRDS(paste0(flder_in, generation, "_diagREDM.RDS"))
  
  file_RDS <- paste0("/Users/morril01/Documents/PhD/GlobalDA/results/results_TMB/simulated_datasets/correlation_beta_intercept_slope/betas_df_", generation, ".RDS")
  if(file.exists(file_RDS)){
    betas_per_model <- readRDS(file_RDS)
  }else{
    files_generation <- all_files_TMB[grepl(generation, all_files_TMB)]
    names_models <- c('diagREDMsinglelambda', 'fullREDMsinglelambda', 'fullREM', 'diagREDM')
    files_generation_by_model <- lapply(names_models, function(i) files_generation[grepl(i, files_generation)])
    
    load_get_betas <- function(i){
      try({ .x <- readRDS(i)
      .betas <- give_betas(.x)
      data.frame(beta_intercept_est=.betas[1,], beta_est=.betas[2,], idx=basename(i)) })
    }
    
    betas_per_model <- lapply(files_generation_by_model, function(j) do.call('rbind', lapply(j, load_get_betas)))
    names(betas_per_model) <- names_models
    saveRDS(betas_per_model, file_RDS)
  }
  runs_fullREM0 <- betas_per_model$fullREM
  runs_fullREDMSL0 <- betas_per_model$fullREDMsinglelambda
  runs_diagREDMSL0 <- betas_per_model$diagREDMsinglelambda
  runs_diagREDM0 <- betas_per_model$diagREDM

  runs_fullREM0$beta_intercept_est <- as.numeric(runs_fullREM0$beta_intercept_est)
  runs_fullREM0$beta_est <- as.numeric(runs_fullREM0$beta_est)
  runs_fullREDMSL0$beta_intercept_est <- as.numeric(runs_fullREDMSL0$beta_intercept_est)
  runs_fullREDMSL0$beta_est <- as.numeric(runs_fullREDMSL0$beta_est)
  runs_diagREDMSL0$beta_intercept_est <- as.numeric(runs_diagREDMSL0$beta_intercept_est)
  runs_diagREDMSL0$beta_est <- as.numeric(runs_diagREDMSL0$beta_est)
  runs_diagREDM0$beta_intercept_est <- as.numeric(runs_diagREDM0$beta_intercept_est)
  runs_diagREDM0$beta_est <- as.numeric(runs_diagREDM0$beta_est)
  
  pdf(paste0("/Users/morril01/Documents/PhD/GlobalDA/results/results_TMB/simulated_datasets/correlation_beta_intercept_slope/", generation, ".pdf"))
  grid.arrange(plot_smooth(runs_fullREM0, title='fullREM'),
               plot_smooth(runs_fullREDMSL0, title='fullREDMSL'),
               plot_smooth(runs_diagREDMSL0, title='diagREDMSL'),
               plot_smooth(runs_diagREDM0, title='diagREDM'), top=generation)
  dev.off()
}

