## To get results we need both the datasets files and the inference results

local=F

if(local){
  rm(list = ls())
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  
  multiple_runs = T
  ## multiple runs
  # generation = "GenerationJnorm"
  # generation = "GenerationK"
  # generation = "GenerationK2"weightedaccuracy_with_d_boxplot.pdf
  # generation = "generationFnorm"
  # generation = "GenerationCnorm"
  # generation = "GenerationJnorm2"
  # generation = "GenerationJnorm3"
  # generation = "GenerationJnormTwoLambdas"
  # generation = "GenerationInoREtwolambdas"
  # generation = "generationHnormtwolambdas"
  # generation = "GenerationMixture1"
  # generation = "GenerationJnormTwoLambdasOneChangingBeta"
  # generation = "GenerationJnormBTwoLambdasOneChangingBeta"
  # generation = "GenerationMixturePCAWG"
  # generation = "GenerationMixturefewersignaturesPCAWG"
  # generation <- 'GenerationMixturefewersignaturespairedPCAWG'
  # generation = "GenerationMixturefewersignaturespairedKidneyRCCPCAWG"
  # generation <- 'GenerationMixturefewersignaturespairedstomachPCAWG'
  # generation <- 'GenerationMixturefewersmallsignaturespairedKidneyRCCPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedProstAdenoCAPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedCNSGBMPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmPancEndocrinePCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmUterusAdenoCAPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGCNSGBMPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmPoissonSigPCAWGColoRectAdenoCAPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmInvPCAWGColoRectAdenoCAPCAWG'
  # generation <- 'GenerationJnormBTwoLambdasOneChangingBeta'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmInvPCAWGCNSGBMPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmObsDMColoRectAdenoCAPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmObsDMCNSGBMPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmObsDMEsoAdenoCAPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmObsDMHeadSCCPCAWG'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmObsDMKidneyChRCCPCAWG'
  # generation <- 'GenerationMixtureallsignaturespairedObsNmCNSGBMPCAWG'
  # generation <- 'GenerationPois'
  # generation <- 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGProstAdenoCAPCAWG'
  generation <- 'GenerationJnormBTwoLambdasOneChangingBeta'
  ##########################################
  # multiple_runs = F
  ## single replicate
  # generation = "generationGnorm"
  # generation = "generationMGnorm"
  # generation = "GenerationInoRE"
  # generation = "GenerationCnorm"
}else{
  library(optparse)
  setwd("3_analysis/simulation_model_assessment/analyse_inference_simulations")
  multiple_runs <- T
  option_list = list(
    make_option(c("--input"), type="character", default=NA,
                help="_results_info.txt for the models that we want to use. Only useful for snakemake; they are not used here", metavar="character"),
    make_option(c("--generation"), type="character", default=NA,
                help="Generation of simulation. This is used", metavar="character")
  );
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  generation <- opt$generation
}

source("../../../2_inference_TMB/helper_TMB.R")
source("../../../1_create_ROO/roo_functions.R")
source("helper_model_assessment.R")

include_singleREDM <- T

library(grid)
library(gridExtra)
#library(ggpubr)
library(reshape2)
library(jcolors)
library(cowplot)
library(ggrepel)
library(dplyr)

####
warning('Note: despite the name, pvals_adj are not adjusted for multiple-testing')

require( tikzDevice )

if(multiple_runs){
  flder_out <- "../../../../results/results_TMB/simulated_datasets/mixed_effects_models_multiple/"
  flder_in <- "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries_multiple/"
}else{
  flder_out <- "../../../../results/results_TMB/simulated_datasets/mixed_effects_models/"
  flder_in <- "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries/"
}


manual = F
if(manual){
  names1 <- paste0(flder_in, generation, "_fullREM_manual.RDS")
  names2 <- paste0(flder_in, generation, "_fullREDMsinglelambda_manual.RDS")
  names3 <- paste0(flder_in, generation, "_diagREDMsinglelambda_manual.RDS")
  names4 <- paste0(flder_in, generation, "_diagREDM_manual.RDS")
  runs_fullREM0 = readRDS(names1)
  runs_fullREDMSL0 = readRDS(names2)
  runs_diagREDMSL0 = readRDS(names3)
  runs_diagREDM0 = readRDS(names3)
}else{
  runs_fullREM0 = readRDS(paste0(flder_in, generation, "_fullREM.RDS"))
  runs_fullREDMSL0 = readRDS(paste0(flder_in, generation, "_fullREDMsinglelambda.RDS"))
  runs_diagREDMSL0 = readRDS(paste0(flder_in, generation, "_diagREDMsinglelambda.RDS"))
  runs_diagREDM0 = readRDS(paste0(flder_in, generation, "_diagREDM.RDS"))
  if(include_singleREDM){
    runs_singleREDM0 = readRDS(paste0(flder_in, generation, "_singleREDM.RDS"))
  }    
}

## match them all (wrt fullREM)
runs_fullREDMSL0 <- runs_fullREDMSL0[match(rownames(runs_fullREM0), rownames(runs_fullREDMSL0)),]
runs_diagREDMSL0 <- runs_diagREDMSL0[match(rownames(runs_fullREM0), rownames(runs_diagREDMSL0)),]
runs_diagREDM0 <- runs_diagREDM0[match(rownames(runs_fullREM0), rownames(runs_diagREDM0)),]
if(include_singleREDM){
  runs_singleREDM0 <- runs_singleREDM0[match(rownames(runs_fullREM0), rownames(runs_singleREDM0)),]
}

## Problem with convergence is acute in DM
table(is.na(runs_fullREM0$beta_est))
table(is.na(runs_fullREDMSL0$beta_est))
table(is.na(runs_diagREDMSL0$beta_est))
table(is.na(runs_diagREDM0$beta_est))

table((runs_fullREM0$converged))
table((runs_fullREDMSL0$converged))
table((runs_diagREDMSL0$converged))
table((runs_diagREDM0$converged))

barplot(c(runs_fullREM0=sum((runs_fullREM0$converged)),
             runs_fullREDMSL0=sum((runs_fullREDMSL0$converged)),
             runs_diagREDMSL0=sum((runs_diagREDMSL0$converged)),
             runs_diagREDM0=sum((runs_diagREDM0$converged))))
image(cbind(runs_fullREM0=((runs_fullREM0$converged)),
          runs_fullREDMSL0=((runs_fullREDMSL0$converged)),
          runs_diagREDMSL0=((runs_diagREDMSL0$converged)),
          runs_diagREDM0=((runs_diagREDM0$converged))))
# pheatmap(cor(cbind(runs_fullREM0=((runs_fullREM0$converged)),
#           runs_fullREDMSL0=((runs_fullREDMSL0$converged)),
#           runs_diagREDMSL0=((runs_diagREDMSL0$converged)),
#           runs_diagREDM0=((runs_diagREDM0$converged)))))

runs_fullREM <- runs_fullREM0[runs_fullREM0$converged,]
runs_fullREDMSL <- runs_fullREDMSL0[runs_fullREDMSL0$converged,]
runs_diagREDMSL <- runs_diagREDMSL0[runs_diagREDMSL0$converged,]
runs_diagREDM <- runs_diagREDM0[runs_diagREDM0$converged,]
if(include_singleREDM){
  runs_singleREDM <- runs_singleREDM0[runs_singleREDM0$converged,]
}

length(runs_fullREM)
length(runs_fullREDMSL)
length(runs_diagREDMSL)
length(runs_diagREDM)

####

system(paste0("mkdir -p ", flder_out, generation, "/summaries/"))

joint_df = cbind.data.frame(fullRE_M=runs_fullREM,
                            fullRE_DMSL=runs_fullREDMSL[match(rownames(runs_fullREM),
                                                                   rownames(runs_fullREDMSL)),],
                            diagRE_DMSL=runs_diagREDMSL[match(rownames(runs_fullREM),
                                                rownames(runs_diagREDMSL)),],
                            diagRE_DM=runs_diagREDM[match(rownames(runs_fullREM),
                                                rownames(runs_diagREDM)),])
if(include_singleREDM){
  joint_df <- cbind(joint_df,
                    singleRE_DM=runs_diagREDM[match(rownames(runs_singleREDM),
                                                  rownames(runs_singleREDM)),])
}
joint_df$names <- sapply(1:nrow(joint_df), function(i) gsub(paste0(joint_df[i,]$fullRE_M.idx_within_dataset, "$"), "", rownames(joint_df)[i]))
joint_df$names <- gsub("[.]$", "", joint_df$names)

sort(unique(gsub("\\..*","",rownames(joint_df))))


try({
pdf(paste0(flder_out, generation, "/summaries/betas_scatterplots.pdf"), height = 2.5)
do.call( 'grid.arrange', c(grobs=lapply(c('fullRE_M', 'fullRE_DMSL', 'diagRE_DMSL', 'diagRE_DM'), function(it_model){
  ggplot(joint_df, aes(x=fullRE_M.beta_true, y=get(paste0(it_model, '.beta_est'))))+geom_point()+theme_bw()+
    geom_abline(slope = 1, intercept = 0, lty='dashed', col='blue')+
    labs(x='True beta', y=paste0('Estimate from ', it_model))+
    annotate("text", label=paste0('rho= ', signif(cor(joint_df[,c('fullRE_M.beta_true')], joint_df[,paste0(it_model, '.beta_est')], use="complete.obs"),
                                                3)),
             x = Inf, y = -Inf, vjust=-0.4, hjust=1.0)
    # theme(text=element_text(family="LM Roman 10", size=20))
}), nrow=1))
dev.off()})


try({pdf(paste0(flder_out, generation, "/summaries/betas_scatterplots_colour.pdf"), height = 2.5)
do.call( 'grid.arrange', c(grobs=lapply(c('fullRE_M', 'fullRE_DMSL', 'diagRE_DMSL', 'diagRE_DM'), function(it_model){
  ggplot(joint_df, aes(x=fullRE_M.beta_true, y=get(paste0(it_model, '.beta_est')), col=fullRE_M.beta_gamma_shape))+geom_point()+theme_bw()+
    geom_abline(slope = 1, intercept = 0, lty='dashed', col='blue')+
    labs(x='True beta', y=paste0('Estimate from ', it_model))+
    annotate("text", label=paste0('rho= ', signif(cor(joint_df[,c('fullRE_M.beta_true')], joint_df[,paste0(it_model, '.beta_est')], use="complete.obs"),
                                                  3)),
             x = Inf, y = -Inf, vjust=-0.4, hjust=1.0)+theme(legend.position = "bottom")
  # theme(text=element_text(family="LM Roman 10", size=20))
}), nrow=1))
dev.off()
})

try({pdf(paste0(flder_out, generation, "/summaries/betas_scatterplots_colour_betagammashape.pdf"), height = 2.5, width=6)
  do.call( 'grid.arrange', c(grobs=lapply(c('fullRE_M', 'fullRE_DMSL', 'diagRE_DMSL', 'diagRE_DM'), function(it_model){
    ggplot(joint_df, aes(x=fullRE_M.beta_gamma_shape, y=get(paste0(it_model, '.beta_est')), col=fullRE_M.beta_gamma_shape))+geom_point()+theme_bw()+
      geom_abline(slope = 1, intercept = 0, lty='dashed', col='blue')+
      labs(x='Beta gamma shape', y=paste0('Estimate from ', it_model))+
      annotate("text", label=paste0('rho= ', signif(cor(joint_df[,c('fullRE_M.beta_gamma_shape')], joint_df[,paste0(it_model, '.beta_est')], use="complete.obs"),
                                                    3)),
               x = Inf, y = -Inf, vjust=-0.4, hjust=1.0)+theme(legend.position = "bottom")
    # theme(text=element_text(family="LM Roman 10", size=20))
  }), nrow=1))
  dev.off()
})

pdf(paste0(flder_out, generation, "/summaries/M_DM_comparison.pdf"), height = 3)
do.call('grid.arrange', list(ggplot(joint_df, aes(x=fullRE_M.beta_true, y=fullRE_M.beta_est, col=fullRE_M.d))+geom_point()+
                               geom_abline(intercept = 0, slope = 1)+theme_bw()+theme(legend.position = "bottom"),
                  ggplot(joint_df, aes(x=fullRE_DMSL.beta_true, y=fullRE_DMSL.beta_est, col=fullRE_M.d))+geom_point()+
                    geom_abline(intercept = 0, slope = 1)+theme_bw()+theme(legend.position = "bottom"), ncol=2))
dev.off()

pdf(paste0(flder_out, generation, "/summaries/M_DM_comparison_only_common.pdf"), height = 3)
do.call('grid.arrange', list(ggplot(joint_df[!is.na(joint_df$fullRE_M.beta_est) & !is.na(joint_df$fullRE_DMSL.beta_est),],
                                    aes(x=fullRE_M.beta_true, fullRE_M.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1)+
                               theme_bw()+theme(legend.position = "bottom"),
                             ggplot(joint_df[!is.na(joint_df$fullRE_M.beta_est) & !is.na(joint_df$fullRE_DMSL.beta_est),],
                                    aes(x=fullRE_DMSL.beta_true, fullRE_DMSL.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1)+
                               theme_bw()+theme(legend.position = "bottom"), ncol=2))
dev.off()


## Read in datasets
cat('Reading datasets\n')
datasets_files = list.files("../../../../data/assessing_models_simulation/datasets/", full.names = TRUE)
datasets_files = datasets_files[grep(pattern = paste0('/multiple_', generation, '_'), datasets_files)]
length(datasets_files)

# match. This used to be commented out up until the 20220130, where I have uncommented it for generation = "GenerationMixturePCAWG". I have modified it, though,
# to accomodate for multiple runs
# datasets_files = datasets_files[match(unique(sapply(rownames(joint_df), function(i) strsplit(i, "_dataset")[[1]][1])), ## what it used to be
#                                       gsub("_dataset.RDS", "", basename(datasets_files)))]
# stop('Problem with getting the same dataset twice. check analysis_aignle_cat_change.R to see how I did it for the updated version')
# datasets_files = datasets_files[match((sapply(rownames(joint_df[joint_df$fullRE_M.idx_within_dataset == 1,]), function(i) strsplit(i, "_dataset")[[1]][1])), ## new version
#                                       (gsub("_dataset[0-9]+.RDS", "", basename(datasets_files))))] ## commented out on 20220218
cat('Matching files\n')
datasets_files = datasets_files[match(unique(joint_df$names),
                                      gsub(".RDS", "", basename(datasets_files)))]
length(datasets_files)

datasets = lapply(datasets_files, readRDS)

if((generation %in% c("GenerationMixturePCAWG", "GenerationMixturefewersignaturesPCAWG", "GenerationMixturefewersignaturespairedPCAWG")) | grepl('GenerationMixturefewersignaturespaired', generation)  | grepl('signaturespaired', generation) ){
  cat('Transforming beta gamma shape from logR to probability')
  datasets <- lapply(datasets, function(i){
    if(i$beta_gamma_shape == -999){
      i$beta_gamma_shape = 0
    }else{
      i$beta_gamma_shape = softmax(c(i$beta_gamma_shape, 0))[1]
    }
    i
  })
}else{
  if(grepl('Mixture', generation)){
    stop('Are you sure you are using probabilities beta_gamma_shape for and not log-ratios?\n')
  }
}
# stop('Below: it shoud not be <unique()>')
# names(datasets) = unique(gsub("_dataset.RDS", "", basename(datasets_files))) ## commeted out on 20220218
names(datasets) = gsub("_dataset.RDS", "", basename(datasets_files))

DA_bool = ( sapply(datasets, function(i) i$beta_gamma_shape) > 0 )

runs_ttest_irl = lapply(datasets_files, function(i)  try(wrapper_run_ttest_ilr(i)))
hist(as.numeric(runs_ttest_irl), breaks=30); table(sapply(runs_ttest_irl, typeof))
runs_ttest_props = lapply(datasets_files, function(i)  try(wrapper_run_ttest_props(i)))
hist(as.numeric(runs_ttest_props), breaks=30); table(sapply(runs_ttest_props, typeof))
table(sapply(runs_ttest_irl, typeof), sapply(runs_ttest_props, typeof))
pvals_runs_HMP = lapply(datasets_files, function(i)  try(wrapper_run_HMP_Xdc.sevsample(i)))
pvals_runs_HMP2 = lapply(datasets_files, function(i)  try(wrapper_run_HMP_Xmcupo.sevsample(i)))
table(sapply(pvals_runs_HMP, typeof), sapply(pvals_runs_HMP2, typeof))
pvals_ttest_ilr = as.numeric(unlist(runs_ttest_irl))
pvals_ttest_ilr_adj = pvals_ttest_ilr
# pvals_perturbation1 =  lapply(datasets_files, function(i){.x <- readRDS(i); ## gives computationally singular systems
# aitchison_perturbation_test(.x$objects_counts, slot_name = "count_matrices_all")})
# pvals_perturbation2 =  lapply(datasets_files, function(i){.x <- readRDS(i); ## gives computationally singular systems
# aitchison_perturbation_test_alt(.x$objects_counts, slot_name = "count_matrices_all")})
pvals_perturbation =  unlist(sapply(lapply(datasets_files, function(i){.x <- readRDS(i);
aitchison_perturbation_test_alt_v2(.x$objects_counts, slot_name = "count_matrices_all")}), `[`, 'pval'))
pvals_permutation = unlist(lapply(datasets_files,  function(i){.x <- readRDS(i);
permutation_test_fun_wrapper(.x$objects_counts, nbootstraps=40)}))
Harris_chi <- sapply(datasets, function(i) iterative_chisqrt_wrapper(i$objects_counts) )

length(pvals_runs_HMP) == length(pvals_perturbation)
length(pvals_runs_HMP) == length(pvals_permutation)

# res_M = readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries/", generation, "_fullREM.RDS"))
# res_DM = readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries/", generation, "_fullREDM.RDS"))

## get the p-values for my models

## remove the last character because it determines what the row of the dataset is
names_datasets_uniq <- sort(unique(gsub(".RDS", "", names(datasets)))) #sort(unique(c(rownames(runs_fullREM0), rownames(runs_fullREDMSL0),
                        #             rownames(runs_diagREDMSL0), rownames(runs_diagREDM0))))
head(names_datasets_uniq)
names(datasets)

colnames(joint_df)

if(multiple_runs){
  pvals_list <- list()
  list_models <- c('fullRE_M', 'fullRE_DMSL', 'diagRE_DMSL', 'diagRE_DM')
  for(str_models in list_models){
    pvals_list[[str_models]] <- joint_df[,paste0(gsub('0$', '', str_models), '.pvals_adj')]
    # pvals_list[[str_models]][(rowSums(joint_df[,paste0(list_models, c('.converged'))]) == 4) %in% c(F, NA)] <- NA
    pvals_list[[str_models]][!(joint_df[,paste0(gsub('0$', '', str_models), '.converged')])] <- NA
    pvals_list[[str_models]] <- pvals_list[[str_models]][joint_df$fullRE_M.idx_within_dataset == 1]
  }
  ## commented out in 12 gen 2022
  # for(str_models in c('fullREM0', 'fullREDMSL0', 'diagREDMSL0', 'diagREDM0')){
  #   assign(gsub('0', '', paste0('pvals_', str_models)),
  #          get(paste0('runs_', str_models))[sapply(names_datasets_uniq,
  #          function(i) grep(i, rownames(get(paste0('runs_', str_models))))[1]),'pvals_adj'])
  #   assign(paste0('names(pvals_', gsub('0', '', str_models), ')'),
  #          gsub("_dataset.*", "", rownames(get(paste0('runs_', str_models)))[unique(get(paste0('runs_', str_models))$idx)]))
  #   
  #   get(paste0('pvals_', gsub('0', '', str_models)))
  #   
  #   ## remove p-vals of runs that didn't converge
  #   if(sum(!sapply(unique(get(paste0('runs_', str_models))$idx),
  #                  function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1])) > 0){
  #     # if there is any non-converged run
  #     cat('\nRemoving p-values of runs that did not converge in ', str_models, '\n')
  #     # get(paste0('pvals_', gsub('0', '', str_models))[!sapply(unique(get(paste0('runs_', str_models))$idx),
  #     #                                                           function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1])])
  #     print(typeof(get(paste0('pvals_', gsub('0', '', str_models)))))
  #     print(length(sapply(unique(get(paste0('runs_', str_models))$idx),
  #                         function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1])))
  #     print(table(sapply(unique(get(paste0('runs_', str_models))$idx),
  #                         function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1])))
  #     print(length(get(paste0('pvals_', gsub('0', '', str_models)))))
  #     ## this below gave problems in only diagDM in GenerationK2, and I don't know why. I am now using the alternative that follows it
  #     # assign(paste0('pvals_', gsub('0', '', str_models))[!sapply(unique(get(paste0('runs_', str_models))$idx),
  #     #                                                            function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1])], NA)
  #     cat('p vals which are going to be set to NA\n')
  #     print(!sapply(unique(get(paste0('runs_', str_models))$idx),
  #                   function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1]))
  #     cat('Current p-vals:')
  #     print(get(paste0('pvals_', gsub('0', '', str_models))))
  #     cat('Changing p vals\n')
  #     assign(paste0('get(pvals_', gsub('0', '', str_models), ')')[!sapply(unique(get(paste0('runs_', str_models))$idx),
  #                                                                         function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1])], NA)
  #     cat('Current p-vals:')
  #     print(get(paste0('pvals_', gsub('0', '', str_models))))
  #     
  #   }
  #   get(paste0('pvals_', gsub('0', '', str_models)))
  #   get(paste0('pvals_', gsub('0', '', (str_models))))[!sapply(unique(get(paste0('runs_', str_models))$idx),
  #                                                              function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1])]
  # }
}else{
  for(str_models in c('fullREM0', 'fullREDMSL0', 'diagREDMSL0', 'diagREDM0')){
    assign(gsub('0', '', paste0('pvals_', str_models)), get(paste0('runs_', str_models))[sapply(unique(get(paste0('runs_', str_models))$idx),
                           function(i) which(get(paste0('runs_', str_models))$idx == i)[1]),'pvals_adj'])
    assign(paste0('names(pvals_', gsub('0', '', str_models), ')'),
           gsub("_dataset.*", "", rownames(get(paste0('runs_', str_models)))[unique(get(paste0('runs_', str_models))$idx)]))
    
    ## remove p-vals of runs that didn't converge
    assign(paste0('pvals_', gsub('0', '', str_models))[!sapply(unique(get(paste0('runs_', str_models))$idx),
           function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1])], NA)
    get(paste0('pvals_', gsub('0', '', (str_models))))[!sapply(unique(get(paste0('runs_', str_models))$idx),
                                         function(i) get(paste0('runs_', str_models))[(get(paste0('runs_', str_models))$idx == i),'converged'][1])]
  }
}

pvals_fullREM <- pvals_list$fullRE_M
pvals_fullREDMSL <- pvals_list$fullRE_DMSL
pvals_diagREDMSL <- pvals_list$diagRE_DMSL
pvals_diagREDM <- pvals_list$diagRE_DM

all(names(pvals_fullREDMSL) == names(pvals_fullREM))
all(names(pvals_fullREDMSL) == names(pvals_diagREDM))

length(pvals_diagREDM)
length(pvals_fullREDMSL)

dim(runs_diagREDM0)
dim(runs_diagREDMSL0)
dim(runs_fullREDMSL0)
dim(runs_fullREM0)

pvals_fullREDMSL
pvals_diagREDM
pvals_diagREDMSL

sort(names(datasets))
sort(unique(gsub("\\..*","",rownames(joint_df))))

cat('Runs ', rownames(runs_fullREM0), '\n')
cat('Pvals ', names(pvals_fullREDMSL), '\n')
cat('Datasets\n ', paste0(sort(names(datasets)), '\n'))
cat('Number of runs ', length(pvals_fullREDMSL), '\n')
cat('Number of datasets ', length(datasets), '\n')
if(length(pvals_fullREDMSL) != length(datasets)){
  stop('The number of runs is not the number of datasets')
}

## p-values from my models are not adjusted for MT (below they are in adj)
pvals_data_frame=cbind.data.frame(pvals_fullREDMSL=pvals_fullREDMSL,
                                  pvals_fullREM=pvals_fullREM,
                                  pvals_diagREDMSL=pvals_diagREDMSL,
                                  pvals_diagREDM=pvals_diagREDM,
                                  ttest_props=unlist(runs_ttest_props),
                                  Harris_chi=Harris_chi,
                                  ttest_ilr_adj=pvals_ttest_ilr_adj,
                                  HMP=unlist(pvals_runs_HMP),
                                  HMP2=unlist(pvals_runs_HMP2),
                                  perturbation=pvals_perturbation,
                                  permutation=pvals_permutation,
                                  true=DA_bool)
head(pvals_data_frame)
print(pvals_fullREDMSL)
print(pvals_fullREM)
print(pvals_diagREDMSL)
print(pvals_diagREDM)
print(runs_diagREDM)

## select only runs that have converged for all models
DA_bool_all_converged <- DA_bool
DA_bool_all_converged[!(!is.na(pvals_fullREDMSL) & !is.na(pvals_fullREM) & !is.na(pvals_diagREDMSL) &
                          !is.na(pvals_diagREDM) & !is.na(Harris_chi) & !is.na(runs_ttest_props) &
                          !is.na(pvals_ttest_ilr_adj) & !is.na(pvals_runs_HMP) & !is.na(pvals_runs_HMP2) &
                          !is.na(pvals_perturbation) & !is.na(pvals_permutation))] <- NA
table(is.na(DA_bool_all_converged))
pvals_data_frame_all_converged = cbind.data.frame(pvals_fullREDMSL=pvals_fullREDMSL[!is.na(DA_bool_all_converged)],
                                                  pvals_fullREM=pvals_fullREM[!is.na(DA_bool_all_converged)],
                                                  pvals_diagREDMSL=pvals_diagREDMSL[!is.na(DA_bool_all_converged)],
                                                  pvals_diagREDM=pvals_diagREDM[!is.na(DA_bool_all_converged)],
                                                  ttest_props=unlist(runs_ttest_props)[!is.na(DA_bool_all_converged)],
                                                  Harris_chi=Harris_chi[!is.na(DA_bool_all_converged)],
                                                  ttest_ilr_adj=pvals_ttest_ilr_adj[!is.na(DA_bool_all_converged)],
                                                  HMP=unlist(pvals_runs_HMP)[!is.na(DA_bool_all_converged)],
                                                  HMP2=unlist(pvals_runs_HMP2)[!is.na(DA_bool_all_converged)],
                                                  perturbation=pvals_perturbation[!is.na(DA_bool_all_converged)],
                                                  permutation=pvals_permutation[!is.na(DA_bool_all_converged)],
                                                  true=DA_bool_all_converged[!is.na(DA_bool_all_converged)])
pvals_data_frame_adj <- adjust_all(pvals_data_frame, method='holm')
pvals_data_frame_all_converged_adj <- pvals_data_frame_adj[!is.na(DA_bool_all_converged),] ## remove after MT correction

cat('Creating <res_all> table')
res_all = give_res_all(pvals_data_frame)
xtable::xtable(res_all)
# xtable::xtable(res_all[,-ncol(res_all)])

res_all <- data.frame(res_all)
res_all$model = rownames(res_all)
res_all

## only results when we haev results for all the models
res_all_common_all_converged <-  give_res_all(pvals_data_frame_all_converged)
res_all_common_all_converged

ggplot(res_all, aes(x=1, y = FPR, col=model))+geom_point()

## group the runs by n, d, etc.
summarise_DA_detection(true = DA_bool, predicted = pvals_fullREDMSL < 0.05)

table(DA_bool, pvals_data_frame$ttest_props <= 0.05)

cat('Creating <put_vals_in_table> table\n')

dim(joint_df)
dim(pvals_data_frame)

match_pvals_params <- match(sapply(rownames(joint_df), function(i) strsplit(i, '_dataset')[[1]][1]),
                            sapply(rownames(pvals_data_frame), function(i) strsplit(i, '_dataset')[[1]][1]))
joint_df <- cbind(joint_df, pvals_data_frame[match_pvals_params,])

joint_df[joint_df$fullRE_M.beta_gamma_shape == 0,'pvals_diagREDM']
joint_df$pvals_diagREDM
ggplot(joint_df[joint_df$fullRE_M.beta_gamma_shape == 0,], aes(x=pvals_diagREDM))+geom_histogram()
ggplot(joint_df[joint_df$fullRE_M.beta_gamma_shape == 0,], aes(x=pvals_fullREDMSL))+geom_histogram()
ggplot(joint_df[joint_df$fullRE_M.beta_gamma_shape == 0,], aes(x=pvals_diagREDMSL))+geom_histogram()

cat('Creating tables as we vary params\n')

if(sum(!is.na(DA_bool_all_converged))>0){
  varying_d <-give_accuracies_with_varying_var('d')
  varying_d_all_converged <-give_accuracies_with_varying_var('d', datasets_arg = datasets[!is.na(DA_bool_all_converged)],
                                                             pvals_data_frame_arg = pvals_data_frame_all_converged)
  varying_n <-give_accuracies_with_varying_var('n')
  varying_n_all_converged <-give_accuracies_with_varying_var('n', datasets_arg = datasets[!is.na(DA_bool_all_converged)],
                                               pvals_data_frame_arg = pvals_data_frame_all_converged)
  varying_n_adj <-give_accuracies_with_varying_var('n', pvals_data_frame_arg = pvals_data_frame_adj)
  varying_n_all_converged_adj <-give_accuracies_with_varying_var('n', datasets_arg = datasets[!is.na(DA_bool_all_converged)],
                                                             pvals_data_frame_arg = pvals_data_frame_all_converged_adj)
  varying_betashape <-give_accuracies_with_varying_var('beta_gamma_shape')
  varying_betashape_all_converged <-give_accuracies_with_varying_var('beta_gamma_shape', datasets_arg = datasets[!is.na(DA_bool_all_converged)],
                                                                     pvals_data_frame_arg = pvals_data_frame_all_converged)
  varying_betashape_adj <-give_accuracies_with_varying_var('beta_gamma_shape', pvals_data_frame_arg = pvals_data_frame_adj)
  varying_betashape_all_converged_adj <-give_accuracies_with_varying_var('beta_gamma_shape', datasets_arg = datasets[!is.na(DA_bool_all_converged)],
                                                                     pvals_data_frame_arg = pvals_data_frame_all_converged_adj)

  try({varying_n_d <-give_accuracies_with_varying_var(var = c('d', 'n'), two_var = T)})
  try({varying_n_betashape <-give_accuracies_with_varying_var(var = c('n', 'beta_gamma_shape'), two_var = T)})
  try({varying_d_betashape <-give_accuracies_with_varying_var(var = c('d', 'beta_gamma_shape'), two_var = T)})
  
  cat('Creating tables as we vary params\n')
  
  ggplot(varying_d, aes(x=d, y = FPR, col=model, group=model))+geom_point()+geom_line()+theme_bw()#+facet_wrap(.~mod, noel)
  ggplot(varying_n, aes(x=n, y = FPR, col=model, group=model))+geom_point()+geom_line()+theme_bw()#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/FPR_with_n.pdf"),
         height = 3.0, width = 4)
  
  cat('Plotting: varying n\n')
  table(DA_bool_all_converged)
  print(datasets[!is.na(DA_bool_all_converged)])
  print(varying_n_all_converged)
  
  ggplot(varying_n_all_converged, aes(x=n, y = FPR, col=model, group=model))+geom_point()+geom_line()+theme_bw()+
    labs(col=FALSE)#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/FPR_with_n_all_converged.pdf"),
         height = 3.0, width = 4)
  
  cat('Plotting: varying betashape\n')
  
  ggplot(varying_betashape, aes(x=beta_gamma_shape, y = FPR, col=model, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model)
  
  cat('Plotting: varying d\n')
  ggplot(varying_d, aes(x=d, y = AUC, col=model, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/AUC_with_d.pdf"),
         height = 3.0, width = 4.5)
  ggplot(varying_betashape, aes(x=beta_gamma_shape, y = AUC, col=model, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/AUC_with_betashape.pdf"),
         height = 3.0, width = 4.5)
  
  cat('Plotting: varying d (2)\n')
  ggplot(varying_d_all_converged, aes(x=d, y = AUC, col=model, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/AUC_with_d_all_converged.pdf"),
         height = 3.0, width = 4.5)
  
  ggplot(varying_n, aes(x=n, y = AUC, col=model, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/AUC_with_n.pdf"),
         height = 3.0, width = 4.5)
  ggplot(varying_n_all_converged, aes(x=n, y = AUC, col=model, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/AUC_with_n_all_converged.pdf"),
         height = 3.0, width = 5.5)
  
  ggplot(varying_betashape, aes(x=beta_gamma_shape, y = AUC, col=model, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model)
  
  ggplot(varying_d, aes(x=d, y = Accuracy, col=model, group=model))+geom_point()+geom_line()+theme_bw()#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_d.pdf"),
         height = 3.0, width = 4.0)
  ggplot(varying_d_all_converged, aes(x=d, y = Accuracy, col=model, group=model))+geom_point()+geom_line()+theme_bw()#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_d_all_converged.pdf"),
         height = 3.0, width = 4.0)
  
  cat('Line 532\n')
  
  ggplot(varying_n, aes(x=n, y = Accuracy, col=model, group=model))+geom_point()+geom_line()+theme_bw()#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_N.pdf"),
         height = 3.0, width = 4.0)
  ggplot(varying_n_all_converged, aes(x=n, y = Accuracy, col=model, group=model))+geom_point()+geom_line()+theme_bw()#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_N_all_converged.pdf"),
         height = 3.0, width = 4.0)
  
  ggplot(varying_betashape, aes(x=beta_gamma_shape+.001, y = Accuracy, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model, nrow=2)+scale_x_continuous(trans = "log10")
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_models.pdf"),
         height = 3.0, width = 6.0)
  
  print(varying_betashape_all_converged$beta_gamma_shape)
  ggplot(varying_betashape_all_converged, aes(x=beta_gamma_shape+.001, y = Accuracy, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model, nrow=2)+scale_x_continuous(trans = "log10")
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_models_all_converged.pdf"),
         height = 3.0, width = 6.0)
  
  ggplot(varying_betashape, aes(x=beta_gamma_shape+.001, y = FPR, group=model))+geom_point()+geom_line()+theme_bw()+facet_wrap(.~model, nrow=2)+scale_x_continuous(trans = "log10")
  # ggsave(paste0(flder_out, generation, "/summaries/accuracy_models.pdf"),
  #        height = 3.5, width = 8)
  ## add fraction of correct classification, which are values I should have for all combinations
  
  print(varying_n)
  
  ggplot(varying_n, aes(x=n, y = WeightedAccuracy, col=model, group=model))+geom_point()+geom_line()+theme_bw()#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_N.pdf"),
         height = 3.0, width = 4.0)
  print(varying_n_all_converged)
  ggplot(varying_n_all_converged, aes(x=n, y = WeightedAccuracy, col=model, group=model))+geom_point()+geom_line()+theme_bw()#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_N_all_converged.pdf"),
         height = 3.0, width = 4.0)
  
  cat('Plotting varying n\n')
  ggplot(varying_n_all_converged, aes(x=n, y = WeightedAccuracy, col=model, group=model))+geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+guides(col=FALSE)+ggtitle(generation)#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_N_all_converged_palette2.pdf"),
         height = 3.0, width = 4.0)
  max_n <- max(varying_n_all_converged$n)
  min_n <- min(varying_n_all_converged$n)
  ggplot(varying_n_all_converged, aes(x=n, y = WeightedAccuracy, col=model, group=model,
                                      lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),
                                      label=model))+
    geom_point()+ geom_line()+theme_bw()+
    geom_label_repel(data = varying_n_all_converged[varying_n_all_converged$n == max_n,],
                     max.overlaps = Inf, aes(x=max(n)), direction = "y", nudge_x=max_n*0.3,force=100,
                     size=3)+
    # geom_text(data = varying_n_all_converged[varying_n_all_converged$n == max_n,],
    #                  aes(x=n*1.2, y=WeightedAccuracy),
    #                  size=3)+
    lims(x=c(min_n, max_n*1.3))+
                     # )+
    scale_color_manual(values=colours_models)+guides(col=FALSE, lty='none')+
    ggtitle(gsub("Generation", "", gsub("generation", "", generation)))#+facet_wrap(.~model)
    # ggtitle(generation)#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_N_all_converged_palette2_names.pdf"),
         height = 3.0, width = 3.5)
  ggplot(varying_n_all_converged, aes(x=n, y = WeightedAccuracy, col=model, group=model,
                                      lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),
                                      label=model))+
    geom_point()+ geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+guides(col=FALSE, lty='none')+
    ggtitle(gsub("Generation", "", gsub("generation", "", generation)))#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_N_all_converged_palette2_lty.pdf"),
         height = 3.0, width = 3.0)

  ggplot(varying_n_all_converged_adj, aes(x=n, y = WeightedAccuracy, col=model, group=model,
                                      lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),
                                      label=model))+
    geom_point()+ geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+guides(col=FALSE, lty='none')+
    ggtitle(gsub("Generation", "", gsub("generation", "", generation)))#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_adj_with_N_all_converged_palette2_lty.pdf"),
         height = 3.0, width = 3.0)

    
  ggplot(varying_n, aes(x=n, y = WeightedAccuracy, col=model, group=model))+geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+labs(col='')+guides(col=FALSE)+ggtitle(generation)#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_N_palette2.pdf"),
         height = 3.0, width = 4.0)
  ggplot(varying_n, aes(x=n, y = WeightedAccuracy, col=model, group=model,
                                      lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),
                                      label=model))+
    geom_point()+ geom_line()+theme_bw()+
    geom_label_repel(data = varying_n[varying_n$n == max_n,],
                     max.overlaps = Inf, aes(x=max(n)), direction = "y", nudge_x=max_n*0.3,force=100,
                     size=3)+
    lims(x=c(min_n, max_n*1.3))+
    # )+
    scale_color_manual(values=colours_models)+guides(col=FALSE, lty='none')+
    ggtitle(gsub("Generation", "", gsub("generation", "", generation)))#+facet_wrap(.~model)
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_N_palette2_names.pdf"),
         height = 3.0, width = 3.5)
  
  ggplot(varying_betashape, aes(x=beta_gamma_shape, y = Accuracy, col=model, group=model, label=model,
                                lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+labs(col='')+guides(col='none', lty='none')+
    ggtitle(generation)+#+facet_wrap(.~model)
    geom_label_repel(data = varying_betashape[varying_betashape$beta_gamma_shape == max(varying_betashape$beta_gamma_shape),],
                   max.overlaps = Inf, aes(x=max(varying_betashape$beta_gamma_shape)), direction = "y", nudge_x=max_n*0.3,force=100,
                   size=3)+
    lims(x=c(min(varying_betashape$beta_gamma_shape), max(varying_betashape$beta_gamma_shape)*1.5))
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_palette2.pdf"),
         height = 3.0, width = 4.0)
  ggplot(varying_betashape_adj, aes(x=beta_gamma_shape, y = Accuracy, col=model, group=model, label=model,
                                lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+labs(col='')+guides(col='none', lty='none')+
    ggtitle(generation)+#+facet_wrap(.~model)
    geom_label_repel(data = varying_betashape_adj[varying_betashape_adj$beta_gamma_shape == max(varying_betashape_adj$beta_gamma_shape),],
                     max.overlaps = Inf, aes(x=max(varying_betashape_adj$beta_gamma_shape)), direction = "y", nudge_x=max_n*0.3,force=100,
                     size=3)+
    lims(x=c(min(varying_betashape_adj$beta_gamma_shape), max(varying_betashape_adj$beta_gamma_shape)*1.5))
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_adj_with_betagammashape_palette2.pdf"),
         height = 3.0, width = 4.0)
  
  ggplot(varying_betashape_all_converged, aes(x=beta_gamma_shape, y = Accuracy, col=model, group=model, label=model,
                                lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+labs(col='')+guides(col='none', lty='none')+
    ggtitle(generation)+#+facet_wrap(.~model)
    geom_label_repel(data = varying_betashape_all_converged[varying_betashape_all_converged$beta_gamma_shape == max(varying_betashape_all_converged$beta_gamma_shape),],
                     max.overlaps = Inf, aes(x=max(varying_betashape_all_converged$beta_gamma_shape)), direction = "y",
                     nudge_x= max(varying_betashape_all_converged$beta_gamma_shape)*0.3,force=100,
                     size=3)+
    lims(x=c(min(varying_betashape_all_converged$beta_gamma_shape), max(varying_betashape_all_converged$beta_gamma_shape)*1.5))
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_all_converged_palette2.pdf"),
         height = 3.0, width = 4.0)
  
  ggplot(varying_betashape_all_converged_adj, aes(x=beta_gamma_shape, y = Accuracy, col=model, group=model, label=model,
                                              lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+labs(col='')+guides(col='none', lty='none')+
    ggtitle(generation)+#+facet_wrap(.~model)
    geom_label_repel(data = varying_betashape_all_converged_adj[varying_betashape_all_converged_adj$beta_gamma_shape == max(varying_betashape_all_converged$beta_gamma_shape),],
                     max.overlaps = Inf, aes(x=max(varying_betashape_all_converged_adj$beta_gamma_shape)), direction = "y",
                     nudge_x= max(varying_betashape_all_converged_adj$beta_gamma_shape)*0.3,force=100,
                     size=3)+
    lims(x=c(min(varying_betashape_all_converged_adj$beta_gamma_shape), max(varying_betashape_all_converged_adj$beta_gamma_shape)*1.5))
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_adj_with_betagammashape_all_converged_palette2.pdf"),
         height = 3.0, width = 4.0)
  
  if((generation %in% c("GenerationMixturePCAWG", "GenerationMixturefewersignaturesPCAWG", "GenerationMixturefewersignaturespairedPCAWG")) | grepl('GenerationMixturefewersignaturespaired', generation) | grepl('signaturespaired', generation)){
    varying_betashape$beta_gamma_shape <- signif(varying_betashape$beta_gamma_shape, 2)
  }
  
  if(grepl('Mixture', generation)){
    title_x <- 'Percentage of mixture'
  }else{
    title_x <- 'Gamma'
  }
  ggplot(varying_betashape, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
                                lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+labs(col='', x=title_x)+guides(col='none', lty='none')+
    ggtitle(generation)+#+facet_wrap(.~model)
    geom_label_repel(data = varying_betashape[varying_betashape$beta_gamma_shape == max(varying_betashape$beta_gamma_shape),],
                     max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
                     nudge_x=max(varying_betashape$beta_gamma_shape)+1,
                     size=3)+
    coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.3))+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
    # lims(x=c(min(varying_betashape$beta_gamma_shape), max(varying_betashape$beta_gamma_shape)*1.5))
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_palette2_factor.pdf"),
         height = 3.0, width = 4.0)
  
  ggplot(varying_betashape_adj, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
                                lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+labs(col='', x=title_x)+guides(col='none', lty='none')+
    ggtitle(generation)+#+facet_wrap(.~model)
    geom_label_repel(data = varying_betashape_adj[varying_betashape_adj$beta_gamma_shape == max(varying_betashape_adj$beta_gamma_shape),],
                     max.overlaps = Inf, aes(x=factor(max(varying_betashape_adj$beta_gamma_shape))), direction = "y",
                     nudge_x=max(varying_betashape_adj$beta_gamma_shape)+1,
                     size=3)+
    coord_cartesian(xlim = c(0, length(unique(varying_betashape_adj$beta_gamma_shape))*1.3))+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
  # lims(x=c(min(varying_betashape_adj$beta_gamma_shape), max(varying_betashape_adj$beta_gamma_shape)*1.5))
  ggsave(paste0(flder_out, generation, "/summaries/accuracy_adj_with_betagammashape_palette2_factor.pdf"),
         height = 3.0, width = 4.0)
  
  
  # ggplot(varying_n_betashape, aes(x=beta_gamma_shape+.001, y = Accuracy, group=model, col=n))+
  #   geom_point()+geom_line()+theme_bw()+facet_wrap(.~model, nrow=2)+scale_x_continuous(trans = "log10")
  # ggplot(varying_n_betashape, aes(x=n, y = Accuracy, group=model, col=beta_gamma_shape))+
  #   geom_point()+geom_line()+theme_bw()+facet_wrap(.~model, nrow=2)+scale_x_continuous(trans = "log10")
  # ggplot(varying_d_betashape, aes(x=beta_gamma_shape+.001, y = Accuracy, group=model, col=d))+
    # geom_point()+geom_line()+theme_bw()+facet_wrap(.~model, nrow=2)+scale_x_continuous(trans = "log10")
  
  varying_n[varying_n$n == 100,]
  
  table(true=DA_bool, diagREDMSL=pvals_diagREDMSL < 0.05)
  table(true=DA_bool, HMP=pvals_runs_HMP < 0.05)
  
  ## plot only the non-DA
  pvals_data_frame[pvals_data_frame$true == F,]
  cat('Plotting varying_betashape\n')
  ggplot(varying_betashape[varying_betashape$beta_gamma_shape == 0,],
         aes(y=FPR, x=factor(model, levels=c("fullREM", "fullREDMSL", "diagREDMSL", "diagREDM", "HMP", "HMP2",
                                             "Harris_chi",  "permutation", "perturbation", "ttest", "ILR"))))+
           geom_point()+labs(x="")+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
  ggsave(paste0(flder_out, generation, "/summaries/FDR_nonDA.pdf"),
         height = 3.0, width = 4.0)
  
  try({
  ggplot(varying_n_betashape[varying_n_betashape$beta_gamma_shape == 0,],
         aes(y=FPR, x=factor(model, levels=c("fullREM", "fullREDMSL", "diagREDMSL", "diagREDM", "HMP", "HMP2",
                                             "Harris_chi",  "permutation", "perturbation", "ttest", "ILR")),
             col=n))+geom_violin()+
    geom_point()+labs(x="")+
    theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
  ggsave(paste0(flder_out, generation, "/summaries/FDR_nonDA_varyingn.pdf"),
         height = 3.0, width = 4.0)
  })
  
  ##-----------------------------------------------------------------##
  
  ## now only including datasets for which we have results for all models
  res_all_common_all_converged
  
  ##-----------------------------------------------------------------##
  ggplot(varying_n[grepl('diag', varying_n$model),],
         aes(y=FPR, x=n, col=model))+geom_point()+geom_line()
  
  varying_n[grepl('diag', varying_n$model),'FPR']
  
  ##' the interesting part is the second list, first column: statistically signif values,
  ##' but with no differential abundance (beta_gamma_shape = 0)
  ##' There is one false positive when n=10  and when n=100. Otherwise there are zero. Then, shouldn't
  ##' the FPR change, and be zero when n is 20 and 50?
  table(unlist(sapply(datasets, `[`, 'n')),
  unlist(sapply(datasets, `[`, 'beta_gamma_shape')),
  pvals_diagREDMSL < 0.05)
  ##-----------------------------------------------------------------##
  
  table(DA_bool, M_est=pvals_fullREM <= 0.05)
  table(DA_bool, DM_est=pvals_fullREDMSL <= 0.05)
  table(DA_bool, ILR_est=pvals_ttest_ilr_adj <= 0.05)
  
  xtable::xtable(table(DA_bool, M_est=pvals_fullREM <= 0.05))
  xtable::xtable(table(DA_bool, DM_est=pvals_fullREDMSL <= 0.05))
  xtable::xtable(table(DA_bool, ILR_est=pvals_ttest_ilr_adj <= 0.05))
  
  head(melt(list(table(DA_bool, M_est=pvals_fullREM <= 0.05),
                 table(DA_bool, DM_est=pvals_fullREDMSL <= 0.05),
                 table(DA_bool, ILR_est=pvals_ttest_ilr_adj <= 0.05))
  ))
  
  
  joint_df$fullRE_M.beta_true
  
  colnames(joint_df)[grepl('onverged', colnames(joint_df))]
  joint_df_converged <- (melt(joint_df, measure.vars = c("fullRE_M.converged", "fullRE_DMSL.converged",
                                       "diagRE_DMSL.converged", "diagRE_DM.converged")))
  ggplot(joint_df_converged, aes(x=diagRE_DM.d, fill=value))+geom_bar()+facet_wrap(.~gsub(".converged", "", variable), nrow=1)+theme_bw()+
    theme(legend.position = "bottom")+labs(x="d")
  ggsave(paste0(flder_out, generation, "/summaries/convergence-varying_d.pdf"),
         height = 3.0, width = 6.0)
  
  cat('Plotting joint_df_converged\n')
  ggplot(joint_df_converged, aes(x=diagRE_DM.n, fill=value))+geom_bar()+facet_wrap(.~gsub(".converged", "", variable), nrow=1)+theme_bw()+
    theme(legend.position = "bottom")+labs(x="n")
  ggsave(paste0(flder_out, generation, "/summaries/convergence-varying_n.pdf"),
         height = 3.0, width = 6.0)
  
  joint_df_converged$value2 <- joint_df_converged$value
  joint_df_converged$value2[is.na(joint_df_converged$value2 )] <- FALSE
  
  joint_df_converged$n = apply(joint_df_converged[,c('diagRE_DM.n', 'diagRE_DMSL.n', 'fullRE_DMSL.n', 'fullRE_M.n')],
                               1, function(i){.x <- unique(i); .x[!is.na(.x)]})
  joint_df_converged$d = apply(joint_df_converged[,c('diagRE_DM.d', 'diagRE_DMSL.d', 'fullRE_DMSL.d', 'fullRE_M.d')],
                               1, function(i){.x <- unique(i); .x[!is.na(.x)]})
  
  cat('Plotting joint_df_converged by summary\n')
  ggplot(joint_df_converged %>% group_by(joint_df_converged$d,
                                         joint_df_converged$n, variable) %>%
           dplyr::summarise(mean=mean(value2)), aes(x=factor(`joint_df_converged$n`), y=`joint_df_converged$d`, fill=mean))+
    geom_tile()+
    facet_wrap(.~gsub(".converged", "", variable), nrow=1)+
    theme_bw()+
    theme(legend.position = "bottom")+labs(x="n")+
    jcolors::scale_fill_jcolors_contin("pal3", reverse = FALSE)+labs(y='d', fill='Fraction of converged runs')
  ggsave(paste0(flder_out, generation, "/summaries/convergence-varying_n_d.pdf"),
         height = 2.5, width = 7.0)
  
  # ggplot(joint_df,
  #        aes(x=fullRE_M.idx_within_dataset, y=(fullRE_M.beta_est), col=fullRE_M.pvals_adj<0.05))+
  #   geom_point(aes(x=fullRE_M.idx_within_dataset, y=fullRE_M.beta_true), shape=4)+
  #   geom_abline(slope = 0, intercept = 0, alpha=0.2)+geom_point(aes(shape=fullRE_M.DA_bool))+
  #   geom_errorbar(aes(ymin=fullRE_M.beta_est-1.96*fullRE_M.beta_stderr, ymax=fullRE_M.beta_est+1.96*fullRE_M.beta_stderr), width=.2,
  #                 position=position_dodge(.9))+
  #   facet_wrap(.~interaction(fullRE_M.idx, fullRE_M.beta_gamma_shape, fullRE_M.DA_bool), scales='free_x', nrow=length(unique(joint_df$fullRE_M.beta_gamma_shape)))+
  #   theme_bw()+theme(legend.position = "bottom")
  #   # scale_colour_viridis_d(option = "plasma")
  #   # scale_colour_manual(values = c("red","#2e8b57", "red", "#2e8b57")) #"#2e8b57"))
  # ggsave(paste0(flder_out, generation, "/summaries/M_betaslopes_confint.pdf"),
  #        height = 14, width = 8)
  
  # ggplot(joint_df,
  #        aes(x=DM.idx_within_dataset, y=(DM.beta_est), col=DM.pvals_adj<0.05))+
  #   geom_point(aes(x=M.idx_within_dataset, y=M.beta_true), shape=4)+
  #   geom_abline(slope = 0, intercept = 0, alpha=0.2)+geom_point(aes(shape=DM.DA_bool))+
  #   geom_errorbar(aes(ymin=DM.beta_est-1.96*DM.beta_stderr, ymax=DM.beta_est+1.96*DM.beta_stderr), width=.2,
  #                 position=position_dodge(.9))+
  #   facet_wrap(.~interaction(DM.idx, DfullRE_M.beta_gamma_shape, DM.DA_bool), scales='free_x', nrow=length(unique(joint_df$DfullRE_M.beta_gamma_shape)))+
  #   theme_bw()+theme(legend.position = "bottom")
  # ggsave(paste0(flder_out, generation, "/summaries/DM_betaslopes_confint.pdf"),
         # height = 14, width = 8)
  
  all(joint_df$fullRE_M.DA_bool == joint_df$DM.DA_bool)
  table(truth=joint_df$fullRE_M.DA_bool, M=joint_df$fullRE_M.pvals_adj <= 0.05)
  table(truth=joint_df$DM.DA_bool, DM=joint_df$DM.pvals_adj <= 0.05)
  
  ## ROC curves
  
  colnames(joint_df)
  
  all(joint_df$fullRE_M.d == joint_df$fullRE_DMSL.d, na.rm = T); all(joint_df$fullRE_M.n == joint_df$fullRE_DMSL.n, na.rm = T); all(joint_df$fullRE_M.beta_gamma_shape == joint_df$fullRE_DMSL.beta_gamma_shape, na.rm = T);  all(joint_df$fullRE_M.DA_bool == joint_df$fullRE_DMSL.DA_bool, na.rm = T); 
  # "fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape", M.pvals_adj", "fullRE_M.DA_bool", "fullREDMSL.pvals_adj" 
  ## select datasets, not betas
  
  cat('Plotting grouping 1\n')
  
  joint_df$M_type = paste0(joint_df$fullRE_M.DA_bool, joint_df$fullRE_M.pvals_adj < 0.05)
  joint_df$DM_type = paste0(joint_df$fullRE_DMSL.DA_bool, joint_df$fullRE_DMSL.pvals_adj < 0.05)
  joint_df_grouping_by_n = joint_df[sapply(unique(joint_df$fullRE_DMSL.idx), function(i) which(joint_df$fullRE_DMSL.idx == i)[1]),] %>%
    dplyr::select(c("fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape", "fullRE_M.pvals_adj", "fullRE_M.DA_bool",
             "fullRE_DMSL.pvals_adj", "M_type", "DM_type" )) %>% 
    group_by(fullRE_M.beta_gamma_shape, fullRE_M.n) %>%
    mutate(sensitivity_DM=sum(DM_type == 'TRUETRUE')/sum(DM_type %in% c('TRUETRUE', 'TRUEFALSE')),
           sensitivity_M=sum(M_type == 'TRUETRUE')/sum(M_type %in% c('TRUETRUE', 'TRUEFALSE')),
           specificity_DM=sum(DM_type == 'FALSEFALSE')/sum(DM_type %in% c('FALSETRUE', 'FALSEFALSE')),
           specificity_M=sum(M_type == 'FALSEFALSE')/sum(M_type %in% c('FALSETRUE', 'FALSEFALSE')))
  
  cat('Plotting grouping 2\n')
  joint_df_grouping_by_d = joint_df[sapply(unique(joint_df$fullRE_DMSL.idx), function(i) which(joint_df$fullRE_DMSL.idx == i)[1]),] %>%
    dplyr::select(c("fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape", "fullRE_M.pvals_adj",
                    "fullRE_M.DA_bool", "fullRE_DMSL.pvals_adj", "M_type", "DM_type" )) %>% 
    group_by(fullRE_M.beta_gamma_shape, fullRE_M.d) %>%
    mutate(sensitivity_DM=sum(DM_type == 'TRUETRUE')/sum(DM_type %in% c('TRUETRUE', 'TRUEFALSE')),
           sensitivity_M=sum(M_type == 'TRUETRUE')/sum(M_type %in% c('TRUETRUE', 'TRUEFALSE')),
           specificity_DM=sum(DM_type == 'FALSEFALSE')/sum(DM_type %in% c('FALSETRUE', 'FALSEFALSE')),
           specificity_M=sum(M_type == 'FALSEFALSE')/sum(M_type %in% c('FALSETRUE', 'FALSEFALSE')))
  
  
  ggplot(droplevels(joint_df_grouping_by_n), aes(x=sensitivity_M, y=1-specificity_M))+
    geom_point()
  ggplot(droplevels(joint_df_grouping_by_d), aes(x=sensitivity_M, y=1-specificity_M))+
    geom_point()
  
  ggplot(droplevels(joint_df_grouping_by_n), aes(x=sensitivity_DM, y=1-specificity_DM,
                                                 col=interaction(fullRE_M.d, fullRE_M.n)))+
    geom_point()
  
  cat('Plotting joint_df_grouping_by_n 1\n')
  
  ggplot(droplevels(joint_df_grouping_by_n), aes(x=fullRE_M.n, y=sensitivity_DM,
                                                 group=fullRE_M.beta_gamma_shape))+
    geom_point()+geom_line()+facet_wrap(.~fullRE_M.beta_gamma_shape)
  ggsave(paste0(flder_out, generation, "/summaries/fullRE_M.beta_gamma_shape_sensitivity.pdf"))
  
  # head(melt(joint_df_grouping_by_n, id.vars = c("fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape", "fullRE_M.pvals_adj", "M_type", "fullRE_M.DA_bool", "sensitivity_M", "specificity_DM")))
  # plot(joint_df_grouping_by_n$sensitivity_M, 1-joint_df_grouping_by_n$specificity_M)
  # 
  # head(melt(joint_df_grouping_by_n[,c("fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape", "fullRE_M.pvals_adj", "M_type", "fullRE_M.DA_bool", "sensitivity_M", "specificity_DM")],
  #           id.vars = c("fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape")))
  
  modify_colnames = function(i){
    colnames(i) = gsub("fullREDMSL.", "", colnames(i))
    colnames(i) = gsub("fullRE_M.", "", colnames(i))
    colnames(i) = gsub("_M", "", colnames(i))
    colnames(i) = gsub("_DM", "", colnames(i))
    i
  }
  
  cat('Plotting joint_df_grouping_by_n_alt\n')
  joint_df_grouping_by_n_alt = rbind(modify_colnames(cbind(joint_df_grouping_by_n[,c("fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape", "fullRE_M.pvals_adj", "M_type", "fullRE_M.DA_bool", "sensitivity_M", "specificity_M")], model="M")),
        modify_colnames(cbind(joint_df_grouping_by_n[,c("fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape", "fullRE_M.pvals_adj", "DM_type", "fullRE_M.DA_bool", "sensitivity_DM", "specificity_DM")], model="DM")))
  
  joint_df_grouping_by_d_alt = rbind(modify_colnames(cbind(joint_df_grouping_by_d[,c("fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape", "fullRE_M.pvals_adj", "M_type", "fullRE_M.DA_bool", "sensitivity_M", "specificity_M")], model="M")),
                                     modify_colnames(cbind(joint_df_grouping_by_d[,c("fullRE_M.d", "fullRE_M.n", "fullRE_M.beta_gamma_shape", "fullRE_M.pvals_adj", "DM_type", "fullRE_M.DA_bool", "sensitivity_DM", "specificity_DM")], model="DM")))
  
  try({
  ggplot(droplevels(joint_df_grouping_by_n_alt), aes(x=d, y=sensitivity, group=interaction(n, model), col=n, shape=model))+
    geom_point()+geom_line()+facet_wrap(.~beta_gamma_shape)
  
  ggplot(droplevels(joint_df_grouping_by_d_alt), aes(x=d, y=sensitivity, group=interaction(n, model), col=n, shape=model))+
    geom_point()+geom_line()+facet_wrap(.~beta_gamma_shape)
  ggplot(droplevels(joint_df_grouping_by_d_alt), aes(x=beta_gamma_shape, y=sensitivity, group=interaction(d, model), col=n, shape=model))+
    geom_point()+geom_line()
  
  ggplot(droplevels(joint_df_grouping_by_d_alt), aes(x=beta_gamma_shape, y=sensitivity, group=interaction(beta_gamma_shape, model), col=model))+
    geom_point()+geom_line()+geom_boxplot()+scale_x_continuous(trans = "log2")
  
  ggplot(droplevels(joint_df_grouping_by_d_alt), aes(x=beta_gamma_shape, y=specificity, group=interaction(beta_gamma_shape, model), col=model))+
    geom_point()+geom_line()+geom_boxplot()+scale_x_continuous(trans = "log2")
  })
  
  ## What hasn't run?
  try({
    cat(paste0(sapply(rownames(joint_df[is.na(joint_df$DM.beta_est),]), function(i) substr(i, start = 1, stop = nchar(i)-1)) %>% unique, sep='" "', collapse=''))
  
  sapply(rownames(joint_df[is.na(joint_df$DM.beta_est),]), function(i) gsub("_dataset", "", substr(i, start = 1, stop = nchar(i)-1))) %>% unique
  
  ggplot(joint_df, aes(x=fullRE_M.d, y=as.numeric(fullRE_M.converged), group=interaction(fullRE_M.n, fullRE_M.d)))+
    geom_jitter(height = 0.1, alpha=0.2, col='#07367d')+geom_violin()+facet_wrap(.~fullRE_M.n)+
    theme_bw()+labs(x='Number of categories (d)', y='Number of successful (1) or unsuccessful (0) convergences')+
    ggtitle('Success in convergence of Multinomial runs')
  ggsave(paste0(flder_out, generation, "/summaries/M_d_n_convergence.pdf"))
  
  ggplot(joint_df, aes(x=fullRE_DMSL.d, y=as.numeric(fullRE_DMSL.converged),
                       group=interaction(fullRE_DMSL.n, fullRE_DMSL.d)))+
    geom_jitter(height = 0.1, alpha=0.2, col='#07367d')+geom_violin()+facet_wrap(.~fullRE_DMSL.n)+
    theme_bw()+labs(x='Number of categories (d)', y='Number of successful (1) or unsuccessful (0) convergences')+
    ggtitle('Success in convergence of Dirichlet-Multinomial runs')
  ggsave(paste0(flder_out, generation, "/summaries/DM_d_n_convergence.pdf"))
  })
  
  # ggplot(joint_df, aes(x=M.d, y=as.numeric(M.converged), group=interaction(M.n, M.d), col=fullRE_M.beta_gamma_shape))+
  #   geom_jitter(height = 0.1, alpha=0.8)+geom_violin()+facet_wrap(.~M.n)+
  #   theme_bw()+labs(x='Number of categories (d)', y='Number of successful (1) or unsuccessful (0) convergences')+
  #   ggtitle('Success in convergence of Multinomial runs')
  
  # ggplot(joint_df, aes(x=fullRE_M.beta_gamma_shape, y=as.numeric(M.converged), group=fullRE_M.beta_gamma_shape, col=M.d))+
  #   geom_violin()+
  #   geom_point()+
  #   # geom_jitter(height = 0.1, alpha=0.8)+
  #   facet_wrap(.~M.n)+
  #   theme_bw()+labs(x='Number of categories (d)', y='Number of successful (1) or unsuccessful (0) convergences')+
  #   ggtitle('Success in convergence of Multinomial runs')
  
  cat('Plotting percentage of successful runs 1\n')
  
  ## get percentage of successful runs
  joint_df_grouping_convergence = joint_df[sapply(unique(joint_df$fullRE_DMSL.idx), function(i) which(joint_df$fullRE_DMSL.idx == i)[1]),] %>%
    dplyr::select(c("fullRE_M.d", "fullRE_M.n", "fullRE_DMSL.converged", "fullRE_M.converged")) %>% 
    group_by(fullRE_M.n, fullRE_M.d) %>%
    mutate(convergence_M=sum(fullRE_M.converged)/(sum(fullRE_M.converged)+sum(!fullRE_M.converged)),
           convergence_DM=sum(fullRE_DMSL.converged)/length(fullRE_DMSL.converged))
  joint_df_grouping_convergence_2 = joint_df[sapply(unique(joint_df$fullRE_DMSL.idx), function(i) which(joint_df$fullRE_DMSL.idx == i)[1]),] %>%
    dplyr::select(c("fullRE_M.d", "fullRE_M.n", "fullRE_DMSL.converged", "fullRE_M.converged", "fullRE_M.beta_gamma_shape")) %>% 
    group_by(fullRE_M.beta_gamma_shape) %>%
    mutate(convergence_M=sum(fullRE_M.converged)/(sum(fullRE_M.converged)+sum(!fullRE_M.converged)),
           convergence_DM=sum(fullRE_DMSL.converged)/length(fullRE_DMSL.converged))
  joint_df_grouping_convergence_2_n = joint_df[sapply(unique(joint_df$fullRE_DMSL.idx), function(i) which(joint_df$fullRE_DMSL.idx == i)[1]),] %>%
    dplyr::select(c("fullRE_M.d", "fullRE_M.n", "fullRE_DMSL.converged", "fullRE_M.converged", "fullRE_M.beta_gamma_shape")) %>% 
    group_by(fullRE_M.beta_gamma_shape, fullRE_M.n) %>%
    mutate(convergence_M=sum(fullRE_M.converged)/(sum(fullRE_M.converged)+sum(!fullRE_M.converged)),
           convergence_DM=sum(fullRE_DMSL.converged)/length(fullRE_DMSL.converged))
  joint_df_grouping_convergence_2_d = joint_df[sapply(unique(joint_df$fullRE_DMSL.idx), function(i) which(joint_df$fullRE_DMSL.idx == i)[1]),] %>%
    dplyr::select(c("fullRE_M.d", "fullRE_M.n", "fullRE_DMSL.converged", "fullRE_M.converged", "fullRE_M.beta_gamma_shape")) %>% 
    group_by(fullRE_M.beta_gamma_shape, fullRE_M.d) %>%
    mutate(convergence_M=sum(fullRE_M.converged)/(sum(fullRE_M.converged)+sum(!fullRE_M.converged)),
           convergence_DM=sum(fullRE_DMSL.converged)/length(fullRE_DMSL.converged))
  
  cat('Plotting percentage of successful runs 2\n')
  try({
  pdf(paste0(flder_out, generation, "/summaries/M_DM_convergence.pdf"))
  grid.arrange(ggplot(joint_df_grouping_convergence, aes(x=fullRE_M.d, y=convergence_M, col=fullRE_M.n, group=fullRE_M.n))+geom_line()+geom_point()+ggtitle('Convergence Multinomial'),
  ggplot(joint_df_grouping_convergence, aes(x=fullRE_M.d, y=convergence_DM, col=fullRE_M.n, group=fullRE_M.n))+geom_line()+geom_point()+ggtitle('Convergence Dirichlet-Multinomial'))
  dev.off()
  })
  
  cat('Plotting percentage of successful runs 3\n')
  
  try({
  grid.arrange(ggplot(joint_df_grouping_convergence_2, aes(x=fullRE_M.beta_gamma_shape, y=convergence_M, group=fullRE_M.beta_gamma_shape))+geom_line()+
                 geom_point()+ggtitle('Convergence Multinomial'),
               ggplot(joint_df_grouping_convergence_2, aes(x=fullRE_M.beta_gamma_shape, y=convergence_DM, group=fullRE_M.beta_gamma_shape))+geom_line()+
                 geom_point()+ggtitle('Convergence Dirichlet-Multinomial'))
  })
  # grid.arrange(ggplot(joint_df_grouping_convergence_2_n, aes(x=fullRE_M.beta_gamma_shape, y=convergence_M, col=M.n, group=interaction(M.n)))+geom_line()+
  #                geom_point()+ggtitle('Convergence Multinomial'),
  #              ggplot(joint_df_grouping_convergence_2_n, aes(x=fullRE_M.beta_gamma_shape, y=convergence_DM, col=M.n,  group=interaction(M.n)))+geom_line()+
  #                geom_point()+ggtitle('Convergence Dirichlet-Multinomial'))
  # grid.arrange(ggplot(joint_df_grouping_convergence_2_d, aes(x=fullRE_M.beta_gamma_shape, y=convergence_M, col=M.d, group=interaction(M.d)))+geom_line()+
  #                geom_point()+ggtitle('Convergence Multinomial'),
  #              ggplot(joint_df_grouping_convergence_2_d, aes(x=fullRE_M.beta_gamma_shape, y=convergence_DM, col=M.d,  group=interaction(M.d)))+geom_line()+
  #                geom_point()+ggtitle('Convergence Dirichlet-Multinomial'))
  # 
  
  ##power
  
  cat('Plotting power\n')
  
  try({
  ggplot(varying_n_betashape %>% dplyr::select(Power, model, n, beta_gamma_shape),
         aes(x=beta_gamma_shape+0.001, col=n, y=Power, group=n))+
    geom_line()+geom_point()+facet_wrap(.~model)+
    theme_bw()+scale_x_continuous(trans = "log10")
  })
  
  cat('Plotting power 2\n')
  
  try({
  ggplot(varying_n_betashape %>% dplyr::select(Power, model, n, beta_gamma_shape),
         aes(x=beta_gamma_shape+0.001, col=model, y=Power, group=model))+
    geom_line()+geom_point()+facet_wrap(.~n, nrow=1)+
    theme_bw()+scale_x_continuous(trans = "log10")+theme(legend.position = "bottom")+
    scale_color_jcolors(palette = "pal8")
  ggsave(paste0(flder_out, generation, "/summaries/power_facet.pdf"), width = 6.5, height = 3)
  })
  
  try({
  ggplot(varying_n, aes(x=Sensitivity, y=1-Specificity, group=model, col=model))+geom_step()+theme_bw()
  })
  
  try({
  plot(varying_n$Sensitivity, varying_n$Specificity)
  })
  
  ## comparison of ilr and test of  proportions when removing the first column
  
  cat('comparison of ilr and test of  proportions when removing the first column\n')
  
  try({
  small_num <- 0.0001
  ggplot(pvals_data_frame, aes(x=ttest_props+small_num,
                               y=ttest_ilr_adj+small_num))+geom_point()+theme_bw()+scale_x_continuous(trans = "log2")+
    scale_y_continuous(trans = "log2")+
    geom_vline(xintercept = (0.05+small_num), lty='dashed', col='blue')+
    geom_hline(yintercept = (0.05+small_num), lty='dashed', col='blue')
  
  runs_ttest_props
  })
  
  
  ## looking at results that have been created using the same set of parameters
  
  ### NOTE! because this is slightly confusing
  ## these are not datasets: the 01 indicates the first beta of dataset 0!
  # rownames(runs_fullREM0)[grep("_dataset04", rownames(runs_fullREM0))]
  # rownames(datasets_files)[grep("_dataset04", rownames(datasets_files))]
  # 
  # ## we don't have any replicates at the moment
  # rownames(runs_fullREM0)[grep("_dataset1", rownames(runs_fullREM0))] ## (no?) runs with dataset #2
  # rownames(runs_fullREM0)[grep("_dataset2", rownames(runs_fullREM0))] ## (no?) runs with dataset #2
  # 
  # 
  # runs_diagREDM0[grep("multiple_GenerationJnorm_50_100_80_7_0.6_NA_NA_NA_dataset", rownames(runs_diagREDM0)),]
  # 
  # rownames(runs_fullREM)
  # table(gsub('.{0,1}$', '', gsub("^.*\\_","", rownames(runs_fullREM))))
  # 
  # ## if there are multiple runs
  # 
  # try({
  #   
  # ggplot(runs_fullREM0[grepl(sub("_[^_]+$", "", rownames(runs_fullREM0)[grep("_dataset1", rownames(runs_fullREM0))][1]),
  #       rownames(runs_fullREM0)),], aes(x=idx, y=beta_est))+geom_point()+facet_wrap(.~idx_within_dataset)
  # })
  # 
  cat('Plotting weighted accuracy of converged\n')
  try({
  a <- ggplot(varying_n_all_converged, aes(x=n, y = WeightedAccuracy, col=model, group=model))+
    geom_point()+geom_line()+theme_bw()+
    scale_color_manual(values=colours_models)+ggtitle(generation)+
    guides(col=guide_legend(nrow=2))+
    theme(legend.position = "bottom")+
    labs(col=NULL)
  
    legend <- cowplot::get_legend(a)
    pdf(paste0(flder_out, "legend_models.pdf"), height = 1, width = 6.5)
    grid.newpage()
    grid.draw(legend)
    dev.off()
  })
  
  ## save results for tests
  cat('Saving objects\n')
  
  object_save <- list(generation=generation,
                      DA_bool=DA_bool,
                      runs_ttest_irl=runs_ttest_irl,
                      runs_ttest_props=runs_ttest_props,
                      pvals_runs_HMP=pvals_runs_HMP,
                      pvals_runs_HMP2=pvals_runs_HMP2,
                      pvals_ttest_ilr=pvals_ttest_ilr,
                      pvals_perturbation=pvals_perturbation,
                      pvals_permutation=pvals_permutation,
                      Harris_chi=Harris_chi,
                      joint_df=joint_df,
                      datasets=datasets,
                      pvals_data_frame=pvals_data_frame)
  
  saveRDS(object = object_save,
          file = paste0("../../../../data/assessing_models_simulation/summaries_synthetic_DA/",
                        generation, ".RDS"))
  
  
  ## we have much fewer datasets than runs
  # DA_bool %>% length
  # nrow(runs_fullREM0)
  # 
  # table(gsub(".*_", "", names(DA_bool)))
  # table(gsub(".*_", "", rownames(runs_fullREM0)))
  # 
  # length(runs_ttest_props)
  # ggplot(varying_n_all_converged, aes(x=n, y = WeightedAccuracy, col=model, group=model))+geom_point()+geom_line()+theme_bw()+
  #   scale_color_manual(values=colours_models)+guides(col=FALSE)+ggtitle(generation)+
  #   geom_label_repel(data = varying_n_all_converged %>% dplyr::filter(n == max(n)),
  #              aes(x=n, y=WeightedAccuracy, label=model), alpha=0.6, size=3)
  # ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_N_all_converged_palette2_annotated.pdf"),
  #        height = 3.0, width = 4.0)
  # # 
  # # ggplot(varying_d, aes(x=model, y = AUC, group=model, col=d))+geom_boxplot()+geom_jitter()+theme_bw()
  # # ggplot(varying_n, aes(x=model, y = AUC, group=model, col=n))+geom_boxplot()+geom_jitter()+theme_bw()
  # 
  # 
  sort_first_col_by_second <- function(x) unlist(x[order(x[,2], decreasing = T),1])
  sort_first_col_by_second(varying_n %>% dplyr::group_by(model) %>% dplyr::summarise(median(WeightedAccuracy)))
  # 
  
  cat('Boxplots of weighted accuracy by n\n')
  ggplot(varying_n, aes(x=factor(model,
                                 levels=sort_first_col_by_second(varying_n %>%
                                         dplyr::group_by(model) %>%
                                         dplyr::summarise(median(WeightedAccuracy, na.rm = T)))),
                        y = WeightedAccuracy, group=model, col=model))+geom_boxplot()+geom_jitter()+theme_bw()+
    scale_color_manual(values=colours_models)+ggtitle(generation)+guides(col=FALSE)+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_n_boxplot.pdf"),
         height = 3.0, width = 4)

  cat('Boxplots of weighted accuracy by d\n')
  ggplot(varying_d, aes(x=factor(model,
                                 levels=sort_first_col_by_second(varying_d %>%
                                                                   dplyr::group_by(model) %>%
                                                                   dplyr::summarise(median(WeightedAccuracy, na.rm = T)))),
                        y = WeightedAccuracy, group=model, col=model))+geom_boxplot()+geom_jitter()+theme_bw()+
    scale_color_manual(values=colours_models)+ggtitle(generation)+guides(col=FALSE)+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ggsave(paste0(flder_out, generation, "/summaries/weightedaccuracy_with_d_boxplot.pdf"),
         height = 3.0, width = 4)
  
  if(generation %in% c("GenerationJnormTwoLambdasOneChangingBeta", "GenerationJnormBTwoLambdasOneChangingBeta", "GenerationPois")){
    ##'  only for GenerationJnormTwoLambdasOneChangingBeta: looking at individual betas
    ##'  we want to see if there are more false negatives in samples where the beta_i is lower
    cat('Plotting single category changes\n')
    
    ggplot(data = runs_diagREDM[runs_diagREDM$beta_true != 0,c('beta_true', 'pvals_adj')],
           aes(x=beta_true, y=pvals_adj))+geom_point()+theme_bw()
    runs_diagREDM$beta_est
  
  
    ## this should be normalised by the abundance of the last category which serves as baseline
    ## reminder: it is always the first log-ratio which is not different from zero
    table(runs_diagREDM[runs_diagREDM$idx_within_dataset == 1,]$beta_true)
    table(runs_diagREDM[runs_diagREDM$idx_within_dataset != 1,]$beta_true)
  
    # only_logR_with_change <- runs_diagREDM[runs_diagREDM$idx_within_dataset == 1,]
    cat('Creating<only_logR_with_change>\n')
    only_logR_with_change <- cbind(joint_df[,grepl('diagRE_DM.', colnames(joint_df))], names=joint_df$names)
    colnames(only_logR_with_change) <- gsub("diagRE_DM.", "", colnames(only_logR_with_change))
    sum(only_logR_with_change$idx_within_dataset == 1)
    length(unique(only_logR_with_change$idx))
    length(unique(only_logR_with_change$names))
    only_logR_with_change <- only_logR_with_change[only_logR_with_change$idx_within_dataset == 1,]
    # only_logR_with_change <- only_logR_with_change[!(is.na(only_logR_with_change$names)),]
    dim(only_logR_with_change)
    length(datasets)
    sum(is.na(only_logR_with_change$names))
    names(datasets)[!(gsub(".RDS", "", names(datasets)) %in% only_logR_with_change$names)]
    only_logR_with_change$names[!(only_logR_with_change$names %in% gsub(".RDS", "", names(datasets)))]
    # plot(density(only_logR_with_change$beta_true))
  
    ##' I would expectL for the same beta gamma shape (or, equivalently, the same beta slope true)
    ##' there are more false negatives the lower beta intercept true
    cat('Plotting <only_logR_with_change>\n')
    ggplot(data = only_logR_with_change,
           aes(x=beta_true, y=beta_intercept_true, col=pvals_adj))+geom_point()+theme_bw()
  
    ggplot(data = only_logR_with_change,
           aes(x=beta_true, y=beta_intercept_true, col=pvals_adj < 0.05))+geom_point()+theme_bw()
    ggplot(data = only_logR_with_change,
           aes(x=factor(beta_true), y=beta_intercept_true, col=pvals_adj < 0.05))+geom_point()+theme_bw()
    
    cat('Adding softmax to <only_logR_with_change>\n')
    
    only_logR_with_change$abundance_first_cat[!is.na(only_logR_with_change$names)] = sapply(datasets[(gsub(".RDS", "", names(datasets)) %in% only_logR_with_change$names)], function(i) softmax(c(i$beta[1,],0))[1])
    only_logR_with_change$cut_beta_intercept_true <- cut(only_logR_with_change$beta_intercept_true, breaks = seq(-1, 1, length.out = 10))
    only_logR_with_change_summarised <- only_logR_with_change %>% group_by(beta_true, cut_beta_intercept_true) %>%
      dplyr::summarise(mean_pvals_adj_signif = mean(pvals_adj < 0.05))
    ggplot(data = only_logR_with_change,
           aes(x=factor(beta_true), y=cut_beta_intercept_true,
               col=pvals_adj < 0.05))+geom_point()+theme_bw()
    ggplot(data = only_logR_with_change,
           aes(x=beta_true, y=cut_beta_intercept_true,
               col=pvals_adj < 0.05))+geom_point()+theme_bw()+scale_x_continuous(trans = "log2")
    ggplot(data = only_logR_with_change_summarised,
           aes(x=factor(beta_true), y=cut_beta_intercept_true,
               fill=mean_pvals_adj_signif))+geom_tile()+theme_bw() ## in GenerationJnormBTwoLambdasOneChangingBeta it doesn't depend on the intercept either
    # ggplot(data = only_logR_with_change,
    #        aes(x=factor(beta_true), y=softmax_intercept,
    #            col=pvals_adj < 0.05))+geom_point()+theme_bw()
    # ggplot(data = only_logR_with_change,
    #        aes(x=factor(beta_true), y=cut_beta_intercept_true_softmax,
    #            col=pvals_adj < 0.05))+geom_point()+theme_bw()
    
    ## we compute the absolute abundance (in proportion) of the first category

    # abundances_in_prob <- lapply(unique(runs_diagREDM$idx), function(idx_dataset_it){
    #   softmax(c(runs_diagREDM[runs_diagREDM$idx == idx_dataset_it,]$beta_intercept_true, 0)) ## abundance of first cat
    # })
    # 
    # abundances_first_cat <- sapply(abundances_in_prob, `[`, 1)
    # abundances_last_cat <- sapply(abundances_in_prob, function(i) i[length(i)]) ## not used
    # 
    # length(abundances_first_cat)
    # dim(only_logR_with_change)
    # only_logR_with_change$abundance_first_cat <- abundances_first_cat
    only_logR_with_change$cut_abundances_first_cat <- cut(only_logR_with_change$abundance_first_cat,
                                                          breaks = seq(0, max(only_logR_with_change$abundance_first_cat, na.rm = T), length.out = 4))
    only_logR_with_change$cut_abundances_first_catV2 <- cut(only_logR_with_change$abundance_first_cat,
                                                            breaks = quantile(only_logR_with_change$abundance_first_cat, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = T), include.lowest = T)
    sum(is.na(only_logR_with_change$cut_abundances_first_catV2)) ## should be 0
    
    only_logR_with_change_summarised_2 <- only_logR_with_change %>%
      group_by(beta_true, cut_abundances_first_cat) %>%
      dplyr::summarise(mean_pvals_adj_signif = mean(pvals_adj < 0.05))
    only_logR_with_change_summarised_2 <- only_logR_with_change_summarised_2[!is.na(only_logR_with_change_summarised_2$cut_abundances_first_cat),]
    
    only_logR_with_change_summarised_V2 <- only_logR_with_change %>%
      group_by(beta_true, cut_abundances_first_catV2) %>%
      dplyr::summarise(mean_pvals_adj_signif = mean(pvals_adj < 0.05))
    only_logR_with_change_summarised_V2 <- only_logR_with_change_summarised_V2[!is.na(only_logR_with_change_summarised_V2$cut_abundances_first_catV2),]
    
    if(F){
      .xxxidx <- which(grepl('multiple_GenerationJnormBTwoLambdasOneChangingBeta_20_100_80_6_0.01_NA_NA_NA_dataset1', names(datasets)))
      names(datasets)[.xxxidx]
      datasets[[.xxxidx]]$beta
      runs_diagREDM[grepl("multiple_GenerationJnormBTwoLambdasOneChangingBeta_20_100_80_6_0.01_NA_NA_NA_dataset11", rownames(runs_diagREDM)),]
      
      a1 <- readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnormBTwoLambdasOneChangingBeta_20_100_80_6_0.01_diagREDM_NA_NA_NA_dataset1.RDS"))
      plot_betas(a1) ## this is differentially abundant but not because of the actually DA run
      createBarplot(give_dummy_rownames(give_dummy_colnames(normalise_rw(datasets[[.xxxidx]]$W))))
      
    }
    
    ggplot(data = only_logR_with_change,
           aes(x=factor(beta_true), y=abundance_first_cat,
               col=pvals_adj < 0.05))+geom_point()+theme_bw()
    ggplot(data = only_logR_with_change,
           aes(x=factor(beta_true), y=abundance_first_cat,
               col=pvals_adj < 0.05))+geom_point()+theme_bw()+facet_wrap(.~n)
    ggplot(data = only_logR_with_change,
           aes(x=factor(beta_true), y=abundance_first_cat,
               col=pvals_adj < 0.05))+geom_point()+theme_bw()+facet_wrap(.~d)
    
    tikz(paste0(flder_out, generation, "/summaries/intercept_and_pvals.tex"),
         height = 3.5, width = 6)
    ggplot(data = only_logR_with_change_summarised_2,
           aes(x=factor(beta_true), y=cut_abundances_first_cat,
               fill=mean_pvals_adj_signif))+geom_tile()+theme_bw()+
      scale_fill_jcolors_contin(palette = "pal2")+
      labs(fill='Fraction of runs of DA', x='Beta slope of changing logR', y='Discretised abundance of DA category')
    # ## in GenerationJnormBTwoLambdasOneChangingBeta it doesn't depend on the intercept either
    dev.off()
  
    ggplot(data = only_logR_with_change_summarised_2,
           aes(x=factor(beta_true), y=cut_abundances_first_cat,
               fill=mean_pvals_adj_signif))+geom_tile()+theme_bw()+
      scale_fill_jcolors_contin(palette = "pal3")+
      labs(fill='Fraction of runs of DA', x='Beta slope of changing logR', y='Discretised abundance of DA category')
    ggsave(paste0(flder_out, generation, "/summaries/intercept_and_pvals.pdf"),
           height = 3.5, width = 6)
   
    ## different categorisation
    ggplot(data = only_logR_with_change_summarised_V2,
           aes(x=factor(beta_true), y=cut_abundances_first_catV2,
               fill=mean_pvals_adj_signif))+geom_tile()+theme_bw()+
      scale_fill_jcolors_contin(palette = "pal3")+
      labs(fill='Fraction of runs of DA', x='Beta slope of changing logR', y='Discretised abundance of DA category')
    ggsave(paste0(flder_out, generation, "/summaries/intercept_and_pvalsV2.pdf"),
           height = 3.5, width = 6)
  
  }
  # 
  # table(only_logR_with_change$beta_true)
  # 
  # runs_diagREDM
  # 
  # ggplot(data = only_logR_with_change,
  #        aes(x=mean(pvals_adj < 0.05, na.rm = T), y=cut_beta_intercept_true))+
  #          geom_point()+theme_bw()
  # 
  # ggplot(data = only_logR_with_change,
  #        aes(x=as.numeric(pvals_adj < 0.05), y=cut_beta_intercept_true))+
  #   geom_violin()+theme_bw()
  # 
  # plot(joint_df$fullRE_M.beta_intercept_true,
  #      joint_df$fullRE_M.beta_true)
  # 
  # plot(joint_df$fullRE_M.beta_intercept_true,
  #      joint_df$fullRE_M.beta_gamma_shape)
  # 
  # ## including all the logR, in GenerationJnormTwoLambdasOneChangingBeta/ GenerationJnormBTwoLambdasOneChangingBeta
  # plot(joint_df$fullRE_M.beta_true,
  #      joint_df$fullRE_M.beta_gamma_shape)
  # 
  # ## including only logR that changes, in GenerationJnormTwoLambdasOneChangingBeta/ GenerationJnormBTwoLambdasOneChangingBeta
  # plot(joint_df[joint_df$fullRE_M.idx_within_dataset == 1,]$fullRE_M.beta_true,
  #      joint_df[joint_df$fullRE_M.idx_within_dataset == 1,]$fullRE_M.beta_gamma_shape)

}else{
  warning('There were no common converged runs\n')
  saveRDS(object = NA,
          file = paste0("../../../../data/assessing_models_simulation/summaries_synthetic_DA/",
                        generation, ".RDS"))
}


