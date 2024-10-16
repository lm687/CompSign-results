
rm(list = ls())

library(TMB)
library(optparse)
library(CompSign)

# source("1_create_ROO/roo_functions.R")
# source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
# source("2_inference_TMB/mm_multinomial/helper_functions.R")
# source("2_inference_TMB/helper_TMB.R")

debug <- F
if(debug){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("../") ## from where snakemake is run
  
  opt <- list()
  opt$input = '../data/roo/Liver-HCC_signaturesMSE_ROO.RDS'
  opt$output = '../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDM_Liver-HCC_signaturesMSE.RDS'
  opt$model = 'diagREDM' 
  opt$feature_type = 'signaturesMSE' 
  opt$optimiser = 'nlminb' 
  opt$simulation_bool = F 
  opt$read_directly = T 
  opt$use_previous_run_startingvals  = T
  
  opt <- list()
  opt$read_directly = T 
  opt$simulation_bool = T
  opt$input = '../data/assessing_models_simulation/datasets/multiple_GenerationMixturewithCCFb_68_NA_500_96_0.2_NA_NA_NA/multiple_GenerationMixturewithCCFb_68_NA_500_96_0.2_NA_NA_NA_dataset1.RDS'
  opt$output = '../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationMixturewithCCFb_68_NA_500_96_0.2_diagREDM_NA_NA_NA/multiple_GenerationMixturewithCCFb_68_NA_500_96_0.2_diagREDM_NA_NA_NA_dataset1.RDS'
  opt$model = 'diagREDM'
  opt$optimiser = 'nlminb'
  opt$use_previous_run_startingvals = T
  opt$feature_type = "signatures"
  
}


option_list = list(
  make_option(c("--model"), type="character", default=NA,
              help="Which model to use for inference", metavar="character"),
  make_option(c("--input"), type="character", default=NA,
              help="Input file with dataset (RDS)", metavar="character"),
  make_option(c("--simulation_bool"), type="logical", default=T,
              help="Is ct the name of the file to read?", metavar="logical"),
  make_option(c("--read_directly"), type="logical", default=T,
              help="Should the opt$input file be read directly with load_PCAWG?", metavar="logical"),
  make_option(c("--feature_type"), type="character", default="signatures",
              help="Type of feature: signatures, signaturesPCAWG, etc", metavar="character"),
  make_option(c("--output"), type="character", default=NA,
              help="Output file in which to write the results of the inference (RDS file)", metavar="character"),
  make_option(c("--optimiser"), type="character", default="optim",
              help="Which optimiser to use", metavar="character"),
  make_option(c("--nonexo_bool"), type="logical", default=F,
              help="Should only nonexogenous signatures be selected?", metavar="logical"),
  make_option(c("--use_previous_run_startingvals"), type="logical", default=F,
              help="Should we use any previous run, if available, to use the estimated values as starting values? The previous run should be called the same as the output, but with _NC.RDS instead of .RDS", metavar="logical"),
  make_option(c("--return_opt_bool"), type="logical", default=F,
              help="Should the RDS file contain the opt (TRUE for LRT, FALSE OTHERWISE)", metavar="logical")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source("../../../CompSign-do-not-modify_do_everything_online/R/DA_functions.R")


## by default
simulation_bool = opt$simulation_bool

if(opt$optimiser == "optim"){
  use_nlminb = F
}else{
  use_nlminb = T
}

system(paste0('mkdir -p ', dirname(opt$output)))

cat('Model:', opt$model, '\n')
cat('Feature type:', opt$feature_type, '\n')
cat('Using nlminb:', use_nlminb, '\n')
cat('Simulation boolean:', simulation_bool, '\n')

if(opt$nonexo_bool | grepl('nonexo', opt$output)){
  opt$model <- gsub("wSBS1SBS5nonexo", "", opt$model)
  opt$model <- gsub("nonexo", "", opt$model)
}
cat('Model:', opt$model, '\n')

stopifnot(opt$model == 'fullREM')


if(opt$model == "fullREM"){
  model = 'fullRE_M'
}

cat(opt$input)
a = readRDS(opt$input)

a$exposures_bins_QP_list$x <- as.matrix(a$exposures_bins_QP_list$x)

a$exposures_bins_QP_list$x <- cbind(a$exposures_bins_QP_list$x,
                                    baselinepatient=as.vector(as.numeric(a$exposures_bins_QP_list$x[,'bin_order'] == 1)))
a$exposures_bins_QP_list$x <- (a$exposures_bins_QP_list$x[,c('mean_ccf', 'baselinepatient')])


time1=Sys.time() 
results_inference = try({wrapper_run_TMB(a$exposures_bins_QP_list, model=model)})
time2=Sys.time() 

saveRDS(object = results_inference, file = opt$output)
