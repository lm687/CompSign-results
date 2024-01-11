## load the dataset, or simulated dataset, and run HiLDA
# load(system.file("extdata/sample.rdata", package="HiLDA"))

rm(list = ls())

library(TMB)
library(optparse)
library(HiLDA)
# library(CompSign)

# source("1_create_ROO/roo_functions.R")
# source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
# source("2_inference_TMB/mm_multinomial/helper_functions.R")
# source("2_inference_TMB/helper_TMB.R")

debug <- F
if(debug){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("../../../")
  opt <- list()
  opt$input = '../data/assessing_models_simulation/datasets/multiple_GenerationMixtureSimulation_100_200_NA_NA_0.5_NA_NA_NA/multiple_GenerationMixtureSimulation_100_200_NA_NA_0.5_NA_NA_NA_dataset9.RDS'
  opt$output = '../data/assessing_models_simulation/inference_results/HiLDA/multiple_GenerationMixtureSimulation_100_200_NA_NA_0.5_HiLDA_NA_NA_NA/multiple_GenerationMixtureSimulation_100_200_NA_NA_0.5_HiLDA_NA_NA_NA_dataset9.RDS'
  opt$model = 'HiLDA'
}

option_list = list(
  make_option(c("--model"), type="character", default=NA,
              help="Which model to use for inference", metavar="character"),
  make_option(c("--input"), type="character", default=NA,
              help="Input file with dataset (RDS)", metavar="character"),
  make_option(c("--output"), type="character", default=NA,
              help="Output file in which to write the results of the inference (RDS file)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tcsm_path = '/Users/lenamorrill/Software/tcsm/'

cat('Model:', opt$model, '\n')


# flder <- "../../CompSign-results/data/assessing_models_simulation/datasets/"
# flder <- list.files(flder, full.names = T)
# flder <- grep(x = flder, pattern = 'GenerationMixtureSimulation', value = T)
# flder_fles <- unlist(sapply(flder, list.files, full.names = T))

# for(i_file in flder_fles){
  # dataset <- readRDS(i_file)
cat(opt$input)
dataset = readRDS(opt$input) ## using this instead of load_PCAWG as we are interested in loading the input matrices for HiLDA and TSMC


stopifnot(dataset$d == ncol(dataset$W))
  
if(opt$model == 'HiLDA'){
  ## with initial values
  cat('\nRunning HiLDA...\n')
  time1=Sys.time() 
  # results_inference <- 'dummy'
  results_inference <- hildaTest(inputG=dataset$W_HilDA, numSig=dataset$d, refGroup=1:(nrow(dataset$W)/2), nIter=1000,
                          localTest=TRUE)
  time2=Sys.time() 
  cat('... end\n')
}else if(opt$model == 'TCSM'){
  ## write files
  arg_TCSM_mutation_count_file = paste0(dirname(opt$output), '/mutation_countTEMP')
  arg_TCSM_num_sigs = ncol(dataset$W)
  arg_TCSM_covariate = paste0(dirname(opt$output), '/', gsub(".RDS", "", basename(opt$input)), 'covariateTEMP')
  out_TCSM_exposures = paste0(dirname(opt$output), '/', gsub(".RDS", "", basename(opt$input)), 'exposuresTEMP')
  out_TCSM_signatures = paste0(dirname(opt$output), '/', gsub(".RDS", "", basename(opt$input)), 'signaturesTEMP')
  out_TCSM_effect = paste0(dirname(opt$output), '/', gsub(".RDS", "", basename(opt$input)), 'effectTEMP')
  out_TCSM_sigma = paste0(dirname(opt$output), '/', gsub(".RDS", "", basename(opt$input)), 'sigmaTEMP')
  out_TCSM_gamma = paste0(dirname(opt$output), '/', gsub(".RDS", "", basename(opt$input)), 'gammaTEMP')
  out_TCSM_script = paste0(dirname(opt$output), '/', gsub(".RDS", "", basename(opt$input)), 'scriptTEMP.sh')
  rownames(dataset$W_TCSM) = paste0('Subsample', 1:nrow(dataset$W_TCSM))
  write.table(dataset$W_TCSM, file = arg_TCSM_mutation_count_file,
              row.names = T, quote = F, sep = '\t', col.names=NA)
  write.table(data.frame(group=rep(c(0, 1), each=nrow(dataset$W_TCSM)/2), row.names = rownames(dataset$W_TCSM)), file = paste0(dirname(opt$output), '/covariateTEMP'),
              row.names = T, quote = F, sep = '\t', col.names=NA)
  
  library(reticulate)
  use_condaenv("tcsm")
  
  ## conda init; conda activate tcsm;
  writeLines(text = paste0(tcsm_path, 'src/./run_tcsm.R  ', arg_TCSM_mutation_count_file, ' ', arg_TCSM_num_sigs,
                           ' -c=', arg_TCSM_covariate, ' --covariates group', ' --exposures=', out_TCSM_exposures, ' --signatures=',
                           out_TCSM_signatures, ' --effect=', out_TCSM_effect, ' --sigma=', out_TCSM_sigma, ' --gamma=', out_TCSM_gamma),
             con = out_TCSM_script)
  
  ## run TCSM and output to temporary files
  time1=Sys.time() 
  system(readLines(out_TCSM_script))
  time2=Sys.time() 
  
  ## read output files and create object
  
  ## remove temporary files
  # for(i in arg_TCSM_covariate...)
  
}else{
  stop('<model> not implemented')
}

# ## need to get the output file
# saveRDS(hildaLocal, "../../CompSign-results/data/assessing_models_simulation/inference_results/HiLDA/")
# saveRDS(time2-time1, "../../CompSign-results/data/assessing_models_simulation/inference_results/HiLDA/")

system(paste0('mkdir -p ', dirname(opt$output)))
saveRDS(object = results_inference, file = opt$output)

saveRDS(object = time2-time1, file = gsub(".RDS", ".time", opt$output))


# }
