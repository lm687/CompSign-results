## load the dataset, or simulated dataset, and run HiLDA
# load(system.file("extdata/sample.rdata", package="HiLDA"))

rm(list = ls())

library(optparse)

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

cat(opt$input, '\n')
dataset = readRDS(opt$input) ## using this instead of load_PCAWG as we are interested in loading the input matrices for HiLDA and TCSM


stopifnot(dataset$d == ncol(dataset$W))
  
if(opt$model == 'HiLDA'){
  library(HiLDA)
  
  ## with initial values
  cat('\nRunning HiLDA...\n')
  time1=Sys.time() 
  # results_inference <- 'dummy'
  results_inference <- hildaTest(inputG=dataset$W_HilDA, numSig=dataset$d, refGroup=1:(nrow(dataset$W)/2), nIter=1000,
                          localTest=TRUE)
  time2=Sys.time() 
  cat('... end\n')
  
  system(paste0('mkdir -p ', dirname(opt$output)))
  saveRDS(object = results_inference, file = opt$output)
  saveRDS(object = time2-time1, file = gsub(".RDS", ".time", opt$output))
  
}if(opt$model == 'HiLDAGlobal'){
  library(HiLDA)
  
  ## with initial values
  cat('\nRunning HiLDA...\n')
  time1=Sys.time() 
  results_inference <- hildaTest(inputG=dataset$W_HilDA, numSig=dataset$d, refGroup=1:(nrow(dataset$W)/2), nIter=1000,
                                 localTest=FALSE)
  time2=Sys.time() 
  cat('... end\n')
  
  system(paste0('mkdir -p ', dirname(opt$output)))
  saveRDS(object = results_inference, file = opt$output)
  saveRDS(object = time2-time1, file = gsub(".RDS", ".time", opt$output))
  
}else if(opt$model == 'TCSM'){
  ## write files
  it_num <- gsub(".RDS", "", gsub("dataset", "", gsub("..*_", "", basename(opt$input))))
  arg_TCSM_mutation_count_file = paste0(dirname(opt$output), '/mutationcount', it_num, 'TEMP')
  arg_TCSM_covariate = paste0(dirname(opt$output), '/', 'covariate', it_num, 'TEMP')
  
  system(paste0('mkdir -p ', dirname(arg_TCSM_mutation_count_file)))

  rownames(dataset$W_TCSM) = paste0('Subsample', 1:nrow(dataset$W_TCSM))
  write.table(dataset$W_TCSM, file = arg_TCSM_mutation_count_file,
              row.names = T, quote = F, sep = '\t', col.names=NA)
  write.table(data.frame(group=rep(c(0, 1), each=nrow(dataset$W_TCSM)/2),
                         row.names = rownames(dataset$W_TCSM)), file = arg_TCSM_covariate,
              row.names = T, quote = F, sep = '\t', col.names=NA)

}else if(opt$model %in% c('diagREDMsigextraction1', 'diagREDMsigextraction2', 'diagREDMsigextraction3')){

  ## Using TMB, but first extracting signatures in different ways
  
  source("2_inference_TMB/helper_TMB.R")
  source("1_create_ROO/roo_functions.R")
  source("1_create_ROO/helper_1_create_ROO.R")

  library(CompSign)
  
  signature_definitions <- read.table("../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv", sep=',', header = T)
  signature_definitions_cats <- paste0(signature_definitions$Type, '_', signature_definitions$SubType)
  signature_definitions <- apply(signature_definitions[,-c(1,2)], 2, function(i) i)
  rownames(signature_definitions) <- signature_definitions_cats
  
  ##' re-extracting signatures, and selecting the active signatures using strategy 1

  ## fit all sigs
  W_QP_0 = extract_sigs(dataset$W_muts, dataset$W, fitWsubsetsigs = F)
  
  if(opt$model %in% c('diagREDMsigextraction1', 'diagREDMsigextraction2')){
      ##' select signatures which are 80% of the total sum of exposures
    if(opt$model == c('diagREDMsigextraction1')){
      percentile_arg = 0.8
    }else if(opt$model == c('diagREDMsigextraction2')){
      percentile_arg = 0.7
    }
    sorted_exposures_sum = sort(colSums(W_QP_0), decreasing = T)
    selected_signatures = names(sorted_exposures_sum)[1:(max(which(cumsum(normalise_rw(sorted_exposures_sum)) <= percentile_arg))+1)]
    stopifnot(sum(W_QP_0[,selected_signatures])/sum(W_QP_0) >= percentile_arg)
  }else if(opt$model == c('diagREDMsigextraction3')){
    ## remove signatures which are too sparse
    max_frac_zeros = 0.2 ## if a signature has exposures of 0 in this fraction of samples or more, do not include
    selected_signatures = names(which((colSums(W_QP_0 == 0)/nrow(W_QP_0)) <= max_frac_zeros))
  }
  W_QP = extract_sigs(dataset$W_muts, dataset$W, subset_sigs = selected_signatures)
  stopifnot(colnames(W_QP) == selected_signatures)
  
  
  ## note that each dataset might have a different number of betas
  dataset$objects_counts@count_matrices_all = split_matrix_in_half(W_QP)
  
  time1=Sys.time()
  results_inference = try(wrapper_run_TMB(object = list(Y = W_QP,
                                                        x= cbind(rep(1, nrow(W_QP)), rep(c(0,1), each=(nrow(W_QP)/2))),
                                                        z=give_z_matrix(nrow(W_QP))),
                                          model = 'diagRE_DM', use_nlminb=T,
                                          return_opt=F))
  time2=Sys.time() 
  
  system(paste0('mkdir -p ', dirname(opt$output)))
  saveRDS(object = results_inference, file = opt$output)
  saveRDS(object = time2-time1, file = gsub(".RDS", ".time", opt$output))
  
}else{
  stop('<model> not implemented')
}

# ## need to get the output file
# saveRDS(hildaLocal, "../../CompSign-results/data/assessing_models_simulation/inference_results/HiLDA/")
# saveRDS(time2-time1, "../../CompSign-results/data/assessing_models_simulation/inference_results/HiLDA/")

# }
