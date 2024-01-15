rm(list = ls())

library(optparse)

debug <- F
if(debug){
  opt <- list()
  opt$filestring = '../data/assessing_models_simulation/inference_results/TCSM/GenerationMixtureSimulation_150_200_NA_NA_-3_TCSM_NA_NA_NA/'
  opt$it = 4
  opt$input = '../data/assessing_models_simulation/inference_results/TCSM/GenerationMixtureSimulation_150_200_NA_NA_-3_TCSM_NA_NA_NA/exposures4TEMP'
  opt$output = '../data/assessing_models_simulation/inference_results/TCSM/multiple_GenerationMixtureSimulation_150_200_NA_NA_-3_TCSM_NA_NA_NA/multiple_GenerationMixtureSimulation_150_200_NA_NA_-3_TCSM_NA_NA_NA_dataset4.RDS'

}

option_list = list(
  make_option(c("--filestring"), type="character", default=NA,
              help="Which model to use for inference", metavar="character"),
  make_option(c("--it"), type="numeric", default=NA,
              help="iteration", metavar="numeric"),
  make_option(c("--input"), type="character", default=NA,
              help="Input file with dataset (RDS)", metavar="character"),
  make_option(c("--output"), type="character", default=NA,
              help="Output file in which to write the results of the inference (RDS file)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

fles_list <- list(covariate=paste0(opt$filestring, '/', 'covariate', opt$it, 'TEMP'),
                  exposures=paste0(opt$filestring, '/', 'exposures', opt$it, 'TEMP'),
                  signatures=paste0(opt$filestring, '/', 'signatures', opt$it, 'TEMP'),
                  effect=paste0(opt$filestring, '/', 'effect', opt$it, 'TEMP'),
                  sigma=paste0(opt$filestring, '/', 'sigma', opt$it, 'TEMP'),
                  gamma=paste0(opt$filestring, '/', 'gamma', opt$it, 'TEMP'),
                  EstSignificance=paste0(opt$filestring, '/', 'EstSignificance', opt$it, 'TEMP'))
first_row_to_colnames <- function(i){
  x <- data.frame(i[-1,])
  colnames(x) <- i[1,]
  x
}

first_col_to_rownames <- function(i){
  x <- data.frame(i[,-1])
  rownames(x) <- i[,1]
  x
}

joint_output_TCSM <- sapply(fles_list, read.table, sep='\t')
joint_output_TCSM$covariate <- first_col_to_rownames(first_row_to_colnames(joint_output_TCSM$covariate))
joint_output_TCSM$exposures <- first_col_to_rownames(first_row_to_colnames(joint_output_TCSM$exposures))
joint_output_TCSM$EstSignificance <- first_col_to_rownames(first_row_to_colnames(joint_output_TCSM$EstSignificance))
saveRDS(joint_output_TCSM, file = opt$output)
for(i in fles_list){
  cat(paste0('rm ', i))
  cat('\n')
  system(paste0('rm ', i))
}