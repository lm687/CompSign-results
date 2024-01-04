setwd(dirname(rstudioapi::getSourceEditorContext()$path))

it=200

## Inverting categories
for(i in 0:(it-1)){
  cat('Iteration: ', i, '/', it)
  flename <- paste0("../../../../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_200_14072_80_6_0_betainterceptPCAWG4_betaslopePCAWG4_covmatPCAWG4/multiple_GenerationJnorm_200_14072_80_6_0_betainterceptPCAWG4_betaslopePCAWG4_covmatPCAWG4_dataset", i, ".RDS")
  x <- readRDS(flename)
  attr(x$objects_counts, "count_matrices_all")[[1]] <-   attr(x$objects_counts, "count_matrices_all")[[1]][,ncol(  attr(x$objects_counts, "count_matrices_all")[[1]]):1]
  attr(x$objects_counts, "count_matrices_all")[[2]] <-   attr(x$objects_counts, "count_matrices_all")[[2]][,ncol(  attr(x$objects_counts, "count_matrices_all")[[2]]):1]
  ## modify for bias assessment
  if(is.null(dim(x$sd_RE))){
    ## is a vector
    x$sd_RE <-   NULL ## as we simply have not simulated the data with the original baseline
  }
  x$beta <- NULL  ## as we simply have not simulated the data with the original baseline
  x$u <- NULL  ## as we simply have not simulated the data with the original baseline
  x$alpha <- NULL  ## as we simply have not simulated the data with the original baseline
  x$alphabar <- x$alphabar[,ncol(x$alphabar):1] ## this we know as it is in the simplex and not LR space
  x$W <- x$W[,ncol(x$W):1]
  flename2 <- gsub('GenerationJnorm', 'GenerationJnormInv', flename)
  system(paste0('mkdir -p ', dirname(flename2)))
  saveRDS(x, file = flename2)
}

## Using most abundant category as baseline

for(i in 0:(it-1)){
  ## Using the category with highest total count as baseline
  cat('Iteration: ', i, '/', it, '\n')
  # flename <- paste0("../../../../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_200_14072_80_6_0_betainterceptPCAWG4_betaslopePCAWG4_covmatPCAWG4/multiple_GenerationJnorm_200_14072_80_6_0_betainterceptPCAWG4_betaslopePCAWG4_covmatPCAWG4_dataset", i, ".RDS")
  # flename <- paste0("../../../../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_200_14072_80_6_0_betainterceptPCAWG4_betaslopezerosPCAWG4_covmatPCAWG4/multiple_GenerationJnorm_200_14072_80_6_0_betainterceptPCAWG4_betaslopezerosPCAWG4_covmatPCAWG4_dataset", i, ".RDS")
  flename <- paste0("../../../../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_200_14072_80_6_0_betainterceptPCAWG4_betaslopeonechangingPCAWG4_covmatPCAWG4/multiple_GenerationJnorm_200_14072_80_6_0_betainterceptPCAWG4_betaslopeonechangingPCAWG4_covmatPCAWG4_dataset", i, ".RDS")

  x <- readRDS(flename)
  max_col <- which.max(colSums(x$W))
  order_cols <- c( (1:ncol(x$W))[-max_col], max_col)
  attr(x$objects_counts, "count_matrices_all")[[1]] <-   attr(x$objects_counts, "count_matrices_all")[[1]][,order_cols]
  attr(x$objects_counts, "count_matrices_all")[[2]] <-   attr(x$objects_counts, "count_matrices_all")[[2]][,order_cols]
  ## modify for bias assessment
  if(is.null(dim(x$sd_RE))){
    ## is a vector
    x$sd_RE <-   NULL ## as we simply have not simulated the data with the original baseline
  }
  x$beta <- NULL  ## as we simply have not simulated the data with the original baseline
  x$u <- NULL  ## as we simply have not simulated the data with the original baseline
  x$alpha <- NULL  ## as we simply have not simulated the data with the original baseline
  x$alphabar <- x$alphabar[,order_cols] ## this we know as it is in the simplex and not LR space
  x$W <- x$W[,order_cols]
  flename2 <- gsub('GenerationJnorm', 'GenerationJnormMax', flename)
  system(paste0('mkdir -p ', dirname(flename2)))
  saveRDS(x, file = flename2)
}
