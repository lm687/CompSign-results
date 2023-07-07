## Creating DA datasets using mixtures, as in Holmes
## here no random effects

rm(list = ls())

library(optparse)
library(uuid)
source("2_inference/helper/helper_DA_stan.R")
source("1_create_ROO/roo_functions.R")
source("1_create_ROO/helper_1_create_ROO.R")
source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("2_inference_TMB/helper_TMB.R") ## for softmax

option_list = list(
  make_option(c("--input"), type="character", default=NA,
              help="Text with small description of the type of simulation being carried out", metavar="character"),
  make_option(c("--d"), type="numeric", default=NA,
              help="Number of features", metavar="numeric"),
  make_option(c("--n"), type="numeric", default=NA,
              help="Number of samples", metavar="numeric"),
  make_option(c("--nlambda"), type="numeric", default=NA,
              help="Parameter lambda for Poisson draws of number of mutations in sample", metavar="numeric"),
  make_option(c("--beta_gamma_shape"), type="numeric", default=NA,
              help="Fraction of additional mixture for differential abundance (from 0 to 1)", metavar="numeric"),
  make_option(c("--lambda"), type="numeric", default=0,                            ## not used
              help="Overdispersion parameter", metavar="numeric"),
  make_option(c("--outfile"), type="character", default=NA,
              help="Output file in which to write the dataset (RDS file)", metavar="character"),
  make_option(c("--beta_intercept_input"), type="character", default=NA, 
              help="Used to specify multinomial proportions for the first (base) mixture population", metavar="character"),
  make_option(c("--beta_slope_input"), type="character", default=NA,
<<<<<<< HEAD
              help="Used to specify multinomial proportions for the first (differential abundance) mixture population", metavar="character"),
=======
              help="Used to specify multinomial proportions for the second (differential abundance) mixture population", metavar="character"),
>>>>>>> b7516544d6581da5bf0a960e309788c1fba6dff6
  make_option(c("--sdRE_input"), type="character", default=NA,
              help="Fixed standard deviations and covariances for RE", metavar="character") 
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## (A) simulate from the model with different parameters
d = opt$d ## number of signatures
n = opt$n ## number of samples
beta_gamma_shape = opt$beta_gamma_shape ## shape parameter for the beta
if(!( (beta_gamma_shape >= 0) & (beta_gamma_shape <= 1) )){stop('beta_gamma_shape bust be between 0 and 1')}

Nm_lambda = opt$nlambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)


sim_props_mix1 = F
if(is.null(opt$beta_intercept_input)){
  cat('Simulating beta intercept\n')
  sim_props_mix1 = T
}else{
  if(is.na(opt$beta_intercept_input)){
    cat('Simulating beta intercept\n')
    sim_props_mix1 = T
  }else{
    if(opt$beta_intercept_input == 'NA'){
      sim_props_mix1 = T
    }else{   
      cat('Reading input beta intercept file')
      beta[1,] = readRDS(opt$beta_intercept_input)
      sim_props_mix1 = F
    }
  }
}

if(sim_props_mix1)   props_mix1 <- MCMCpack::rdirichlet(1, rep(1/d, d))


sim_props_mix2 = F
if(is.null(opt$beta_slope_input)){
  cat('Simulating beta slope\n')
  sim_props_mix2 = T
}else{
  if(is.na(opt$beta_slope_input)){
    sim_props_mix2 = T
  }else{
    if(opt$beta_slope_input == 'NA'){
      sim_props_mix2  = T	
    }else{
      sim_props_mix2 = F
      cat('Reading input beta slope file')
      beta[2,] = readRDS(opt$beta_slope_input)
    }
  }
}

if(sim_props_mix2)  props_mix2 <- MCMCpack::rdirichlet(1, rep(1/d, d))

## create alpha
Nm = rpois(n*2, lambda = Nm_lambda) # number of mutations per sample \in N^{2*n}
## split the total counts for each of the two populations
Nm_mix2 <- round(Nm*beta_gamma_shape)
## the first group is a pure population
Nm_mix2[1:n] <- 0
Nm_mix1 <- Nm-Nm_mix2

## draw counts
W <- t(sapply(1:(n*2), function(idx) rmultinom(n = 1, size = Nm_mix1[idx], prob = props_mix1)+
                                      rmultinom(n = 1, size = Nm_mix2[idx], prob = props_mix2)))

dim(W)

## Save as object so that we can perform the inference
objects_counts <- new("exposures_cancertype",
                      cancer_type="simulated data",
                      type_classification = "simulated two group",
                      number_categories = 2,
                      id_categories = c('sim1', 'sim2'),
                      active_signatures = "absent; simulation",
                      count_matrices_all = list(give_dummy_col_names(give_dummy_row_names(W[1:n,])), 
                                                give_dummy_col_names(give_dummy_row_names(W[(n+1):(2*n),]))),
                      count_matrices_active = list(list(), list()),
                      sample_names = rownames(give_dummy_row_names(W[1:n,])),
                      modification = "none",
                      is_null_active = TRUE,
                      is_empty="Non-empty"
)

sd_RE <- NA
lambda <- NA
X_sim <- NA
Z_sim <- NA
u <- NA
lambdas <- NA
alphabar <- NA
alpha <- NA
beta <- NA

uuid = uuid::UUIDgenerate()
write.table(paste0("3_analysis/helper/table_simulation_params_", uuid, ".txt"), append = FALSE, x = cbind('d', 'n', 'beta_gamma_shape', 'sd_RE', 'Nm_lambda'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(paste0("3_analysis/helper/table_simulation_params_", uuid, ".txt"), append = TRUE, x = cbind(d, n, beta_gamma_shape, sd_RE, Nm_lambda), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)


saveRDS(list(objects_counts=objects_counts, d=d, n= n, beta_gamma_shape=beta_gamma_shape, sd_RE=sd_RE, lambda=lambda, Nm_lambda=Nm_lambda,
             X_sim = X_sim, beta = beta, Z_sim = Z_sim, u = u, lambdas = lambdas, alphabar = alphabar, alpha = alpha, Nm = Nm, W = W),
        file = opt$outfile)


