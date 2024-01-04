## modifying the probabilities underlying a multinomial distribution with perturbation, instead of using
## gaussian noise - otherwise identical to Generation K.

rm(list = ls())

library(MCMCpack)
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
              help="Shape parameter for gamma distribution for beta (i.e. slope coefficient for changes in exposure between groups)", metavar="numeric"),
  make_option(c("--lambda"), type="numeric", default=0,
              help="Overdispersion parameter", metavar="numeric"), ## not used
  make_option(c("--beta_intercept_input"), type="character", default=NA,
              help="Fixed intercept for the betas", metavar="character"), ## not used
  make_option(c("--beta_slope_input"), type="character", default=NA,
              help="Fixed slope for the betas", metavar="character"), ## not used
  make_option(c("--sdRE_input"), type="character", default=NA,
              help="Fixed slope for the betas", metavar="character"), ## not used
  make_option(c("--outfile"), type="character", default=NA,
              help="Output file in which to write the dataset (RDS file)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## (A) simulate from the model with different parameters
d = opt$d # opt$Nk ## number of signatures
n = opt$n #opt$Ns ## number of samples
beta_gamma_shape = opt$beta_gamma_shape # parameter for the gaussian noise

## simulate proportions at random
props <- MCMCpack::rdirichlet(n = n, alpha = rep(2/d, d))

## shift for the second group
perturbation_sigs <- rgamma(n = d, rate=beta_gamma_shape, shape =beta_gamma_shape)
shifts_samples <- matrix(sapply(perturbation_sigs, function(pert_it){
  sapply(rnorm(n = n, mean = 1, sd=pert_it/3), function(i) max(c(0.0001,i)))
}), ncol=d)
props2 <- normalise_rw(props*shifts_samples)

#plot(compositions::rcomp(rbind(props, props2)), col=c(rep(1, n), rep(2, n)))

W <- t(apply(rbind(props, props2), 1, function(prop_it) rmultinom(n = 1, size = rpois(1, opt$nlambda), prob = prop_it)))

X_sim = matrix(NA, nrow=2, ncol=2*n)
## the samples are split into two groups
X_sim[1,] = 1
X_sim[2,] = rep(c(0,1), each=n)

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
lambdas <- NA
alphabar <- NA
alpha <- NA
Nm <- NA

uuid = uuid::UUIDgenerate()
write.table(paste0("3_analysis/helper/table_simulation_params_", uuid, ".txt"), append = FALSE, x = cbind('d', 'n', 'beta_gamma_shape', 'sd_RE', 'Nm_lambda'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(paste0("3_analysis/helper/table_simulation_params_", uuid, ".txt"), append = TRUE, x = cbind(d, n, beta_gamma_shape, sd_RE, opt$nlambda), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

saveRDS(list(objects_counts=objects_counts, d=d, n= n, beta_gamma_shape=beta_gamma_shape, sd_RE=sd_RE, lambda=lambda, Nm_lambda=opt$nlambda,
             X_sim = X_sim, beta = beta, Z_sim = NA, u = NA, lambdas = lambdas, alphabar = alphabar, alpha = alpha, Nm = Nm, W = W),
        # file = paste0("../data/assessing_models_simulation/datasets/ ", uuid, ".RDS"))
        file = opt$outfile)


