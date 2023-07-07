## Creating dataset to assess the performance of models for parameter recovery
## Same as Generation C norm, but with random effects for all log-ratios

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
              help="Shape parameter for gamma distribution for beta (i.e. slope coefficient for changes in exposure between groups)", metavar="numeric"),
  make_option(c("--lambda"), type="numeric", default=0,
              help="Overdispersion parameter", metavar="numeric"),
  make_option(c("--beta_intercept_input"), type="character", default=NA,
              help="Fixed intercept for the betas", metavar="character"),
  make_option(c("--beta_slope_input"), type="character", default=NA,
              help="Fixed slope for the betas", metavar="character"),
  make_option(c("--sdRE_input"), type="character", default=NA, ## not used
              help="Fixed slope for the betas", metavar="character"), ## not used
  make_option(c("--outfile"), type="character", default=NA,
              help="Output file in which to write the dataset (RDS file)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## (A) simulate from the model with different parameters
d = opt$d # opt$Nk ## number of signatures
n = opt$n #opt$Ns ## number of samples
beta_gamma_shape = opt$beta_gamma_shape  ##opt$hyperparam_shape ## shape parameter for the beta

sd_RE = runif(n = d-1, min = 0, max = 1) ## standard deviation for random effects

lambda = rep(opt$lambda, 2) ## overdispersion scalars. Lower value -> higher overdispersion
Nm_lambda = opt$nlambda ## opt$Nm_lambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)

## Group effects
## covariate matrix
X_sim = matrix(NA, nrow=2, ncol=2*n)
## the samples are split into two groups
X_sim[1,] = 1
X_sim[2,] = rep(c(0,1), each=n)
beta = matrix(0, nrow=2, ncol=d-1)
beta[1,] = runif(n = d-1, min = -1, max = 1)
if(beta_gamma_shape == 0){
  ## if non-differentially abundant, make it truly non-differentially abundant, i.e. exactly zero
  beta[2,] = 0
}else{
  beta[2,] = rnorm(n = d-1, mean = beta_gamma_shape, sd = .6) ## for the slope coefficients
}

## Random effects
Z_sim0 = matrix(0, nrow=n, ncol=n)
diag(Z_sim0) = 1
Z_sim = t(rbind(Z_sim0, Z_sim0))

## independent random effects
u = sapply(sd_RE, rnorm, n = n, mean = 0)

## lambdas: overdispersion
lambdas = c(rep(lambda[1], n), rep(lambda[2], n))

## create alpha
alphabar = softmax( cbind(t(X_sim)%*%beta + t(Z_sim)%*%u, 0) )
alpha = alphabar * lambdas

## Create the counts
Nm = rpois(n*2, lambda = Nm_lambda) # number of mutations per sample \in N^{2*n}
W = matrix(NA, nrow = 2*n, ncol = d)
for(l in 1:(2*n)){
  W[l,] = HMP::Dirichlet.multinomial(Nrs = Nm[l], shape = alpha[l,])
}

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

uuid = uuid::UUIDgenerate()
write.table(paste0("3_analysis/helper/table_simulation_params_", uuid, ".txt"), append = FALSE, x = cbind('d', 'n', 'beta_gamma_shape', 'sd_RE', 'Nm_lambda'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(paste0("3_analysis/helper/table_simulation_params_", uuid, ".txt"), append = TRUE, x = cbind(d, n, beta_gamma_shape, sd_RE, Nm_lambda), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

saveRDS(list(objects_counts=objects_counts, d=d, n= n, beta_gamma_shape=beta_gamma_shape, sd_RE=sd_RE, lambda=lambda, Nm_lambda=Nm_lambda,
             X_sim = X_sim, beta = beta, Z_sim = Z_sim, u = u, lambdas = lambdas, alphabar = alphabar, alpha = alpha, Nm = Nm, W = W),
        # file = paste0("../data/assessing_models_simulation/datasets/ ", uuid, ".RDS"))
        file = opt$outfile)


