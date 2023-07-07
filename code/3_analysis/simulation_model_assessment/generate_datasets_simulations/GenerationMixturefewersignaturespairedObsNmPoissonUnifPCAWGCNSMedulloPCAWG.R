## with gaussian noise

rm(list = ls())

library(optparse)
library(uuid)

debug <- F
if(debug){
  opt <- list()
  opt$input = ''
  opt$n = 50
  opt$beta_gamma_shape = -8
  opt$d = 7
  opt$nlambda = 200
  opt$lambda = - 80
  
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("../../../")
  getwd()
  
}else{
  
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
                help="Used to specify multinomial proportions for the second (differential abundance) mixture population", metavar="character"),
    make_option(c("--sdRE_input"), type="character", default=NA,
                help="Fixed standard deviations and covariances for RE", metavar="character") 
  );
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
}

source("2_inference/helper/helper_DA_stan.R")
source("1_create_ROO/roo_functions.R")
source("1_create_ROO/helper_1_create_ROO.R")
source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("2_inference_TMB/helper_TMB.R") ## for softmax

## (A) simulate from the model with different parameters
d = opt$d ## number of signatures
n = opt$n ## number of samples
beta_gamma_shape = opt$beta_gamma_shape ## shape parameter for the beta
cat('beta_gamma_shape: ', beta_gamma_shape, '\n')
cat('beta_gamma_shape in softmax: ', softmax(c(beta_gamma_shape,0)), '\n')

Nm_lambda = opt$nlambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)
# Nm_lambdas <- rpois(n = n*2, lambda = Nm_lambda)

ct_name <- "CNS-Medullo"

HMP_estimates <- readRDS("../data/restricted/pcawg/pcawg_HMP_estimation_results.RDS")
HMP_estimates <- HMP_estimates[[ct_name]]

ct1 <- load_PCAWG(ct_name, typedata = "signaturesPCAWG", simulation = F, path_to_data="../data/")

## select only the four most prevalent signatures
selected_sigs <- names(sort(colSums(ct1$Y), decreasing = T)[1:4])
d <- length(selected_sigs)
ct1 <- give_subset_sigs_TMBobj(ct1, colnames(ct1$Y)[!(colnames(ct1$Y) %in% selected_sigs)])

ct2 <- ct1$Y[ct1$x[,2] == 1,]
ct1 <- ct1$Y[ct1$x[,2] == 0,]

stopifnot(all(rownames(ct2) == rownames(ct1)))

length(colnames(ct1))
colnames(ct2)

## proportions of the mixture
if((beta_gamma_shape == -999)){
  cat('Is infinite\n')
  props_mix1 <- c(1,0)
}else{
  props_mix1 <- softmax(c(0, beta_gamma_shape))
}
print(props_mix1)

## create alpha

# print(Nm_lambdas[1:n])
if(nrow(ct1) < (2*n)){
  warning('Cancer type has fewer observations than runs. Sampling with replacement\n')
  sample_ids <- sample(1:(nrow(ct1)), size = n, replace = T)
}else{
  sample_ids <- sample(1:(nrow(ct1)), size = n)
}


W <- rbind(t(sapply(1:n, function(id_patient){
  rmultinom(1, sum(ct1[sample_ids[id_patient],]), normalise_rw(ct1[sample_ids[id_patient],]))
})),
t(sapply(c(1:n), function(id_patient){
  rmultinom(1, round(sum(ct2[sample_ids[id_patient],])*props_mix1[1]), normalise_rw(ct1[sample_ids[id_patient],]))+
    rmultinom(1, round(sum(ct2[sample_ids[id_patient],])*props_mix1[2]), normalise_rw(ct2[sample_ids[id_patient],]))
})))

W <- rbind(W[1:(n),]+rpois((n)*(d), lambda = 10),
           W[(1+n):(2*n),]+rpois((n)*(d), lambda = 20))

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
             X_sim = X_sim, beta = beta, Z_sim = Z_sim, u = u, lambdas = lambdas, alphabar = alphabar, alpha = alpha, Nm = rowSums(W), W = W),
        file = opt$outfile)


