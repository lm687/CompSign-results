## Creating DA datasets using mixtures, as in Holmes
## using data from PCAWG, and mixing the two groups

rm(list = ls())

library(optparse)
library(uuid)
library(HiLDA) ## for object class MutationFeatureData
source("helper.R")
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

source("1_create_ROO/roo_functions.R")
source("1_create_ROO/helper_1_create_ROO.R")
source("2_inference_TMB/helper_TMB.R") ## for softmax

signature_definitions <- read.table("../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv", sep=',', header = T)
signature_definitions_cats <- paste0(signature_definitions$Type, '_', signature_definitions$SubType)
signature_definitions <- apply(signature_definitions[,-c(1,2)], 2, function(i) i)
rownames(signature_definitions) <- signature_definitions_cats

## (A) simulate from the model with different parameters
# d = opt$d ## number of signatures
n = opt$n ## number of samples
beta_gamma_shape = opt$beta_gamma_shape ## shape parameter for the beta
cat('beta_gamma_shape: ', beta_gamma_shape, '\n')
cat('beta_gamma_shape in softmax: ', softmax(c(beta_gamma_shape,0)), '\n')

Nm_lambda = opt$nlambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)
# Nm_lambdas <- rpois(n = n*2, lambda = Nm_lambda)

ct1 <- load_PCAWG("Lymph-CLL", typedata = "signaturesPCAWG", simulation = F, path_to_data="../data/", override_warning_X_Z=T)

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
colnames(W) <- colnames(ct1)

## W are the exposures. From these exposures, get mutations.
W_muts <- simulate_mutations_from_signatures(W) ## from helper.R
stopifnot(rownames(signature_definitions) == rownames(W))

#----------------------------------------------------------------#
## Re-extract signatures for the fitting the CompSign models. Signature extraction as done for cancer samples
W_QP = t(sapply(1:nrow(W_muts), function(i){
  QPsig(signatures.ref = signature_definitions[,colnames(W)], tumour.ref = rep(colnames(W_muts), W_muts[i,]))
}))
W_QP[W_QP<0] <- 0 ## some values migth be extremely low but negative

## get new exposures
W_QP = t(sapply(1:nrow(W_QP), function(j) table(factor(sample(x = 1:ncol(W_QP), size = rowSums(W)[j], replace = T, prob = W_QP[j,]),
                                                     levels=1:ncol(W_QP)))))
colnames(W_QP) <- colnames(W)
stopifnot(rowSums(W_QP) == rowSums(W))
#----------------------------------------------------------------#

#----------------------------------------------------------------#
## Preparing object for HilDA
six_substitutions_levels <- unique(gsub("_.*", "", rownames(signature_definitions)))
nucleotides_levels <- c('A', 'C', 'G', 'T')
hilda_obj <- matrix_mutations_to_HilDA(t(W_muts))
stopifnot(sum(hilda_obj@countData[3,]) == sum(W_QP))
#----------------------------------------------------------------#

#----------------------------------------------------------------#
## Preparing object for TCSM
W_muts_TCSM <- W_muts
colnames(W_muts_TCSM) <- sapply(colnames(W_muts_TCSM), function(i) paste0(substr(i, 5, 5), '[', substr(i, 1, 3), ']', substr(i, 7, 7)))
#----------------------------------------------------------------#

#----------------------------------------------------------------#
## Save as object so that we can perform the inference
objects_counts <- new("exposures_cancertype",
                      cancer_type="simulated data",
                      type_classification = "simulated two group",
                      number_categories = 2,
                      id_categories = c('sim1', 'sim2'),
                      active_signatures = "absent; simulation",
                      count_matrices_all = list(give_dummy_col_names(give_dummy_row_names(W_QP[1:n,])), 
                                                give_dummy_col_names(give_dummy_row_names(W_QP[(n+1):(2*n),]))),
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

d = ncol(W)
uuid = uuid::UUIDgenerate()
write.table(paste0("3_analysis/helper/table_simulation_params_", uuid, ".txt"), append = FALSE, x = cbind('d', 'n', 'beta_gamma_shape', 'sd_RE', 'Nm_lambda'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(paste0("3_analysis/helper/table_simulation_params_", uuid, ".txt"), append = TRUE, x = cbind(d, n, beta_gamma_shape, sd_RE, Nm_lambda), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

cat('Outfile: ', opt$outfile, '\n')
saveRDS(list(objects_counts=objects_counts, W=W, W_muts=W_muts,W_QP=W_QP, W_HilDA = hilda_obj, W_TCSM=W_muts_TCSM,
             d=d, n= n, beta_gamma_shape=beta_gamma_shape, sd_RE=sd_RE, lambda=lambda, Nm_lambda=Nm_lambda,
             X_sim = X_sim, beta = beta, Z_sim = Z_sim, u = u, lambdas = lambdas, alphabar = alphabar, alpha = alpha, Nm = rowSums(W), W = W),
        file = opt$outfile)


#----------------------------------------------------------------#