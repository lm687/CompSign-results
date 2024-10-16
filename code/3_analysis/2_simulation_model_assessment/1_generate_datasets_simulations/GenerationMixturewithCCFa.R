## Generating a dataset that can be used to test the differences between CompSign and TrackSig
## generate both a vcf file (for TrackSig) and an RDS matrix that CompSign can work with

## One single changepoint. This changepoint is always more or less in the same ccf quantile
# Signature extraction with mutSigExtractor

rm(list = ls())

library(optparse)
library(uuid)

debug <- F
if(debug){

  opt <- list(beta_gamma_shape=4)
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("../../../../code/") ## from where snakemake is run

}else{
  
  option_list = list(
    make_option(c("--input"), type="character", default=NA,
                help="Text with small description of the type of simulation being carried out", metavar="character"),
    make_option(c("--n"), type="numeric", default=NA,
                help="Number of samples", metavar="numeric"),
    make_option(c("--beta_gamma_shape"), type="numeric", default=NA,
                help="Fraction of additional mixture for differential abundance (from 0 to 1)", metavar="numeric"),
    make_option(c("--itnum"), type="numeric", default=NA,
                help="Index for replicate", metavar="numeric"),
    make_option(c("--lambda"), type="numeric", default=0,
                help="Overdispersion parameter", metavar="numeric"),
    make_option(c("--outfile"), type="character", default=NA,
                help="Output file in which to write the dataset (RDS file)", metavar="character"),
    ##### below: not used #####
    make_option(c("--nlambda"), type="numeric", default=NA,
                help="Parameter lambda for Poisson draws of number of mutations in sample", metavar="numeric"),
    make_option(c("--d"), type="numeric", default=NA,
                help="Number of features", metavar="numeric"),
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
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
source("2_inference_TMB/helper_TMB.R")
source("1_create_ROO/roo_functions.R")
source("1_create_ROO/helper_1_create_ROO.R")
source("../../../CompSign-do-not-modify_do_everything_online/R/helper_functions.R")
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
datasetgeneration <- 'GenerationMixturewithCCFa'
ct <- "Lymph-BNHL"
ccf_quantile = 0.6
dirichlet_precision_for_ccf_changepoint = opt$lambda #500
##----------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
##' Including all signatures, even with artefacts (SBS_SIGNATURE_PROFILES_V2 and SBS_SIGNATURE_PROFILES_V3 from
##' mutSigExtractor do not include artefactual signatures)
sigs_cosmic0 <- read.table(paste0( "../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
                           stringsAsFactors = FALSE, sep = ',', header = TRUE)
rownames(sigs_cosmic0) <- paste0(substr(sigs_cosmic0$SubType, 1, 1),'[',
                                 sigs_cosmic0$Type, ']', substr(sigs_cosmic0$SubType, 3, 3))
sigs_cosmic0 <- sigs_cosmic0[-c(1,2)];
sigs_cosmic <- colnames(sigs_cosmic0)
##-----------------------------------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Selecting data from one cancer type

cancer_type_trinucleotides <- load_PCAWG(ct, typedata = 'nucleotidesubstitution3', path_to_data = "../data/", override_warning_X_Z = T)
cancer_type_active_sigs <- load_PCAWG(ct, typedata = 'signaturesPCAWG', path_to_data = "../data/", override_warning_X_Z = T)
active_sigs <- colnames(cancer_type_active_sigs$Y)

n = nrow(cancer_type_trinucleotides$Y)/2
d = ncol(cancer_type_trinucleotides$Y)

cancer_type_trinucleotidesY = split_matrix_in_half(cancer_type_trinucleotides$Y)
samples <- unique(rownames(cancer_type_trinucleotides$Y))
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## use the ccfs from the real data
## but use the clonal/subclonal trinucleotides that we already have

# a <- read.table("../../../../../../do-not-use/GlobalDA_2023/restricted_data_from_GlobalDA/pcawg/repository_1567600367.tsv", sep = '\t')
muts_ccf <- sapply(samples, function(sample_it) read.table(paste0("../../../../DM/do-not-use/GlobalDA_2023/restricted_data_from_GlobalDA/pcawg/consensus_subclonal_reconstruction_mutccf_20170325/", sample_it, "_mutation_ccf.txt"), sep='\t', header = T), simplify = F)
muts_ccf <- sapply(muts_ccf, function(i) i[order(i$ccf, decreasing = T),], simplify = F) ## sort, even though this will be split into two groups according to the changepoint later
# vcfs <- sapply(samples, function(sample_it) read.table(paste0("../../../../../../do-not-use/GlobalDA_2023/restricted_data_from_GlobalDA/pcawg/pcawg_restricted_snv/", sample_it, ".consensus.20160830.somatic.snv_mnv.vcf.gz"), sep='\t', header = T), simplify = F)
# vcfs$`00b9d0e6-69dc-4345-bffd-ce32880c8eef`

## mutations are as found in the dataset, with the same ccf
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## select one ccf changepoint. they are all very tightly around ccf_quantile (non-dispersed Dirichlet)
ccf_quantile_per_sample <- MCMCpack::rdirichlet(n = length(muts_ccf), alpha = c(ccf_quantile, 1-ccf_quantile)*dirichlet_precision_for_ccf_changepoint)[,1]
ccf_quantile_per_sample

ccf_changepoints_per_sample <- sapply(1:length(samples), function(samples_it)  quantile(muts_ccf[[samples_it]]$ccf, na.rm = T,
                                                         probs = ccf_quantile_per_sample[[samples_it]]) )

## We select the mutations before and after the changepoint, as we are going to draw from those
mutations_per_group = sapply(1:length(samples), function(samples_ct){
  x = list(muts_ccf[[samples_ct]][muts_ccf[[samples_ct]]$ccf < ccf_changepoints_per_sample[samples_ct],][,c(1,2,4)],
                                                                           muts_ccf[[samples_ct]][muts_ccf[[samples_ct]]$ccf >= ccf_changepoints_per_sample[samples_ct],][,c(1,2,4)])
  return(sapply(x, function(x2) x2[!is.na(x2$ccf),], simplify = F)) ## remove mutations which are NA
  }, simplify = F)
names(mutations_per_group) <- samples


## before the changepoint, all mutations are drawn from the pool of mutations of the first group

## after the changepoint, the ccfs from the cancer type are kept, but the mutations are sampled (with probability:
pi_mixtures = c(1/exp(opt$beta_gamma_shape), 1-1/exp(opt$beta_gamma_shape))
## from the first group and from the second group respectively. The higher beta_gamma_shape, the higher pi_mixtures[1],
## and the clearer the DA. At beta_gamma_shape=0, all the mixture is from the first group (no DA)

mutations_per_group = sapply(samples, function(it_sample){
  sapply(1:length(mutations_per_group[[it_sample]]), function(group){
  if(group == 1){
    muts = sample(x = names(cancer_type_trinucleotidesY[[group]][it_sample,]),
           prob = normalise_rw(cancer_type_trinucleotidesY[[group]][it_sample,]),
           size = nrow(mutations_per_group[[it_sample]][[group]]), replace = T)
    ## all mutations come from the same group (the first one)
  }else if(group == 2){
    ## mutations come from a mixture
    muts_per_mixture = round(pi_mixtures*nrow(mutations_per_group[[it_sample]][[group]]))
    stopifnot(sum(muts_per_mixture) == nrow(mutations_per_group[[it_sample]][[group]]))
    muts1 = sample(x = names(cancer_type_trinucleotidesY[[1]][it_sample,]),
                  prob = normalise_rw(cancer_type_trinucleotidesY[[1]][it_sample,]),
                  size = muts_per_mixture[1], replace = T)
    muts2 = sample(x = names(cancer_type_trinucleotidesY[[group]][it_sample,]),
                   prob = normalise_rw(cancer_type_trinucleotidesY[[group]][it_sample,]),
                   size = muts_per_mixture[2], replace = T)
    muts <- sample(c(muts1, muts2), replace = F, size = nrow(mutations_per_group[[it_sample]][[group]])) ## shuffle
  }else{
    stop
  }

  cbind(mutations_per_group[[it_sample]][[group]], mutation=muts)
  }, simplify=F)
}, simplify=F)
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## object for running run_TMB

# Split the matrices below into 100-mutation bins
binSize=100
mutations_per_group_bins100 <- sapply(names(mutations_per_group), function(sample_it){
  f_mat = do.call('rbind', mutations_per_group[[sample_it]])
  f_to_split = rep(1:(1+(round(nrow(f_mat)/binSize))), binSize) ## the bin they will belong to
  res = split(f_mat, f=f_to_split[1:nrow(f_mat)]) ## the last bin might be quite empty, <100 muts.
  return(res)
})
  
give96mat <- function(mutations_per_group_arg, equal_groups_per_patient=T){
  if(equal_groups_per_patient){
    ## in cases where there are 2 groups per patient, for instance
    nlevels = length(mutations_per_group_arg[[1]]) ## 2 if two groups, as many as bins otherwise
    W96 = do.call('rbind', sapply(1:nlevels, function(group_it){do.call('rbind', t(sapply(names(mutations_per_group_arg), function(sample_it){
      table(factor(mutations_per_group_arg[[sample_it]][[group_it]]$mutation, levels=rownames(sigs_cosmic0)))
    }, simplify = F)))}, simplify = F))
    rownames(W96) <- rep(names(mutations_per_group_arg), nlevels)
  }else{
    ## when there are a different number of bins per patient
    nlevels_per_patient = sapply(mutations_per_group_arg, length)
    W96 = sapply(names(mutations_per_group_arg), function(sample_it){
      sapply(1:nlevels_per_patient[[sample_it]], function(group_it){
        table(factor(mutations_per_group_arg[[sample_it]][[group_it]]$mutation, levels=rownames(sigs_cosmic0)))},
        simplify = F)}, simplify = F)
    ## right now it's nested with patients -> bins.
    W96 <- do.call('rbind', sapply(W96, function(i) do.call('rbind', i)))
    rownames(W96) <- make.unique(rep(names(mutations_per_group_arg), sapply(mutations_per_group_arg, length)))
    warning("Note that the output of this function is different than if there are two groups")
  }
  return(W96)
}

W96 <- give96mat(mutations_per_group)
W96bins <- give96mat(mutations_per_group_bins100, equal_groups_per_patient = F)


W96_split = split_matrix_in_half(W96)
objects_counts_in_list = list(Y=do.call('rbind', W96_split),
                              x=cbind(1, rep(c(0,1), sapply(W96_split, nrow))),
                              z=give_z_matrix(n_times_2 = length(names(mutations_per_group))*2))


## extract signatures for each of the two groups
W_twogroup_QP_list <- extract_sigs_TMB_obj(dataset_obj_trinucleotide = objects_counts_in_list, subset_signatures = active_sigs,
                          signature_definition = sigs_cosmic0, signature_fitting_method = 'QP')

objects_counts <- new("exposures_cancertype",
                      cancer_type=paste0('Simulated_from_', ct),
                      type_classification = "simulated two group",
                      number_categories = 2,
                      id_categories = c('CCFhigh', 'CCFlow'),
                      active_signatures = "absent; simulation",
                      count_matrices_all = split_matrix_in_half(W_twogroup_QP_list$Y),
                      count_matrices_active = split_matrix_in_half(W_twogroup_QP_list$Y),
                      sample_names = names(mutations_per_group),
                      modification = "none",
                      is_null_active = TRUE,
                      is_empty="Non-empty")

Z_matrix_one_obs_per_patient <- diag(length(mutations_per_group_bins100))
Z_matrix_bins <- do.call('rbind', sapply(1:length(mutations_per_group_bins100), function(i) Z_matrix_one_obs_per_patient[rep(i, length(mutations_per_group_bins100[[i]])),]))

objects_counts_in_list_bins = list(Y=W96bins,
                              x=cbind(baseline=1, ## baseline
                                      bin_order=unlist(sapply(sapply(mutations_per_group_bins100, length), function(j) 1:j)), ## bin (order)
                                      max_ccf=unlist(sapply(mutations_per_group_bins100, function(i) sapply(i, function(j) max(j$ccf)))),## max ccf of each bin
                                      min_ccf=unlist(sapply(mutations_per_group_bins100, function(i) sapply(i, function(j) min(j$ccf)))),## min ccf of each bin
                                      mean_ccf=unlist(sapply(mutations_per_group_bins100, function(i) sapply(i, function(j) mean(j$ccf)))),## mean ccf of each bin
                                      median_ccf=unlist(sapply(mutations_per_group_bins100, function(i) sapply(i, function(j) median(j$ccf))))## median ccf of each bin
                                             ),
                              z=Z_matrix_bins)

## extract signatures for the bins, too
W_bins_QP_list <- extract_sigs_TMB_obj(dataset_obj_trinucleotide = objects_counts_in_list_bins,
                                       subset_signatures = active_sigs,
                                       signature_definition = sigs_cosmic0, signature_fitting_method = 'QP')

##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## save: 
## 1. clonal and subclonal exposures
## 2. exposures for bins of 100 mutations
## 3. the original dataset with ccf information
cat('Outfile: ', opt$outfile, '\n')

to_save = list(objects_counts=objects_counts, 
               exposures_twogroup_QP_list=W_twogroup_QP_list,
               exposures_bins_QP_list=W_bins_QP_list,
               simulated_mutations=mutations_per_group,
               objects_counts_in_list=objects_counts_in_list,
               objects_counts_in_list_bins=objects_counts_in_list_bins,
               W_bins_QP_list='see W_bins',
               W_QP=NULL,
               W_HilDA = NULL, W_TCSM=NULL,
               mutations_per_group_ccf = 'see <W_muts>',
               d=d, n= n, beta_gamma_shape=opt$beta_gamma_shape, sd_RE=NULL, lambda=opt$lambda, Nm_lambda=NULL,
               X_sim = NULL, beta = NULL, Z_sim = NULL, u = NULL, lambdas = NULL, alphabar = NULL, alpha = NULL,
               Nm = NULL)
        
system(paste0('mkdir -p ', dirname(opt$outfile)))
# output = paste0("../data/assessing_models_simulation/multiple_", datasetgeneration, "_", n, "_NA_", dirichlet_precision_for_ccf_changepoint, "_", d,
# "_", opt$beta_gamma_shape, "_NA_NA_NA/multiple_", datasetgeneration, "_", n, "_NA_", dirichlet_precision_for_ccf_changepoint, "_", d, "_", opt$beta_gamma_shape, "_NA_NA_NA_dataset", opt$itnum, ".RDS")
saveRDS(to_save, opt$outfile)
##----------------------------------------------------------------------------##

