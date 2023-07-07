rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggrepel)
library(Ternary)
library(MCMCpack)
library(dplyr)
library(gridExtra)
library(parallel)
library(TMB)


TMB::compile("../../2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda.cpp", "-std=gnu++17")
dyn.load(dynlib("../../2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda"))
TMB::compile("../../2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda2.cpp", "-std=gnu++17")
dyn.load(dynlib("../../2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda2"))
TMB::compile("../../2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda.cpp", "-std=gnu++17")
dyn.load(dynlib("../../2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda"))
TMB::compile("../../2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial"))
TMB::compile("../../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial"))
TMB::compile("../../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomialsinglelambda.cpp", "-std=gnu++17")
dyn.load(dynlib("../../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomialsinglelambda"))
TMB::compile("../../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomialsinglelambda2.cpp", "-std=gnu++17")
dyn.load(dynlib("../../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomialsinglelambda2"))
TMB::compile("../../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial_singlecov.cpp", "-std=gnu++17")
dyn.load(dynlib("../../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial_singlecov"))
TMB::compile("../../2_inference_TMB/mm_multinomial/fullRE_ME_halfdirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../2_inference_TMB/mm_multinomial/fullRE_ME_halfdirichletmultinomial"))
TMB::compile("../../2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial"))
TMB::compile("../../2_inference_TMB/mm_multinomial/fullRE_ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../2_inference_TMB/mm_multinomial/fullRE_ME_multinomial"))


source("../../2_inference_TMB/helper_TMB.R")
source("../../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")

enough_samples = read.table("../../../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
ct <- enough_samples[5]
ct
nonexogenous = read.table("../../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t", comment.char = "#", fill = F)

# rm(list = c('fullRE_DM', 'fullRE_DM_ns1', 'fullRE_DM_nonexo', 'diagRE_DM', 'sparseRE_DM'))
# fullRE_M <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/fullRE_M_", ct, "_signatures.RDS"))
# fullRE_M
# fullRE_DM <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_DMSL_", ct, "_signatures.RDS"))
# fullRE_DM
# fullRE_DM_ns1 <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_DMSL_", ct, "_nucleotidesubstitution1.RDS"))
# fullRE_DM_ns1
# fullRE_DM_nonexo_old <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_DM_", ct, "_signatures.RDS"))
# fullRE_DM_nonexo_old
# fullRE_DM_nonexo <- readRDS(paste0("../../data/../pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_DMSL_", ct, "_signatures.RDS"))
# fullRE_DM_nonexo
# sparseRE_DM_nonexo <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_nonexo_DMSL_", ct, "_signatures.RDS"))
# sparseRE_DM_nonexo
# diagRE_DM <- readRDS(paste0("../../data/../pcawg_robjects_cache/tmb_results/nlminb/diagRE_DMSL_", ct, "_signatures.RDS"))
# diagRE_DM
# sparseRE_DM <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_DMSL2_", ct, "_signatures.RDS"))
# sparseRE_DM


# wald_TMB_wrapper(fullRE_M, verbatim = F)
# wald_TMB_wrapper(fullRE_DM_ns1, verbatim = F)
# wald_TMB_wrapper(diagRE_DM)
# wald_TMB_wrapper(sparseRE_DM)
# wald_TMB_wrapper(fullRE_DM_nonexo)


diagRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/diagRE_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

diagRE_DMDL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/diagRE_DM_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMDL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/fullRE_DM_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_M_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_M_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_M <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/fullRE_M_", ct, "_signatures.RDS")))
}, simplify = F)

diagRE_M <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/diagRE_M_", ct, "_signatures.RDS")))
}, simplify = F)

sparseRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_DMSL2_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMDL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fulLRE_nonexo_DM_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMDL_sortednonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fulLRE_sortednonexo_DM_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMSL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

diagRE_DMSL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/diagRE_nonexo_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMSL_SBS1 <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_SBS1baseline_DMSL_", ct, "signatures.RDS")))
}, simplify = F)

fullRE_halfDM <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_halfDM_", ct, "signatures.RDS")))
}, simplify = F)

fullRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

diagRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/diagRE_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

sparseRE_DMSL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_nonexo_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMSL_nonexo[sapply(sparseRE_DMSL_nonexo, typeof) == "character"]

pvals_fullRE_M <- sapply(fullRE_M, function(i) try(wald_TMB_wrapper(i)))
pvals_fullRE_M <- p.adjust(pvals_fullRE_M)
pvals_diagRE_DM <- sapply(diagRE_DMSL, function(i) try(wald_TMB_wrapper(i)))
pvals_diagRE_DM <- p.adjust(pvals_diagRE_DM)
pvals_DM <- sapply(sparseRE_DMSL, function(i) try(wald_TMB_wrapper(i)))
pvals_DM <- p.adjust(pvals_DM)
pvals_DMnonexo <- sapply(sparseRE_DMSL_nonexo, function(i) try(wald_TMB_wrapper(i)))
pvals_DMnonexo <- p.adjust(pvals_DMnonexo)

sort(pvals_DM)
sort(pvals_DMnonexo)

num_samples <- sapply(enough_samples, function(ct){
    .xx <- try(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures"),
                                       sigs_to_remove = unique(nonexogenous$V1)))
    try(nrow(.xx$Y)/2)
})

effect_size <- sapply(enough_samples, function(ct){
  .xx <- try(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures"),
                                     sigs_to_remove = unique(nonexogenous$V1)))
  try(sum((normalise_rw(.xx$Y[1:(nrow(.xx$Y)/2),]) - normalise_rw(.xx$Y[(1+nrow(.xx$Y)/2):(nrow(.xx$Y)),]))**2)/(nrow(.xx$Y)/2))
})

length(python_like_select_name(fullRE_DM_nonexo$par.fixed, grep_substring = "beta"))
length(python_like_select_name(fullRE_DM_nonexo_old$par.fixed, grep_substring = "beta"))
length(python_like_select_name(sparseRE_DM_nonexo$par.fixed, grep_substring = "beta"))


ggplot(cbind.data.frame(pvals_DM=pvals_DM, pvals_DMnonexo=pvals_DMnonexo,
                        num_samples=as.numeric(num_samples),
                        ct=enough_samples),
       aes(x=log(pvals_DM), y=log(pvals_DMnonexo), col=num_samples,
           label=ct))+geom_point()+
  geom_hline(yintercept = log(0.05))+geom_vline(xintercept = log(0.05))+
  geom_label_repel()

ggplot(cbind.data.frame(pvals_DM=pvals_DM, pvals_DMnonexo=pvals_DMnonexo,
                        num_samples=as.numeric(num_samples),
                        ct=enough_samples),
       aes(x=(pvals_DM), y=(pvals_DMnonexo), col=num_samples,
           label=ct))+geom_point()+
  geom_hline(yintercept = (0.05))+geom_vline(xintercept = (0.05))+
  geom_label_repel()

plot(pvals_fullRE_M, pvals_DM)

ggplot(cbind.data.frame(effect_size=as.numeric(effect_size), minlogpval=-log(as.numeric(pvals_DM)),
                        label=enough_samples), aes(x=effect_size, y=minlogpval, label=label))+
  geom_point()+geom_label_repel(alpha=0.2)


ct <- "Breast-AdenoCA"

## Only for some CT is the cov matrix positive semi-definite!

ct <- enough_samples[5]
ct <- "Prost-AdenoCA"
subset_sigs_sparse_cov_idx_nonexo <- read.table("../../../current/subset_sigs_sparse_cov_idx_nonexo.txt", stringsAsFactors = F, fill = T)
sparseRE_DMSL_nonexo$`Breast-AdenoCA`



ct <- enough_samples[10]

pca_sim <- give_sim_from_estimates(ct, "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T)

give_sim_from_estimates(enough_samples[2], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T)
give_sim_from_estimates(enough_samples[3], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T, sig_of_interest = 'SBS2')
give_sim_from_estimates(enough_samples[4], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T, sig_of_interest = 'SBS8')
give_sim_from_estimates(enough_samples[5], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T, sig_of_interest = 'SBS8')
# give_sim_from_estimates(enough_samples[5], "signatures", sigs_to_remove=unique(nonexogenous$V1), model = "fullRE_M", bool_give_PCA = T, sig_of_interest = 'SBS8')

ggsave(paste0("../../results/results_TMB/pcawg/simulation_from_estimates/", ct, ".pdf"))
# counts = apply(alpha, 1, function(i) 

give_sim_from_estimates(enough_samples[5], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T, sig_of_interest = 'SBS8')

CNS_Medullo <- give_subset_sigs_TMBobj(load_PCAWG(enough_samples[5], "signatures", path_to_data = "../../../data/"),
                                sigs_to_remove = unique(nonexogenous$V1))

give_barplot_from_obj(CNS_Medullo, T)

# mvtnorm::rmvnorm(n = 100, mean = rep(0,dmin1), sigma = cov_matb)



amalgamation1 <- as(normalise_rw(all_probs[,10:12]), 'matrix')
rownames(amalgamation1)[rownames(amalgamation1) == ""] = paste0('Sim', 1:n_sim)
bad_idx <- as.numeric(which( (is.na(rowSums(amalgamation1))) | !(rowSums(amalgamation1) == 1)))
amalgamation1 <- amalgamation1[-bad_idx,]
df_pca_amalgamation1 <- df_pca[-bad_idx,]



aa <- sapply(c('Observed', 'Simulated'), function(arg_group){
  am <- amalgamation1[df_pca$col[-bad_idx] == arg_group,]
  dens <- TernaryDensity(am, resolution = 10L)
  plot.new()
  TernaryPlot(atip = colnames(am)[1], btip = colnames(am)[2], ctip = colnames(am)[3],
              grid.lines = 0, grid.col = NULL)
  ColourTernary(dens, direction = 1)
  # ColourTernary(amalgamation, spectrum=rainbow(256))
  TernaryPoints(am, col = factor(df_pca$col[-bad_idx]), pch = '.', cex=5)
  TernaryDensityContour(am, resolution = 30L)
})  

aa[1]
aa[2]

eigen(cov_mat)$values
eigen(cov_matb)$values

#obj_nonexo$z %*% 

table(sapply(fullRE_M_nonexo, function(i) try(i$pdHess)))
table(sapply(fullRE_M, function(i) try(i$pdHess)))
table(sapply(diagRE_M, function(i) try(i$pdHess)))
table(sapply(diagRE_DMSL, function(i) try(i$pdHess))) ## very good
table(sapply(fullRE_DMSL_nonexo, function(i) try(i$pdHess))) ## half-half
table(sapply(diagRE_DMSL_nonexo, function(i) try(i$pdHess))) ## very good

fullRE_DM_Prost_AdenoCA <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_DMSL_", "Prost-AdenoCA", "_signatures.RDS"))
fullRE_DM_Prost_AdenoCA$pdHess
fullRE_DM_Prost_AdenoCA_nonexo <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_DM_Prost-AdenoCA_signatures.RDS"))
fullRE_DM_Prost_AdenoCA_nonexo$pdHess
fullRE_M_Prost_AdenoCA_nonexo <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_M_Prost-AdenoCA_signatures.RDS"))
fullRE_M_Prost_AdenoCA_nonexo$pdHess

give_sim_from_estimates(ct = "Prost-AdenoCA", typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS8', model = "fullRE_M")
ct = "Prost-AdenoCA"
typedata = "signatures"
bool_nonexo = T
sigs_to_remove=unique(nonexogenous$V1)
bool_give_PCA = T
sig_of_interest = 'SBS8'
model = "fullRE_M"

fullRE_DM_Prost_AdenoCA

give_sim_from_estimates(ct = enough_samples[2], typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_M")

give_sim_from_estimates(ct = enough_samples[9], typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_M")


give_sim_from_estimates(ct = enough_samples[10], typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_M")

give_sim_from_estimates(ct = "Lung-SCC", typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_DM")

give_sim_from_estimates(ct = "Lung-SCC", typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_M")


fullRE_DMSL_nonexo$`Lymph-CLL`$pdHess
fullRE_M_nonexo$`Lymph-CLL`$pdHess

Lymph_CLL_obj <- give_subset_sigs_TMBobj(load_PCAWG(ct = "Lymph-CLL", typedata = "signatures", path_to_data = "../../../data/"),
                        sigs_to_remove = unique(nonexogenous$V1))
give_ranked_plot_simulation(tmb_fit_object = fullRE_DMSL_nonexo$`Lymph-CLL`,
                            data_object = Lymph_CLL_obj, print_plot = T, nreps = 20, model = "DMSL", integer_overdispersion_param = 1)

good_plots = give_ranked_plot_simulation(tmb_fit_object = fullRE_DMSL_nonexo$`Lymph-CLL`,
                                         data_object = Lymph_CLL_obj, print_plot = F, nreps = 20, model = "DMSL", integer_overdispersion_param = 1)
good_plots_unsorted = lapply(list(good_plots), function(i) lapply(i, function(j) cbind.data.frame(sorted_value=as.vector(j), rank_number=1:length(j)) ))
give_interval_plots_2(df_rank = good_plots_unsorted[[1]], data_object = Lymph_CLL_obj,
                      loglog = F, title = 'Lymph_CLL')

## here there's a proble!! we're using DMSL with the M model
good_plots_M = give_ranked_plot_simulation(tmb_fit_object = fullRE_DMSL_nonexo$`Lymph-CLL`,
                                         data_object = Lymph_CLL_obj, print_plot = F, nreps = 20, model = "M")
good_plots_unsorted_M = lapply(list(good_plots_M), function(i) lapply(i, function(j) cbind.data.frame(sorted_value=as.vector(j), rank_number=1:length(j)) ))
give_interval_plots_2(df_rank = good_plots_unsorted_M[[1]], data_object = Lymph_CLL_obj,
                      loglog = F, title = 'Lymph_CLL')

pdf("../../../results/results_TMB/pcawg/assessing_models/ranked_plots/Lymph_CLL_DMSL_vs_M.pdf",
    height = 3.5, width = 5)
grid.arrange(give_interval_plots_2(df_rank = good_plots_unsorted_M[[1]], data_object = Lymph_CLL_obj,
                                   loglog = F, title = 'Lymph CLL (M)')+guides(col=guide_legend(title="In confidence interval")),
             give_interval_plots_2(df_rank = good_plots_unsorted[[1]], data_object = Lymph_CLL_obj,
                                   loglog = F, title = 'Lymph CLL (DM)')+guides(col=guide_legend(title="In confidence interval")),
             nrow=1
)
dev.off()

# tmb_fit_object = fullRE_DMSL_nonexo$`Lymph-CLL`
# data_object = Lymph_CLL_obj
# print_plot = T
# nreps = 20
# model = "DMSL"
# integer_overdispersion_param = 1
# 
# unction(tmb_fit_object, data_object, print_plot = T, nreps = 1, model, integer_overdispersion_param){
#   
#   if(model == 'M'){
#     ## theta is always going to be the same. Only replicate the draws from the multinomial
#     
#     sim_theta = simulate_from_M_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
#                                     x_matrix = data_object$x, z_matrix = data_object$z)
#     
#     sim_counts = t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) )
#     
#     if(nreps>1){
#       sim_counts = replicate(nreps, t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) ), simplify = F)
#     }
#   }else if(model %in% c('DM', 'DMSL')){
#     give_sim_data = function(){
#       if(model == 'DM'){
#         sim_theta = simulate_from_DM_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
#                                          x_matrix = data_object$x, z_matrix = data_object$z, integer_overdispersion_param=integer_overdispersion_param)
#       }else if(model == 'DMSL'){
#         sim_theta = simulate_from_DMSL_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
#                                            x_matrix = data_object$x, z_matrix = data_object$z, integer_overdispersion_param=integer_overdispersion_param)
#       }
#       sim_counts = t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) )
#       return(sim_counts)
#     }    
#     if(nreps>1){
#       sim_counts = replicate(nreps, give_sim_data(), simplify = F)
#     }else{
#       sim_counts = give_sim_data()
#     }
#   }else{
#     stop('Specify a correct model')
#   }
#   
#   stopifnot(all(dim(data_object$Y) == dim(sim_counts)))
#   if(print_plot)  plot(sort(data_object$Y), sort(sim_counts))
#   
#   return(sim_counts)
#   
# }
# 
# 
# tmb_fit_object = tmb_fit_object
# full_RE = T
# x_matrix = data_object$x
# z_matrix = data_object$z
# integer_overdispersion_param=integer_overdispersion_param
# 
# function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, integer_overdispersion_param){
#   dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
#   overdispersion_lambda = rep(integer_overdispersion_param*exp(python_like_select_name(tmb_fit_object$par.fixed, "log_lambda")), nrow(x_matrix))
#   if(full_RE){
#     re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
#     ntimes2 = nrow(z_matrix)
#     logRmat = z_matrix %*% re_mat + 
#       x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
#     sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* softmax(c(logRmat[l,], 0)))))
#   }else{
#     sim_thetas = softmax(cbind(sapply(1:dmin1,
#                                       function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
#                                  give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
#     sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* sim_thetas[l,])))
#   }
#   return(sim_thetas)
# }



ct <- enough_samples[1]

for(ct in enough_samples){
  pdf(paste0("../../../results/results_TMB/pcawg/summaries_betas/", gsub("[.]", "_", ct), "_allsigs.pdf"), height = 2, width = 9)
  grid.arrange(plot_betas(fullRE_M[[ct]])+ggtitle(paste0(ct, '\n fullRE_M')),
  plot_betas(diagRE_M[[ct]])+ggtitle(paste0(ct, '\n diagRE_M')),
  plot_betas(fullRE_DMSL[[ct]])+ggtitle(paste0(ct, '\n fullRE_DMSL')),
  plot_betas(diagRE_DMSL[[ct]])+ggtitle(paste0(ct, '\n diagRE_DMSL')),
  plot_betas(sparseRE_DMSL[[ct]])+ggtitle(paste0(ct, '\n sparseRE_DMSL')), nrow=1)
  dev.off()
  
  pdf(paste0("../../../results/results_TMB/pcawg/summaries_betas/", gsub("[.]", "_", ct), "_nonexo_allsigs.pdf"), height = 2, width = 9)
  grid.arrange(plot_betas(fullRE_M_nonexo[[ct]])+ggtitle(paste0(ct, '\n fullRE_M_nonexo')),
  plot_betas(fullRE_DMSL_nonexo[[ct]])+ggtitle(paste0(ct, '\n fullRE_DMSL_nonexo')),
  plot_betas(diagRE_DMSL_nonexo[[ct]])+ggtitle(paste0(ct, '\n diagRE_DMSL_nonexo')),
  plot_betas(sparseRE_DMSL_nonexo[[ct]])+ggtitle(paste0(ct, '\n sparseRE_DMSL_nonexo')),
                                                 nrow=1)
  dev.off()
}

for(ct in gsub("[.]", "_", enough_samples)){
  cat("\\begin{figure}[h]")
  cat("\\includegraphics[width=\\textwidth]{/Users/morril01/Documents/PhD/GlobalDA/results/results_TMB/pcawg/summaries_betas/", ct, "_allsigs.pdf}\n", sep = "")
  cat("\\includegraphics[width=\\textwidth]{/Users/morril01/Documents/PhD/GlobalDA/results/results_TMB/pcawg/summaries_betas/", ct, "_nonexo_allsigs.pdf}", sep = "")
  cat("\\caption{", gsub("_", " ", ct), "}", sep = "")
  cat("\\end{figure}\n\n\n")
}


fullRE_M[[ct]]
diagRE_M[[ct]]
fullRE_DMSL[[ct]]
diagRE_DMSL[[ct]]
sparseRE_DMSL[[ct]]

fullRE_M_nonexo[[ct]]
fullRE_DMSL_nonexo[[ct]]
diagRE_DMSL_nonexo[[ct]]
sparseRE_DMSL_nonexo[[ct]]

.dm1 <- wrapper_run_TMB(model = "fullREDMsinglelambda",
                object = load_PCAWG(ct = ct, typedata = "signatures", path_to_data = "../../../data/"))
.dm1
.dm1_sorted <- wrapper_run_TMB(model = "fullREDMsinglelambda",
                        object = sort_columns_TMB(load_PCAWG(ct = ct, typedata = "signatures", path_to_data = "../../../data/")))
.dm1_sorted
give_subset_sigs_TMBobj(load_PCAWG(ct = i[1,1], typedata = i[1,2]),
                        sigs_to_remove = unique(nonexogenous$V1))
.obj <- sort_columns_TMB(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures", path_to_data = "../../../data/"),
                                         sigs_to_remove = unique(nonexogenous$V1)))
.dm1_nonexo <- wrapper_run_TMB(model = "fullREDMsinglelambda",
                        object = .obj, use_nlminb = T)
.dm1_nonexo

ct <- enough_samples[2]
give_barplot_from_obj(load_PCAWG(ct = ct, typedata = "signatures",
                                 path_to_data = "../../../data/"), legend_on = T)
.dm1_DMDL <- wrapper_run_TMB(model = "fullRE_DM",
                        object = sort_columns_TMB(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures",
                                            path_to_data = "../../../data/"),
                                            sigs_to_remove = c('NANANANA'))),
                        use_nlminb = T, smart_init_vals = T)
.dm1_DMDL
fill_covariance_matrix(arg_d = length(python_like_select_name(.dm1_DMDL$par.fixed, 'beta'))/2,
                       arg_entries_var = python_like_select_name(.dm1_DMDL$par.fixed, 'logs_sd_RE'),
                       arg_entries_cov = python_like_select_name(.dm1_DMDL$par.fixed, 'cov_par_RE'))
saveRDS(object = .dm1_DMDL, file=paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/", "fullRE_DMDL_sorted20210524_",
                                paste0(ct, "signatures", collapse = "_"), ".RDS"))

plot_betas(.dm1_nonexo, names_cats = paste0(colnames(.obj$Y)[-ncol(.obj$Y)], '/', colnames(.obj$Y)[ncol(.obj$Y)]))+labs(x='')


mclapply(enough_samples, function(ct){
             x = wrapper_run_TMB(model = "fullREDMsinglelambda", object=sort_columns_TMB_SBS1(load_PCAWG(ct = ct, typedata = "signatures", path_to_data = "../../../data/")),
                                 use_nlminb = T)
             x
             saveRDS(object = x, file=paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/", "fullRE_SBS1baseline_DMSL_",
                                             paste0(ct, "signatures", collapse = "_"), ".RDS"))
})

fullRE_DMSL_SBS1_betas <- lapply(fullRE_DMSL_SBS1, function(i){
  .x <- try(give_betas(i)[2,])
  if((typeof(.x) == 'character')){
    .x <- NA
  }else{
    .sum_i = summary(i)
    .x <- t(python_like_select_rownames(.sum_i, 'beta')[c(F,T),])
  }
  .x
})
for(i in 1:length(fullRE_DMSL_SBS1)){
  if(!(typeof(fullRE_DMSL_SBS1[[i]]) == 'character')){
    .nmes <- colnames(sort_columns_TMB_SBS1(load_PCAWG(ct = enough_samples[i], typedata = "signatures", path_to_data = "../../../data/"))$Y)
    colnames(fullRE_DMSL_SBS1_betas[[i]]) = paste0(.nmes[-length(.nmes)], '/', .nmes[length(.nmes)])
      
  }
}

fullRE_DMSL_SBS1_betas_all <- lapply(1:length(fullRE_DMSL_SBS1_betas), function(i) try(data.frame(ct=names(fullRE_DMSL_SBS1_betas[i]), beta=t(fullRE_DMSL_SBS1_betas[[i]]),
                                                                                                               logR=colnames(fullRE_DMSL_SBS1_betas[[i]]))))
fullRE_DMSL_SBS1_betas_all <- do.call('rbind', fullRE_DMSL_SBS1_betas_all[sapply(fullRE_DMSL_SBS1_betas_all, typeof) == 'list'])

fullRE_DMSL_SBS1_betas_all[!grepl("/SBS1$", fullRE_DMSL_SBS1_betas_all$logR),]
## select only those with SBS1 as baseline
fullRE_DMSL_SBS1_betas_all <- fullRE_DMSL_SBS1_betas_all[grepl("/SBS1$", fullRE_DMSL_SBS1_betas_all$logR),]

fullRE_DMSL_SBS1_betas_all$phHess <- sapply(fullRE_DMSL_SBS1, function(i) try(i$pdHess))[match(fullRE_DMSL_SBS1_betas_all$ct, names(fullRE_DMSL_SBS1))]

ggplot(fullRE_DMSL_SBS1_betas_all, aes(x=ct, col=logR, y=beta.Estimate))+geom_point()+
  facet_wrap(.~logR, scales = "free_x", nrow=5)+
  geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std..Error`, ymax=`beta.Estimate`+`beta.Std..Error`), width=.1)+
  geom_hline(yintercept = 0, col='blue', lty='dashed')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
ggsave("../../../results/results_TMB/pcawg/all_betas_all_ct.pdf", width = 25, height = 15)

multiple_obs_SBS1_betas <- fullRE_DMSL_SBS1_betas_all %>% dplyr::select(logR) %>% table > 2
ggplot(fullRE_DMSL_SBS1_betas_all %>% filter(logR %in% names(multiple_obs_SBS1_betas[multiple_obs_SBS1_betas])), aes(x=ct, col=logR, y=beta.Estimate))+geom_point()+
  facet_wrap(.~logR, scales = "free_x", nrow=5)+
  geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std..Error`, ymax=`beta.Estimate`+`beta.Std..Error`), width=.1)+
  geom_hline(yintercept = 0, col='blue', lty='dashed')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))


fullRE_DMSL_SBS1_betas_all$exogenous = (gsub("/.*", "", fullRE_DMSL_SBS1_betas_all$logR) %in% nonexogenous$V1)
ggplot(fullRE_DMSL_SBS1_betas_all, aes(x=logR, y=beta.Estimate, col=exogenous))+geom_point()+
  facet_wrap(.~ct, scales = "free_x", nrow=5)+
  geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std..Error`, ymax=`beta.Estimate`+`beta.Std..Error`), width=.1)+
  geom_hline(yintercept = 0, col='blue', lty='dashed')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
ggsave("../../../results/results_TMB/pcawg/all_betas_all_ct_byct.pdf", width = 12, height = 12)

#-------------------------------------------------------------------------------------#

mclapply(enough_samples, function(ct){
  x = wrapper_run_TMB(model = "fullREhalfDM", object=sort_columns_TMB(load_PCAWG(ct = ct, typedata = "signatures", path_to_data = "../../../data/")),
                      use_nlminb = T)
  x
  saveRDS(object = x, file=paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/", "fullRE_halfDM_",
                                  paste0(ct, "signatures", collapse = "_"), ".RDS"))
})
sapply(fullRE_halfDM, function(i) try(i$pdHess))

#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#

list_models <- c( 'diagRE_M', 'fullRE_M',
                  'diagRE_DMDL','fullRE_halfDM', 'fullRE_DMDL', 
                  'diagRE_DMSL','sparseRE_DMSL', 'fullRE_DMSL', 'fullRE_DMSL_SBS1',
                  'fullRE_M_nonexo','diagRE_DMSL_nonexo','sparseRE_DMSL_nonexo', 'fullRE_DMSL_nonexo',
                  'fullRE_DMDL_nonexo', 'fullRE_DMDL_sortednonexo')

all_summaries <- lapply(lapply(list_models, get), function(i){
    give_summary_of_runs2(i, long_return = T)})
names(all_summaries) <- list_models

ggplot(melt(all_summaries), aes(x=factor(L1, levels=list_models), y=value, fill=L2))+geom_tile()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))


# load_PCAWG(ct = "Skin-Melanoma.mucosal", typedata = "signatures", path_to_data = "../../../data/")
  
ct <- "Bone-Osteosarc"
obj <- sort_columns_TMB(load_PCAWG(ct = ct, typedata = "signatures",
                                   path_to_data = "../../../data/"))
sortedM <- wrapper_run_TMB(model = "fullRE_M",
                           object = obj)
sortedM

## use params from ME M for ME DM
dmin1 <- ncol(obj$Y)-1

sortedDM <- wrapper_run_TMB(model = "fullRE_DM",
                           object = obj,
                           smart_init_vals = F, use_nlminb = T,
                           initial_params = list(
                             beta = matrix(python_like_select_name(sortedM$par.fixed, 'beta'), nrow=2),
                             u_large = matrix(sortedM$par.random, ncol=dmin1),
                             logs_sd_RE=python_like_select_name(sortedM$par.fixed, 'logs_sd_RE'),
                             cov_par_RE = python_like_select_name(sortedM$par.fixed, 'cov_par_RE'),
                             log_lambda = matrix(c(2,2))))
sortedDM

non_duplicated_rows <- function(i){
  rownames(i)[duplicated(rownames(i))] = paste0(rownames(i)[duplicated(rownames(i))], "_2")
  i
}
createBarplot(normalise_rw(non_duplicated_rows(Lymph_CLL_obj$Y)), order_labels = names(sort(rowSums(non_duplicated_rows(Lymph_CLL_obj$Y)), decreasing = F)))


