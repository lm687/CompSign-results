## Comments on the cancer-specific results of PCAWG

##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../2_inference_TMB/helper_TMB.R")
source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
source("../3_analysis/recovery_COSMIC_signatures/recover_COSMIC_signatures.R")

library(TMB)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(dplyr)
library(jcolors)
library(viridis)
library(mutSigExtractor)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
enough_samples
nonexogenous = read.table("../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t",
                          comment.char = "#", fill = F)
nonexogenouswSBS1SBS5 = read.table("../../data/cosmic/exogenous_signatures_SBS_withSBS1SBS5.txt", sep = "\t",
                                   comment.char = "#", fill = F)
sigs_cosmic0 <- read.table(paste0( "../../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
                          stringsAsFactors = FALSE, sep = ',', header = TRUE)
rownames(sigs_cosmic0) <- paste0(substr(sigs_cosmic0$SubType, 1, 1),'[',
                                 sigs_cosmic0$Type, ']', substr(sigs_cosmic0$SubType, 3, 3))
sigs_cosmic0 <- sigs_cosmic0[-c(1,2)];
sigs_cosmic <- colnames(sigs_cosmic0)

nucleotide_colours_logR <- c('C$>$A/T$>$G'= '#3cb371', 'C$>$G/T$>$G'= '#90ee90', 'C$>$T/T$>$G'= '#66cdaa',
  'T$>$A/T$>$G'= '#cd5c5c', 'T$>$C/T$>$G'= '#f4a460')
nucleotide_colours <- c('C>A' = '#3cb371', 'C>G'= '#90ee90', 'C>T'= '#66cdaa',
                        'T>A'= '#cd5c5c', 'T>C'= '#f4a460', 'T>G'='red')
nucleotide_colours_dollar <- c('C$>$A' = '#3cb371', 'C$>$G'= '#90ee90', 'C$>$T'= '#66cdaa',
                        'T$>$A'= '#cd5c5c', 'T$>$C'= '#f4a460', 'T$>$G'='red')

pcawg_meta <- read.table("../../data/restricted/pcawg/repository_1567600367.tsv", sep = "\t", head=T)
metastatic_samples <- pcawg_meta[grepl('Metastatic', pcawg_meta$Specimen.Type),]
pcawg_meta2 <- read.table("../../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt", sep = "\t", head=T)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
read_info <- function(ct){
  .x <- list(fullRE_M_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct, "_signaturesPCAWG.RDS"))),
             fullRE_DMSL_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambda_", ct, "_signaturesPCAWG.RDS"))),
             fullRE_M_nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullRE_DMSL_nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             diagRE_DMDL_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS"))),
             diagRE_DMDL_nonexo_SP =  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             diagRE_DMDL_wSBS1SBS5nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMwSBS1SBS5nonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMnoscaling_SP_nonexo =  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMnoscaling_SP_nonexo_subsets_and_amalgamations <- try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMonefixedlambdanonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMonefixedlambda2nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambda2nonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMonefixedlambdanonexo_SPSaA = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWGSaA.RDS"))),
             fullREM_MSE = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct, "_signaturesMSE.RDS"))),
             fullREDM_MSE = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", ct, "_signaturesMSE.RDS"))),
             fullREDM_nucleotide1 = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", ct, "_nucleotidesubstitution1.RDS"))),
             diagREDM_MSE = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_", ct, "_signaturesMSE.RDS"))),
             dataset_all_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../data/", load_all_sigs = T, override_warning_X_Z = T),
             dataset_active_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../data/", override_warning_X_Z = T),
             dataset_nucleotidesubstitution1 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution1", path_to_data = "../../data/", override_warning_X_Z = T),
             dataset_nucleotidesubstitution3 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution3", path_to_data = "../../data/", override_warning_X_Z = T),
             dataset_nucleotidesubstitution3MSE = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution3MSE", path_to_data = "../../data/", override_warning_X_Z = T),
             dataset_active_sigs_MSE = load_PCAWG(ct = ct, typedata = "signaturesMSE", path_to_data = "../../data/", load_all_sigs = F, override_warning_X_Z = T),
             DMM = list(z_DMM=lapply(1:8, function(k) try(read.table(paste0("../../data/roo_for_DMM_SPpcawg/DMM_output/", ct, "_signaturesPCAWG_all", k, "_dmm.z"), sep = ',', skip = 1))),
                        fit_DMM = lapply(1:8, function(k) try(read.table(paste0("../../data/roo_for_DMM_SPpcawg/DMM_output/", ct, "_signaturesPCAWG_all", k, "_dmm.fit"), sep = ' '))))
               )
  .x$dataset_nonexo <- give_subset_sigs_TMBobj(.x$dataset_active_sigs, nonexogenous$V1)
  .x$dataset_nonexoSBS1SBS5 <- give_subset_sigs_TMBobj(.x$dataset_active_sigs, nonexogenouswSBS1SBS5$V1)
  .x$colnames_notsorted_SP <- try(colnames(.x$dataset_active_sigs$Y))
  .x$logR_notsorted_SP <- try(vector_cats_to_logR(.x$colnames_notsorted_SP))
  
  .x$colnames_nonexo_notsorted_SP <- try(colnames(.x$dataset_nonexo$Y))
  .x$logR_nonexo_notsorted_SP <- try(vector_cats_to_logR(.x$colnames_nonexo_notsorted_SP))
  
  ## nonexo, with SBS1 and SBS5
  .x$colnames_wSBS1SBS5nonexo_notsorted_SP <- try(colnames(.x$dataset_nonexoSBS1SBS5$Y))
  .x$logR_wSBS1SBS5nonexo_notsorted_SP <- try(vector_cats_to_logR(.x$colnames_wSBS1SBS5nonexo_notsorted_SP))
  return(.x)
}
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
read_info_list <- lapply(enough_samples, read_info)
names(read_info_list) <- enough_samples
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../3_analysis/helper/pcawg.colour.palette.R")
pcawg_palette <- pcawg.colour.palette(x = gsub("\\..*", "", names(read_info_list)),  scheme = "tumour.subtype")
names(pcawg_palette) <- names(read_info_list)
##-----------------------------------------------------------------------------------------------------##


plot_betas_scatter <- function(TMB_obj, names_cats=NULL, softmax=F){
  betas_mat <- t(matrix(python_like_select_name(TMB_obj$par.fixed, 'beta'), nrow=2))
  if(softmax){
    betas_mat <- t(softmax(t(rbind(betas_mat, 0))))
    names_cats <- c(sapply(names_cats, function(i) strsplit(i, '/')[[1]][1]),
                    strsplit(names_cats[1], '/')[[1]][2])
  }
  betas_df <- data.frame(betas_mat)
  if(!is.null(names_cats)){
    betas_df$logR = names_cats
  }else{
    betas_df$logR <- 1:nrow(betas_df)
  }
  a <- ggplot(betas_df, aes(x=X1, y=X2, label=logR))+labs(x='Beta intercept', y='Beta slope')
  if(!TMB_obj$pdHess){
    a <- a+geom_point(col='red')
  }else{
    a <- a+geom_point()
  }
  a <- a + geom_label_repel(size=3)+theme_bw()
  return(a)
}

give_mycolours <- function(n){
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  myColors <- col_vector[1:n]
}

plot_scatter <- function(m){
  m <- data.frame(m)
  names_sigs <- colnames(m)
  colnames(m) <- c('x1', 'x2')
  ggplot(m, aes(x=x1, y=x2))+geom_point()+theme_bw()+labs(x=names_sigs[1], y=names_sigs[2])
}

plot_covariance_mat <- function(TMB_obj, names_cats=NULL, title='Covariances', arg_cluster_rows=T,
                                arg_cluster_cols=T, lims=NULL){
  .x <- L_to_cov(python_like_select_name(TMB_obj$par.fixed, 'cov_par_RE'),
                 d=length(python_like_select_name(TMB_obj$par.fixed, 'beta'))/2)
  if(!is.null(names_cats))  colnames(.x) <- rownames(.x) <- names_cats
  if(!is.null(lims)){
    require(RColorBrewer)
    breaksList = seq(lims[1], lims[2], length.out = 10)
    hm <- pheatmap::pheatmap(.x, main=title, cluster_rows = arg_cluster_rows,
                             cluster_cols = arg_cluster_cols,
                             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                             breaks = breaksList)
  }else{
    hm <- pheatmap::pheatmap(.x, main=title, cluster_rows = arg_cluster_rows, cluster_cols = arg_cluster_cols, )
  }
  return(arrangeGrob(hm[[4]]))
}

give_plot_differential_precision <- function(TMB_obj){
  ovrdisp <- cbind.data.frame( plot_lambdas(TMB_obj, return_df=T, plot=F))
  if(TMB_obj$pdHess){
    col='red'
  }else{
    col='black'
  }
  
  ovrdisp[ovrdisp$name == 'Lambda 1','name'] = 'Early'
  ovrdisp[ovrdisp$name == 'Lambda 2','name'] = 'Late'
  
  ggplot(ovrdisp, aes(x=name,  y=`Estimate`))+
    geom_point(col=col)+
    geom_errorbar(aes(ymin=`Estimate`-`Std..Error`,
                      ymax=`Estimate`+`Std..Error`), width=.1)+
    theme_bw()+
    labs(x='', y='Log lambda')
}

give_DMM_plot <- function(TMB_obj){
  if(all(sapply(read_info_list[[ct]]$DMM$z_DMM, typeof) == 'character')){
    ggplot()+theme_bw()
  }else{
    ggplot(cbind(melt(normalise_rw(TMB_obj$dataset_active_sigs$Y)),
                 grouptiming=rep(c('Clonal','Subclonal'), each=dim(TMB_obj$dataset_active_sigs$Y)[1]/2),
                 groupDMM=paste0('DMM ', apply(TMB_obj$DMM$z_DMM[[which.min(unlist(sapply(TMB_obj$DMM$fit_DMM, `[`, 2)))]][,2:3], 1, which.max))),
           aes(x=Var1, y=value, fill=Var2))+geom_bar(stat = "identity")+facet_wrap(.~grouptiming+groupDMM, ncol=2)+
      scale_fill_manual(name = "grp",values = give_mycolours(ncol(read_info_list[[ct]]$dataset_active_sigs$Y)))+theme_bw()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  }
}

make_rownames_unique <- function(i){
  rownames(i)[duplicated(rownames(i))] <- paste0(rownames(i)[duplicated(rownames(i))], '_2')
  i
}

logR_nucleotidesubs1 <- vector_cats_to_logR(sort(unlist( sapply(c('C', 'T'),
                                                                function(base){ paste0(base, '>', c('A', 'C', 'G', 'T')[! (c('A', 'C', 'G', 'T') == base)]) }))))

betas_df_to_softmax <- function(df_betas){
  .x <- softmax(c(df_betas$Estimate, 0))
  names(.x) <- c(gsub("/.*", "", df_betas$LogR), gsub("*./", "", df_betas$LogR[1]))
  .x
}


give_boxplot_normalised_for_sig <- function(TMB_obj, sig){
  if(length(sig) == 1){
    ggplot(reshape2::melt(split(normalise_rw(TMB_obj$Y[,sig]), f  = TMB_obj$x[,2])),
           aes(x=L1, y=value))+geom_boxplot(outlier.shape = NA)+geom_jitter()+theme_bw()+ggtitle(sig)+
      labs(x='Group', y='Exposure')
  }
}

##-----------------------------------------------------------------------------------------------------##

# if(opt$model == "fullREM"){
# TMB::compile("../2_inference_TMB/mm_multinomial/fullRE_ME_multinomial.cpp",  "-std=gnu++17")
# dyn.load(dynlib("../2_inference_TMB/mm_multinomial/fullRE_ME_multinomial"))
# # }else if(opt$model == "diagREM"){
# TMB::compile("2_inference_TMB/mm_multinomial/diagRE_ME_multinomial.cpp",  "-std=gnu++17")
# dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_ME_multinomial"))
# mod_model_name = "diagRE_M"
# # }else if(opt$model == "fullREDM"){
TMB::compile("../2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("../2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial"))
# mod_model_name = "fullRE_DM"
# }else if(opt$model == "diagREDM"){
TMB::compile("../2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("../2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial"))
# mod_model_name = "diagRE_DM"
# # }else if(opt$model =="fullREDMsinglelambda"){
# TMB::compile("../2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
# dyn.load(dynlib("../2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda"))
# mod_model_name = "fullREDMsinglelambda"
# # use_nlminb=T
# # }else if(opt$model =="diagREDMsinglelambda"){
# TMB::compile("2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
# dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda"))
# mod_model_name = "diagREDMsinglelambda"
# # use_nlminb=T
# # }else if(opt$model =="FEDMsinglelambda"){
# TMB::compile("2_inference_TMB/mm_multinomial/FE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
# dyn.load(dynlib("2_inference_TMB/mm_multinomial/FE_dirichletmultinomial_single_lambda"))
# mod_model_name = "FEDMsinglelambda"
# # }else if(opt$model =="fullREDMnoscaling"){
# TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomialnoscaling.cpp",  "-std=gnu++17")
# dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomialnoscaling"))
# mod_model_name = "fullREDMnoscaling"
TMB::compile("../2_inference_TMB/mm_multinomial/singleRE_dirichlet_multinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("../2_inference_TMB/mm_multinomial/singleRE_dirichlet_multinomial"))


fullRE_DMDL_features <- list()
amalgamated_extra <- list()
fullRE_DMDL_extra <- list()
diagRE_DMDL_extra <- list()
mutsigexposures <- list()
fullDM_MSE <- list()
fullDM_MSE_nonexo <- list()
diagDM_MSE <- list()
diagDM_MSE_nonexo <- list()
new_sigs <- list()
fullDM_newsigs <- list()
fits_sigs <- list()

##-----------------------------------------------------------------------------------------------------##

## no metastatic files
for( ct in enough_samples){
  cat(ct, '\n')
  print(table(rownames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% metastatic_samples$Object.ID))
  print(table(rownames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% gsub("\\.consensus.*", "", metastatic_samples$File.Name)[1:3]))
}


##-----------------------------------------------------------------------------------------------------##
## reports for all ct
for(ct in enough_samples){
  .a <- give_barplot_from_obj(read_info_list[[ct]]$dataset_active_sigs, legend_on = T, legend_bottom = T, plot=F, arg_title='')
  pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report1.pdf"), onefile=FALSE, height = 12)
  print(cowplot::plot_grid(cowplot::plot_grid(.a),
                     cowplot::plot_grid(plot_betas(read_info_list[[ct]]$fullRE_DMSL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP, title = 'fullRE_DMSL_SP'),
                                        plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP, title = 'diagRE_DMDL_SP'),
                                        plot_betas(read_info_list[[ct]]$diagRE_DMDL_nonexo_SP, names_cats = read_info_list[[ct]]$logR_nonexo_notsorted_SP, title = 'diagRE_DMDL_nonexo_SP'),
                                        plot_betas(read_info_list[[ct]]$diagRE_DMDL_wSBS1SBS5nonexo_SP, names_cats = read_info_list[[ct]]$logR_wSBS1SBS5nonexo_notsorted_SP, title = 'diagRE_DMDL_wSBS1SBS5nonexo_SP')),
                     cowplot::plot_grid(plot_betas_scatter(read_info_list[[ct]]$diagRE_DMDL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP),
                                        plot_betas_scatter(read_info_list[[ct]]$diagRE_DMDL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP, softmax=T)),
                     cowplot::plot_grid(cowplot::plot_grid(createBarplot(normalise_rw(non_duplicated_rows(read_info_list[[ct]]$dataset_active_sigs$Y)),
                                                      order_labels = names(sort(rowSums(non_duplicated_rows(read_info_list[[ct]]$dataset_active_sigs$Y)),
                                                                                decreasing = F)), remove_labels=T)+guides(fill=F)+ggtitle('Sorted by mut. toll'),
                                                      give_plot_differential_precision(read_info_list[[ct]]$diagRE_DMDL_SP), ncol=1),
                                        plot_covariance_mat(read_info_list[[ct]]$fullRE_DMSL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP)),
                     nrow=4, rel_heights = c(2,2,1,2)))
  dev.off()
  Sys.sleep(1)
}

for(ct in enough_samples){
  pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report3MSE.pdf"), onefile=F, height = 6)
  cowplot::plot_grid(plot_grid(plot_betas(read_info_list[[ct]]$fullREDM_MSE,
                                names_cats = vector_cats_to_logR(colnames(read_info_list[[ct]]$dataset_active_sigs_MSE$Y))),
                               give_plot_differential_precision(read_info_list[[ct]]$fullREDM_MSE)+ggtitle('Lambda fullRE DMDL'), nrow=1, rel_widths = c(2,1)),
                     give_barplot_from_obj(read_info_list[[ct]]$dataset_active_sigs_MSE, legend_on = T, legend_bottom = T, plot=F, arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     nrow=2)
  dev.off()
}
##-----------------------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------------------##
ct=enough_samples[1]
names_sigs = read_info_list[[ct]]$colnames_notsorted_SP
abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))
rel_path = "../../"
return_dataframe = T
sigextraction = sig_extraction_method
exposures_input = NULL

give_all_bleeding_summary <- function(sig_extraction_method){
  
  all_bleeding <- lapply(enough_samples, function(ct){
    .x <- give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                             abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                             rel_path = "../../", return_dataframe = T, sigextraction = sig_extraction_method)
    .x$ct <-  ct
    .x
  })

  
  all_bleeding <- do.call('rbind', all_bleeding)
  
  all_bleeding$bool = (all_bleeding$exposuressigsQP > all_bleeding$exposures)
  all_bleeding$bool[all_bleeding$exposuressigsQP  == all_bleeding$exposures] <- NA ## neither increase nor decrease
  table(is.na(all_bleeding$bool))
  table(all_bleeding$bool)
  # all_bleeding_summary$bool[all_bleeding_summary$bool] = 'Overestimate'
  # all_bleeding_summary$bool[all_bleeding_summary$bool == 'FALSE'] = 'Underestimate'
  
  all_bleeding_summary <- all_bleeding %>% group_by(sig, ct) %>% dplyr::summarise(mean_change=mean(bool, na.rm = T))
  all_bleeding_summary$shape = all_bleeding_summary$mean_change > 0.65
  all_bleeding_summary
}

all_bleeding_summary <- give_all_bleeding_summary('QP')
all_bleeding_summary_MutationalPatterns <- give_all_bleeding_summary('MutationalPatterns')
all_bleeding_summary_sigminer <- give_all_bleeding_summary('sigminer')

tikzDevice::tikz(file ="../../results/results_TMB/pcawg/reports_per_cancer_type/bleeding_signatures.tex",
                 width=7, height = 8)
ggplot((all_bleeding_summary), aes(x=ct, y=factor(sig, levels=rev(gtools::mixedsort(unique(all_bleeding_summary$sig)))),
                                   fill=mean_change))+geom_tile()+
  # scale_fill_viridis(begin = 0, end = 1)+
  scale_fill_gradientn(colours = c('red', 'yellow', 'blue'), values = c(0,0.5, 1))+
  # jcolors::scale_fill_jcolors_contin("pal3", reverse = TRUE, bias = 0.2)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(data = all_bleeding_summary[all_bleeding_summary$shape == T,])+
  labs(x='Cancer types', y='Signatures', fill='Fraction of overestimates over changed exposures')+
  theme(legend.position = "bottom")+guides(fill=FALSE)
dev.off()

tikzDevice::tikz(file ="../../results/results_TMB/pcawg/reports_per_cancer_type/bleeding_signatures_MutationalPatterns.tex",
                 width=7, height = 8)
ggplot((all_bleeding_summary_MutationalPatterns), aes(x=ct, y=factor(sig, levels=rev(gtools::mixedsort(unique(all_bleeding_summary_MutationalPatterns$sig)))),
                                   fill=mean_change))+geom_tile()+
  # scale_fill_viridis(begin = 0, end = 1)+
  scale_fill_gradientn(colours = c('red', 'yellow', 'blue'), values = c(0,0.5, 1))+
  # jcolors::scale_fill_jcolors_contin("pal3", reverse = TRUE, bias = 0.2)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(data = all_bleeding_summary_MutationalPatterns[all_bleeding_summary_MutationalPatterns$shape == T,])+
  labs(x='Cancer types', y='Signatures', fill='Fraction of overestimates over changed exposures')+
  theme(legend.position = "bottom")+guides(fill=FALSE)
dev.off()


tikzDevice::tikz(file ="../../results/results_TMB/pcawg/reports_per_cancer_type/bleeding_signatures_sigminer.tex",
                 width=7, height = 8)
ggplot((all_bleeding_summary_sigminer), aes(x=ct, y=factor(sig, levels=rev(gtools::mixedsort(unique(all_bleeding_summary_sigminer$sig)))),
                                            fill=mean_change))+geom_tile()+
  # scale_fill_viridis(begin = 0, end = 1)+
  scale_fill_gradientn(colours = c('red', 'yellow', 'blue'), values = c(0,0.5, 1))+
  # jcolors::scale_fill_jcolors_contin("pal3", reverse = TRUE, bias = 0.2)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(data = all_bleeding_summary_sigminer[all_bleeding_summary_sigminer$shape == T,])+
  labs(x='Cancer types', y='Signatures', fill='Fraction of overestimates over changed exposures')+
  theme(legend.position = "bottom")+guides(fill=FALSE)
dev.off()


all_bleeding_summary_df=cbind.data.frame(QP=all_bleeding_summary,
                 MutationalPatterns=all_bleeding_summary_MutationalPatterns[,3:4],
                 sigminer=all_bleeding_summary_sigminer[,3:4])
all_bleeding_summary_df_melt <- (melt(all_bleeding_summary_df[,c("QP.sig", "QP.ct", "QP.mean_change", "MutationalPatterns.mean_change", "sigminer.mean_change" )]))

ggplot(all_bleeding_summary_df, aes(x=QP.mean_change, y=MutationalPatterns.mean_change,
                                    col=QP.sig))+
  geom_point()+theme_bw()

get_ordered_sig <- function(i){
  i[order(i[,2]),1]
}
ggplot(all_bleeding_summary_df_melt, aes(x=factor(QP.sig, levels=all_bleeding_summary_df_melt %>% group_by(QP.sig) %>% summarise(median(value)) %>% get_ordered_sig %>% unlist
), y=value, group=interaction(QP.sig,variable), col=gsub(".mean_change", "", variable)))+
  geom_boxplot()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x='Signature', y='Fraction of samples with overestimated exposure', col='Fitting method')+
  theme(legend.position = "bottom")
ggsave("../../results/results_TMB/pcawg/reports_per_cancer_type/bleeding_comparison_methods_summary.pdf",
       width = 10, height = 5)
ggplot(all_bleeding_summary_df_melt, aes(x=factor(QP.ct, levels=all_bleeding_summary_df_melt %>% group_by(QP.ct) %>% summarise(median(value)) %>% get_ordered_sig %>% unlist
), y=value, group=interaction(QP.ct,variable)))+
  geom_boxplot()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x='Cancer type', y='Fraction of samples with overestimated exposure')
ggsave("../../results/results_TMB/pcawg/reports_per_cancer_type/bleeding_comparison_methods_summary_ct.pdf",
       width = 10, height = 4)

##-----------------------------------------------------------------------------------------------------##

## Bone-Osteosarc
ct <- "Bone-Osteosarc"
## in https://www.nature.com/articles/s41588-019-0390-2 they show two pops, one with and one without SBS3
createBarplot(normalise_rw(non_duplicated_rows(read_info_list[[ct]]$dataset_active_sigs$Y)),
              order_labels = names(sort(non_duplicated_rows(read_info_list[[ct]]$dataset_active_sigs$Y)[,'SBS3'],
                                        decreasing = F)), remove_labels=T)+ggtitle('Sorted by SBS3')
## run without the minor signatures


View(read_info_list[[ct]]$dataset_active_sigs$Y)
fullDM_no_small_sigs <- wrapper_run_TMB(obje<- = give_subset_sigs_TMBobj(read_info_list[[ct]]$dataset_active_sigs,
                           sigs_to_remove = c('SBS13', 'SBS17a', 'SBS17b', 'SBS30')),
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullDM_no_small_sigs
wald_TMB_wrapper(fullDM_no_small_sigs, fail_non_converged = F)
fullDM_no_small_sigs_nonexo <- wrapper_run_TMB(object = give_subset_sigs_TMBobj(read_info_list[[ct]]$dataset_active_sigs,
                                                                         sigs_to_remove = unique(c(c('SBS13', 'SBS17a', 'SBS17b', 'SBS30'),
                                                                                            nonexogenous$V1))),
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullDM_no_small_sigs_nonexo
wald_TMB_wrapper(fullDM_no_small_sigs_nonexo, fail_non_converged = F)

plot_scatter(read_info_list[[ct]]$dataset_active_sigs$Y[,c('SBS17a', 'SBS17b')])
plot_scatter(read_info_list[[ct]]$dataset_active_sigs$Y[,c('SBS2', 'SBS3')])


fullRE_DMDL_features =  wrapper_run_TMB(object = read_info_list[[ct]]$dataset_nucleotidesubstitution1,
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(cowplot::plot_grid(plot_betas(fullRE_DMDL_features, names_cats = c(vector_cats_to_logR(sort(unlist( sapply(c('C', 'T'),
                          function(base){ paste0(base, '>', c('A', 'C', 'G', 'T')[! (c('A', 'C', 'G', 'T') == base)]) })))))),
                          plot_betas(fullDM_no_small_sigs, names_cats = vector_cats_to_logR(colnames(give_subset_sigs_TMBobj(read_info_list[[ct]]$dataset_active_sigs,
                                                                                                                             sigs_to_remove = c('SBS13', 'SBS17a', 'SBS17b', 'SBS30'))$Y)), title = "Removing various sigs")),
                         give_DMM_plot(read_info_list[[ct]]),
                         give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                                            rel_path = "../../"), nrow=3, rel_heights = c(1,2,2)))
dev.off()


mutsigexposures[[ct]] = load_PCAWG(ct = ct, typedata = "signaturesMSE",
                                   path_to_data = "../../data/", load_all_sigs = F)

dim(mutsigexposures[[ct]]$Y)
dim(read_info_list[[ct]]$dataset_active_sigs$Y)
createBarplot(normalise_rw(non_duplicated_rows(mutsigexposures[[ct]]$Y)),
              order_labels = names(sort(non_duplicated_rows(mutsigexposures[[ct]]$Y)[,'SBS3'],
                                        decreasing = F)), remove_labels=T)+ggtitle('Sorted by SBS3')

read_info_list[[ct]]$fullREDM_MSE
read_info_list[[ct]]$diagREDM_MSE

# fullDM_MSE[[ct]] <- wrapper_run_TMB(object = mutsigexposures[[ct]],
#                                         model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
# saveRDS(object = fullDM_MSE[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_MSE', '.RDS'))
# 
# saveRDS(object = fullDM_MSE[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_MSEnonexo', '.RDS'))
# 
# .a <- give_barplot_from_obj(mutsigexposures[[ct]], legend_on = T, legend_bottom = T, plot=F, arg_title='')
# 
# pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report3MSE.pdf"), onefile=FALSE, height = 12)
# print(cowplot::plot_grid(cowplot::plot_grid(.a),
#                          plot_betas(fullDM_MSE[[ct]], names_cats = c(vector_cats_to_logR(colnames(mutsigexposures[[ct]]$Y)))),
#                          nrow = 2))
# dev.off()
# 
# wald_TMB_wrapper(fullDM_MSE[[ct]])
read_info_list[[ct]]$dataset_active_sigs_MSE
new_sigs[[ct]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                     subset_signatures = c('SBS1', 'SBS13', 'SBS2', 'SBS3', 'SBS5', 'SBS8', 'SBS40'))

new_sigs[[paste0(ct, '_1')]] <- read_info_list[[ct]]$dataset_active_sigs
new_sigs[[paste0(ct, '_2')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c('SBS1','SBS2', 'SBS5', 'SBS8', 'SBS13',  'SBS40'))
new_sigs[[paste0(ct, '_3')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c('SBS1','SBS2', 'SBS3', 'SBS5', 'SBS8', 'SBS13'))
new_sigs[[paste0(ct, '_4')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c('SBS1', 'SBS2', 'SBS5', 'SBS6', 'SBS13', 'SBS30'))
new_sigs[[paste0(ct, '_5')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c('SBS1', 'SBS2', 'SBS5', 'SBS13', 'SBS30'))
new_sigs[[paste0(ct, '_6')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c('SBS1', 'SBS2', 'SBS5', 'SBS13'))
new_sigs[[paste0(ct, '_allsigs')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                       subset_signatures = colnames(SBS_SIGNATURE_PROFILES_V3))
pheatmap::pheatmap(new_sigs[[paste0(ct, '_allsigs')]]$Y)

fits_sigs[[ct]] <- lapply(1:6, function(it) compare_signaturefit_to_data(new_sigs[[paste0(ct, paste0('_', it))]],
                                                read_info_list[[ct]]$dataset_nucleotidesubstitution3, sigs_cosmic0))

plot(unlist(sapply(fits_sigs[[ct]], `[`, 'cossim')), type='h', ylim = c(0,1))
plot(unlist(sapply(fits_sigs[[ct]], `[`, 'rss')), type='l')

give_plot_fits_wrapper(lapply(1:6, function(it) new_sigs[[paste0(ct, paste0('_', it))]]),
                       read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                       sigs_cosmic0)
ggsave(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, '_fits_reextractionsubsetsigs', '.pdf'), height = 3, width = 7)

plot_dataset_obj_mean_variance <- function(dataset_obj){
  ggplot(data.frame(sd=apply(dataset_obj$Y, 2, sd), mean= colSums(dataset_obj$Y),
             sig=colnames(dataset_obj$Y)), aes(x=mean, y=sd, label=sig))+
    geom_point()+geom_label_repel(max.overlaps = Inf)+scale_x_continuous(trans = "log2")+theme_bw()
}

plot_dataset_obj_mean_variance(new_sigs[[ct]])
plot_dataset_obj_mean_variance(new_sigs[[paste0(ct, '_allsigs')]])


remove_na <- function(i) i[!is.na(i)]

give_ggplot_sig_cor_TMB(new_sigs[[ct]], read_info_list[[ct]]$dataset_active_sigs_MSE)
colnames(new_sigs[[ct]]$Y)
colnames(read_info_list[[ct]]$dataset_active_sigs_MSE$Y)

diagDM_newsigs[[ct]] <- wrapper_run_TMB(object = new_sigs[[ct]],
                                        model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_newsigs[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_2')]],
                                        model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_newsigs[[paste0(ct, '_3')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_3')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_newsigs[[paste0(ct, '_4')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_4')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_newsigs[[paste0(ct, '_5')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_5')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_newsigs[[paste0(ct, '_6')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_6')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report4sigreextract.pdf"),
    onefile=FALSE, height = 12)
cowplot::plot_grid(cowplot::plot_grid(
  give_barplot_from_obj(new_sigs[[ct]], legend_on = T, legend_bottom = T, plot=F,
                                         arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
  plot_betas(diagDM_newsigs[[ct]], names_cats = vector_cats_to_logR(colnames(new_sigs[[ct]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_2')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_2')]],
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_2')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_3')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_3')]],
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_3')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_4')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_4')]],
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_4')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_5')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_5')]],
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_6')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_6')]],
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_6')]]$Y)))),
  nrow=6)
dev.off()

## Only 1,2,3,5,8
new_sigs[[paste0(ct, '_select1')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                           subset_signatures = c('SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS8', 'SBS13', 'SBS30', 'SBS40'))
new_sigs[[paste0(ct, '_select1')]] <- give_amalgamated_exposures_TMBobj((new_sigs[[paste0(ct, '_select1')]]),
                           list_groupings = c(list(c('SBS2', 'SBS13')),
                                              as.list(colnames(new_sigs[[paste0(ct, '_select1')]]$Y)[!(colnames(new_sigs[[paste0(ct, '_select1')]]$Y)
                                               %in% c('SBS2', 'SBS13'))])))
diagDM_newsigs[[paste0(ct, '_select1')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_select1')]],
                                                            model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)

new_sigs[[paste0(ct, '_select1')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                           subset_signatures = c('SBS1', 'SBS2', 'SBS13', 'SBS30', 'SBS5'))
fullDM_newsigs[[paste0(ct, '_select1')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_select1')]],
                                                            model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_newsigs[[paste0(ct, '_select1')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_select1')]],
                                                            model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(fullDM_newsigs[[paste0(ct, '_select1')]], names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_select1')]]$Y)))

plot_betas(diagDM_newsigs[[paste0(ct, '_select1')]], names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_select1')]]$Y)))

##-----------------------------------------------------------------------------------------------------##
ct <- "Breast-AdenoCA"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 15000)))


fullRE_DMDL_features[[ct]] =  wrapper_run_TMB(object = read_info_list[[ct]]$dataset_nucleotidesubstitution1,
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

fullRE_DMDL_nonhypermutated =  wrapper_run_TMB(object =give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = fullRE_DMDL_nonhypermutated, file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))
fullRE_DMDL_nonhypermutated <- readRDS( file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))

saveRDS(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut), file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'ROO_nonhypermutated', '.RDS'))

reorder_cols_first_to_last_TMBobj <- function(TMB_obj){
  .x <- colnames(TMB_obj$Y)
  TMB_obj$Y <- cbind(TMB_obj$Y[,-1], TMB_obj$Y[,1])
  colnames(TMB_obj$Y) <- c(.x[-1], .x[1])
  TMB_obj
}

amalgamated_extra[[ct]] <- reorder_cols_first_to_last_TMBobj(give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                  list_groupings = c(list(c('SBS9', 'SBS17a', 'SBS40')),
                                                          as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS9', 'SBS17a', 'SBS40'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
diagRE_DMDL_extra[[ct]]

plot(density(read_info_list[[ct]]$dataset_active_sigs$Y[,'SBS37']))

amalgamated_extra[[paste0(ct, '_2')]] <- reorder_cols_first_to_last_TMBobj(give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                                                               list_groupings = c(list(c('SBS9', 'SBS17a', 'SBS40', 'SBS1')),
                                                                                                                  list(c('SBS5', 'SBS17b', 'SBS37', 'SBS41', 'SBS3')),
                                                                                                                  as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS1', 'SBS9', 'SBS17a', 'SBS40', 'SBS5', 'SBS17b', 'SBS37', 'SBS41', 'SBS3'))]))))
diagRE_DMDL_extra[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_2')]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)

plot_betas(diagRE_DMDL_extra[[paste0(ct, '_2')]], names_cats = vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_2')]]$Y)))

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(cowplot::plot_grid(plot_betas(fullRE_DMDL_features[[ct]],
                                                       names_cats = logR_nucleotidesubs1),
                                            plot_betas(fullRE_DMDL_nonhypermutated, names_cats = read_info_list[[ct]]$logR_notsorted_SP, title = "Removing hypermut samples"), rel_widths = c(2,3)),
                         cowplot::plot_grid(createBarplot(normalise_rw(non_duplicated_rows(read_info_list[[ct]]$dataset_active_sigs$Y)),
                                                          order_labels = names(sort(non_duplicated_rows(read_info_list[[ct]]$dataset_active_sigs$Y)[,'SBS37'])),
                                                          remove_labels=T)+guides(fill=F)+ggtitle('Sorted by SBS37'),
                                            plot_betas(diagRE_DMDL_extra[[ct]],
                                    names_cats = vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='diagRE_DMDL amalgamation')),
                         give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                                            rel_path = "../../"),
                         plot_grid(give_plot_bleeding(names_sigs =  c('SBS9', 'SBS40'),
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS9', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                            rel_path = "../../", scales_x_free = T),
                         give_plot_bleeding(names_sigs =  c('SBS9', 'SBS17a', 'SBS40'),
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS9','SBS17a', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                            rel_path = "../../", scales_x_free = T), ncol=2),
                         nrow=4, rel_heights = c(1,2,2, 1)))
dev.off()

plot(colSums(read_info_list[[ct]]$dataset_active_sigs$Y),
colSums(read_info_list[[ct]]$dataset_active_sigs_MSE$Y))
plot(as.vector(read_info_list[[ct]]$dataset_active_sigs$Y),
     as.vector(read_info_list[[ct]]$dataset_active_sigs_MSE$Y))

new_sigs[[paste0(ct, '_1')]] <- read_info_list[[ct]]$dataset_active_sigs
new_sigs[[paste0(ct, '_2')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS9", "SBS13",
                                                                           "SBS17a", "SBS17b", "SBS18",  "SBS37",  "SBS41" ))
new_sigs[[paste0(ct, '_2MSE')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3MSE,
                                                     subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS9", "SBS13",
                                                                           "SBS17a", "SBS17b", "SBS18",  "SBS37",  "SBS41" ))
new_sigs[[paste0(ct, '_3')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS9", "SBS13",
                                                                           "SBS18",  "SBS37",  "SBS41" ))
new_sigs[[paste0(ct, '_4')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS9", "SBS13",
                                                                           "SBS18",  "SBS41" ))
new_sigs[[paste0(ct, '_5')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS13",
                                                                           "SBS18",  "SBS41" ))
new_sigs[[paste0(ct, '_6')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS9", "SBS13",
                                                                           "SBS18",  "SBS41" ))
new_sigs[[paste0(ct, '_7')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS9",
                                                                           "SBS18",  "SBS41" ))
new_sigs[[paste0(ct, '_8')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1", "SBS13", "SBS3", "SBS5", "SBS8", "SBS9",
                                                                           "SBS18",  "SBS41" ))
new_sigs[[paste0(ct, '_9')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1", "SBS2", "SBS13", "SBS5", "SBS8", "SBS9",
                                                                           "SBS18",  "SBS41" ))
# new_sigs[[paste0(ct, '_10')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
#                                                      subset_signatures = c("SBS2", "SBS3", "SBS5", "SBS8", "SBS9",
#                                                                            "SBS13", "SBS18",  "SBS41" ))
new_sigs[[paste0(ct, '_10')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13",
                                                                           "SBS18",  "SBS41" ))
new_sigs[[paste0(ct, '_11')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                      subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13",
                                                                            "SBS18" ))
new_sigs[[paste0(ct, '_12')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                      subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13",
                                                                            "SBS41" ))
new_sigs[[paste0(ct, '_13')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                      subset_signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13" ))
plot(unlist(sapply(fits_sigs[[ct]], `[`, 'cossim')), type='h', ylim = c(0,1))
plot(unlist(sapply(fits_sigs[[ct]], `[`, 'rss')), type='l')


give_plot_fits_wrapper(lapply(1:10, function(it) new_sigs[[paste0(ct, paste0('_', it))]]),
                       read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                       sigs_cosmic0)+   guides(col=guide_legend(ncol=2,byrow=TRUE))
give_plot_fits_wrapper(lapply(c(1:6, 9:13), function(it) new_sigs[[paste0(ct, paste0('_', it))]]),
                       read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                       sigs_cosmic0)+   guides(col=guide_legend(ncol=2,byrow=TRUE))
ggsave(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, '_fits_reextractionsubsetsigs', '.pdf'), height = 3, width = 7)

diagDM_newsigs[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_2')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(diagDM_newsigs[[paste0(ct, '_2')]],  file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_reextractionsubsetsigs2', '.RDS'))
diagDM_newsigs[[paste0(ct, '_2MSE')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_2MSE')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)

diagDM_newsigs[[paste0(ct, '_3')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_3')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(diagDM_newsigs[[paste0(ct, '_3')]],  file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_reextractionsubsetsigs3', '.RDS'))

diagDM_newsigs[[paste0(ct, '_4')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_4')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(diagDM_newsigs[[paste0(ct, '_4')]],  file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_reextractionsubsetsigs4', '.RDS'))

diagDM_newsigs[[paste0(ct, '_10')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_10')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(diagDM_newsigs[[paste0(ct, '_10')]],  file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_reextractionsubsetsigs4', '.RDS'))

diagDM_newsigs[[paste0(ct, '_11')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_11')]],
                                                       model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(diagDM_newsigs[[paste0(ct, '_11')]],  file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_reextractionsubsetsigs4', '.RDS'))


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report4sigreextract.pdf"),
    onefile=FALSE, height = 13, width = 7.5)
cowplot::plot_grid(
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_1')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, title = 'PCAWGS sigs',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_1')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_2')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_2')]], title = 'Removing SBS40',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_2')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_2MSE')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_2MSE')]], title = 'Removing SBS40 (MSE extraction)',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_2MSE')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_3')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_3')]], title = 'Removing SBS40 and SBS17a/b',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_3')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_4')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_4')]], title = 'Removing SBS40, SBS17a/b, and SBS37',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_4')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_10')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_10')]], title = 'Removing SBS40, SBS17a/b, SBS37, SBS8, SBS9',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_10')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_11')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_11')]], title = 'Removing SBS40, SBS17a/b,  SBS37, SBS8, SBS9, SBS41',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_11')]]$Y)))),
  nrow=7)
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "CNS-GBM"

hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 15000)))
fullRE_DMDL_nonhypermutated =  wrapper_run_TMB(object =give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                               model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_nonhypermutated =  wrapper_run_TMB(object =give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = fullRE_DMDL_nonhypermutated, file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))
saveRDS(object = diagRE_DMDL_nonhypermutated, file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_nonhypermutated', '.RDS'))

plot_betas(diagRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP)

amalgamated_extra[[ct]] <- give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                                                               list_groupings = c(list(c('SBS11', 'SBS30')),
                                                                                                                  as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS11', 'SBS30'))])))
fullRE_DMDL_nonhypermutated_amalgamated =  wrapper_run_TMB(object =give_subset_samples_TMBobj( amalgamated_extra[[ct]], samples_to_remove = hypermut),
                                               model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(plot_betas(fullRE_DMDL_nonhypermutated_amalgamated, vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y))),
                         give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                                            rel_path = "../../"),
                         cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS11'),
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS11'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                            rel_path = "../../"),
                                            give_plot_bleeding(names_sigs = c('SBS5', 'SBS30'),
                                                               abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS30'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                               rel_path = "../../"),
                                            give_plot_bleeding(names_sigs = c('SBS5', 'SBS37'),
                                                               abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS37'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                               rel_path = "../../"), nrow=1
                                            ),
                         nrow=3, rel_heights = c(1,2, 1)))
dev.off()

colnames(read_info_list[[ct]]$dataset_active_sigs$Y)
new_sigs[[paste0(ct, '_1')]] <- read_info_list[[ct]]$dataset_active_sigs
new_sigs[[paste0(ct, '_2')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                      subset_signatures = c("SBS1",  "SBS5",  "SBS11", "SBS30", "SBS37"))
new_sigs[[paste0(ct, '_5')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS11", "SBS30", "SBS37"))
new_sigs[[paste0(ct, '_3')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS5",  "SBS11", "SBS30"))
new_sigs[[paste0(ct, '_4')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS5",  "SBS11"))
give_plot_fits_wrapper(lapply(1:5, function(it) new_sigs[[paste0(ct, paste0('_', it))]]),
                       read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                       sigs_cosmic0)+   guides(col=guide_legend(ncol=2,byrow=TRUE))
ggsave(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, '_fits_reextractionsubsetsigs', '.pdf'), height = 3, width = 7)

diagDM_newsigs[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_2')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_newsigs[[paste0(ct, '_3')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_3')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_newsigs[[paste0(ct, '_4')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_4')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_newsigs[[paste0(ct, '_5')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_5')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report4sigreextract.pdf"),
    onefile=FALSE, height = 11, width = 7.5)
cowplot::plot_grid(
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_1')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, title = 'PCAWG sigs',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_1')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_2')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_2')]], title = 'Removing SBS40',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_2')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_3')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_3')]], title = 'Removing SBS40 or 37',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_3')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_4')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_4')]], title = 'Removing SBS40 or 37',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_4')]]$Y)))),
  cowplot::plot_grid(give_barplot_from_obj(new_sigs[[paste0(ct, '_5')]], legend_on = T, legend_bottom = T, plot=F,
                                           arg_title='', only_normalised = T, nrow_plot = 1, levels_signatures=sigs_cosmic),
                     plot_betas(diagDM_newsigs[[paste0(ct, '_5')]], title = 'Removing SBS40 or 37 or 5',
                                names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)))),
  nrow=5)
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "CNS-Medullo"

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
                         give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                                            rel_path = "../../"),
                         cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS18'),
                                                               abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS11'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                               rel_path = "../../"),
                                            give_plot_bleeding(names_sigs = c('SBS8', 'SBS30'),
                                                               abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS30'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                               rel_path = "../../"), nrow=1
                         ),
                         nrow=2, rel_heights = c(2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "CNS-PiloAstro"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 600)))
fullRE_DMDL_nonhypermutated =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                               model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                                                               list_groupings = c(list(c('SBS19', 'SBS23')),
                                                                                                                  as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS19', 'SBS23'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
wald_TMB_wrapper(fullRE_DMDL_extra[[ct]])

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y))),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS18'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS11'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"),
  #                    give_plot_bleeding(names_sigs = c('SBS8', 'SBS30'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS30'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"), nrow=1
  # ),
  nrow=2, rel_heights = c(1, 2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "ColoRect-AdenoCA"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 1e5)))
fullRE_DMDL_nonhypermutated =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                               model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_nonhypermutated

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS28', 'SBS40', 'SBS45')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS28', 'SBS40', 'SBS45'))]))))
colSums(amalgamated_extra[[ct]]$Y)
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  plot_betas(fullRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                   abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                   rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS45', 'SBS40'),
                                      abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS45', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                                            rel_path = "../../"),
                   give_plot_bleeding(names_sigs = c('SBS28', 'SBS40'),
                                      abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS28', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                                            rel_path = "../../")),
  nrow=3, rel_heights = c(1, 2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Eso-AdenoCA"

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  # plot_betas(fullRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS3', 'SBS40'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
  give_plot_bleeding(names_sigs = c('SBS28', 'SBS17b'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS28', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../")),
  nrow=3, rel_heights = c(2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Head-SCC"

strange_sample <- names(sort((read_info_list[[ct]]$dataset_active_sigs$Y[,'SBS33']), dec=T)[1])
diagRE_DMDL_nostrange_sample =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = strange_sample),
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_nostrange_sample

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  plot_betas(diagRE_DMDL_nostrange_sample, read_info_list[[ct]]$logR_notsorted_SP, title = "without c0812962-a345-48b7-aec0-01336c2d1eed"),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  
  plot_scatter(matrix(read_info_list[[ct]]$dataset_active_sigs$Y[,c('SBS33')]/rowSums(read_info_list[[ct]]$dataset_active_sigs$Y), ncol=2)[!grepl(strange_sample, unique(rownames(read_info_list[[ct]]$dataset_active_sigs$Y))),])+labs(x='Clonal SBS33', y='Subclonal SBS33'),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS33'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS33'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS5', 'SBS16'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS16'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS3', 'SBS4'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS4'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS17a', 'SBS40'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS17a', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../")),
  nrow=4, rel_heights = c(1, 2, 1, 1)))
dev.off()

plot_scatter(read_info_list[[ct]]$dataset_active_sigs$Y[,c('SBS1', 'SBS33')])+scale_y_continuous(trans = "log2")+
  scale_x_continuous(trans = "log2")



##-----------------------------------------------------------------------------------------------------##

ct <- "Kidney-ChRCC"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 4000)))
diagRE_DMDL_nonhypermutated =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_nonhypermutated

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS17a', 'SBS17b', 'SBS5', 'SBS29')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS17a', 'SBS17b', 'SBS5', 'SBS29'))]))))
diagRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(diagRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP, title='Without hypermutated sample'),
                     plot_betas(diagRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Amalgamation without hypermutated sample'), nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # plot_scatter(matrix(read_info_list[[ct]]$dataset_active_sigs$Y[,c('SBS33')]/rowSums(read_info_list[[ct]]$dataset_active_sigs$Y), ncol=2)[!grepl(strange_sample, unique(rownames(read_info_list[[ct]]$dataset_active_sigs$Y))),])+labs(x='Clonal SBS33', y='Subclonal SBS33'),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS17a'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS17a'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS5', 'SBS17b'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS1', 'SBS17a'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS17a'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS1', 'SBS17b'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../")),
  nrow=3, rel_heights = c(1, 2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Kidney-RCC.clearcell"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 10000)))
diagRE_DMDL_nonhypermutated =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_nonhypermutated

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(diagRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP, title='Without hypermutated sample'),
                     # plot_betas(diagRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Amalgamation without hypermutated sample'),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # plot_scatter(matrix(read_info_list[[ct]]$dataset_active_sigs$Y[,c('SBS33')]/rowSums(read_info_list[[ct]]$dataset_active_sigs$Y), ncol=2)[!grepl(strange_sample, unique(rownames(read_info_list[[ct]]$dataset_active_sigs$Y))),])+labs(x='Clonal SBS33', y='Subclonal SBS33'),
  # cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS17a'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS17a'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"),
  #                    give_plot_bleeding(names_sigs = c('SBS5', 'SBS17b'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"),
  #                    give_plot_bleeding(names_sigs = c('SBS1', 'SBS17a'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS17a'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"),
  #                    give_plot_bleeding(names_sigs = c('SBS1', 'SBS17b'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../")),
  nrow=3, rel_heights = c(1, 2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Kidney-RCC.papillary"

# kidney papillary needs SBS13 removed, because it's in extremely low quantity and probably many zeros

give_plot_bleeding_wrapper <- function(sigs_bleed, ...){
  give_plot_bleeding(names_sigs = sigs_bleed,
                   abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( sigs_bleed, read_info_list[[ct]]$colnames_notsorted_SP)],
                   rel_path = "../../", resort_sig_labels = F, ...)
}

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj( read_info_list[[ct]]$dataset_active_sigs,
                                                              list_groupings = c(list(c('SBS13', 'SBS2'),
                                                                                      c('SBS5', 'SBS29'),
                                                                                       c( 'SBS40', 'SBS22')),
    as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in%
                                                                  c('SBS13', 'SBS2', 'SBS22', 'SBS29', 'SBS5', 'SBS40'))]))))
colSums(amalgamated_extra[[ct]]$Y)
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
saveRDS(object = fullRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_subset', '.RDS'))
wald_TMB_wrapper(fullRE_DMDL_extra[[ct]])

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS13', 'SBS5')),
                     give_plot_bleeding_wrapper( c('SBS13', 'SBS40')), ncol=2),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS29', 'SBS5')),
                     give_plot_bleeding_wrapper( c('SBS29', 'SBS40')), ncol=2),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
                                title='Amalgamation'),
                     nrow=1),
  nrow=4, rel_heights = c(2, 1, 1, 2)))
dev.off()


##-----------------------------------------------------------------------------------------------------##

ct <- "Liver-HCC"


colnames(dataset_subset$Y)
ncol(dataset_subset$Y)
colnames(read_info_list[[ct]]$dataset_active_sigs$Y)
ncol(read_info_list[[ct]]$dataset_active_sigs$Y)

typedata='signatures'
dataset_subset <- give_subset_sigs_TMBobj( read_info_list[[ct]]$dataset_active_sigs,
                  sigs_to_remove = c('SBS14', 'SBS17a', 'SBS17b', 'SBS18', 'SBS24', 'SBS28', 'SBS31', 'SBS35', 'SBS53', 'SBS54', 'SBS9'))
diagRE_DMDL_subset =  wrapper_run_TMB(object = dataset_subset,
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_subset
saveRDS(object = diagRE_DMDL_subset, file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_subset', '.RDS'))

## there are many active signatures
load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../data/")
length(python_like_select_name(read_info_list[[ct]]$fullRE_M_SP$par.fixed, 'beta'))/2
length(python_like_select_name(diagRE_DMDL_subset$par.fixed, 'beta'))/2


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(plot_betas(diagRE_DMDL_subset, vector_cats_to_logR(colnames(dataset_subset$Y)), title='Subset of sigs'),
                     nrow=1),
  nrow=2, rel_heights = c(2, 1)))
dev.off()

new_sigs[[ct]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                       subset_signatures = c('SBS1', 'SBS4', 'SBS5', 'SBS6', 'SBS12', 'SBS16', 'SBS19', 'SBS22', 'SBS29', 'SBS30', 'SBS40'))

diagDM_newsigs[[ct]] <- wrapper_run_TMB(object = new_sigs[[ct]],
                                        model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[ct]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[ct]]$Y)))

## comparison of betas between two tmb runs


compare_betas_tmb(tmb_obj_1 = read_info_list[[ct]]$diagRE_DMDL_SP,
                  names_cats1 = (colnames(read_info_list[[ct]]$dataset_active_sigs$Y)),
                  tmb_obj_2 = diagDM_newsigs[[ct]] ,
                  names_cats2 = (colnames(new_sigs[[ct]]$Y)))

## removing SBS5, which is found in low abundance
new_sigs[[paste0(ct, '_2')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                       subset_signatures = c('SBS1', 'SBS4', 'SBS6', 'SBS12', 'SBS16', 'SBS19', 'SBS22', 'SBS29', 'SBS30', 'SBS40'))

diagDM_newsigs[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_2')]],
                                        model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_2')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_2')]]$Y)))

## removing SBS40
new_sigs[[paste0(ct, '_4')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c('SBS1', 'SBS4', 'SBS5', 'SBS6', 'SBS12', 'SBS16', 'SBS19', 'SBS22', 'SBS29', 'SBS30'))

diagDM_newsigs[[paste0(ct, '_4')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_4')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_4')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_4')]]$Y)))

## further to SBS5, removing SBS40
new_sigs[[paste0(ct, '_3')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c('SBS1', 'SBS4', 'SBS6', 'SBS12', 'SBS16', 'SBS19', 'SBS22',
                                                                           'SBS29', 'SBS30'))

diagDM_newsigs[[paste0(ct, '_3')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_3')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_3')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_3')]]$Y)))

## removing SBS6 and SBS40
new_sigs[[paste0(ct, '_5')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c('SBS1', 'SBS4', 'SBS5', 'SBS12', 'SBS16', 'SBS19', 'SBS22', 'SBS29', 'SBS30'))

rbind(compare_signaturefit_to_data(read_info_list[[ct]]$dataset_active_sigs, read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                             signature_defs = sigs_cosmic0),
      compare_signaturefit_to_data(new_sigs[[paste0(ct, '_5')]], read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                   signature_defs = sigs_cosmic0))

par(mfrow=c(1,2))
compare_signaturefit_to_data(read_info_list[[ct]]$dataset_active_sigs, read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                             signature_defs = sigs_cosmic0, plot = T)
compare_signaturefit_to_data(new_sigs[[paste0(ct, '_5')]], read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                             signature_defs = sigs_cosmic0, plot = T)


fullDM_newsigs[[paste0(ct, '_5')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_5')]],
                                                      model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = fullDM_newsigs[[paste0(ct, '_5')]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_subset5', '.RDS'))

diagDM_newsigs[[paste0(ct, '_5')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_5')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_5')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)),
           sort_by_slope = T, title = paste0(ct, ' (subset of signatures)'))
ggsave(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_betas_subset_no_SBS40.pdf"),
       onefile=FALSE, height = 3, width = 5)

tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_betas_allsigsdiagREDM.tex"),
                 height = 1.8, width = 5.5)
plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, names_cats =  vector_cats_to_logR(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)),
           sort_by_slope = T, title = paste0(ct, ' (subset)'))
dev.off()

tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_betas_subset_no_SBS40.tex"),
                 height = 1.8, width = 3)
plot_betas(diagDM_newsigs[[paste0(ct, '_5')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)),
           sort_by_slope = T, title = paste0(ct, ' (subset)'))
dev.off()


tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_betasfullREDMDL_subset_no_SBS40.tex"),
                 height = 1.8, width = 3)
plot_betas(fullDM_newsigs[[paste0(ct, '_5')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)),
           sort_by_slope = T, title = paste0(ct, ' (subset)'))
dev.off()

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_barplot_subset_no_SBS40.pdf"),
       onefile=FALSE, height = 3, width = 5)
give_barplot_from_obj(new_sigs[[paste0(ct, '_5')]], only_normalised = T, levels_signatures=sigs_cosmic, nrow_plot = 1,
                      title = '', title_facets = c(NA, NA, 'Clonal', 'Subclonal'))
dev.off()

tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_barplot_subset_no_SBS40.tex"),
                 height = 1.8, width = 3)
give_barplot_from_obj(new_sigs[[paste0(ct, '_5')]], only_normalised = T, levels_signatures=sigs_cosmic, nrow_plot = 1,
                      title = '', title_facets = c(NA, NA, 'Clonal', 'Subclonal'))
dev.off()


compare_betas_tmb(tmb_obj_1 = diagDM_newsigs[[paste0(ct, '_4')]],
                  names_cats1 = colnames(new_sigs[[paste0(ct, '_4')]]$Y),
                  tmb_obj_2 = diagDM_newsigs[[paste0(ct, '_5')]],
                  names_cats2 = colnames(new_sigs[[paste0(ct, '_5')]]$Y),
                  names_groups = c('No SBS40', 'No SBS40, no SBS6'))

compare_betas_tmb(tmb_obj_1 = diagDM_newsigs[[paste0(ct, '_3')]],
                  names_cats1 = colnames(new_sigs[[paste0(ct, '_3')]]$Y),
                  tmb_obj_2 = diagDM_newsigs[[paste0(ct, '_5')]],
                  names_cats2 = colnames(new_sigs[[paste0(ct, '_5')]]$Y),
                  names_groups = c('No SBS40, no SBS5', 'No SBS40, no SBS6'))

compare_betas_tmb(tmb_obj_1 = read_info_list[[ct]]$diagRE_DMDL_SP,
                  names_cats1 = colnames(read_info_list[[ct]]$dataset_active_sigs$Y),
                  tmb_obj_2 = fullDM_newsigs[[paste0(ct, '_5')]],
                  names_cats2 = colnames(new_sigs[[paste0(ct, '_5')]]$Y),
                  names_groups = c('Original signatures', 'Subset'), include_missing_as_inf = T)+
  guides(alpha='none')

##-----------------------------------------------------------------------------------------------------##

ct <- "Lung-SCC" 

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS1', 'SBS5', 'SBS8')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS1', 'SBS5', 'SBS8'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)


plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)))

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS1', 'SBS4'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS4'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../"),
  give_plot_bleeding(names_sigs = c('SBS1', 'SBS5'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../"), ncol=2),
  nrow=2, rel_heights = c(2, 1)))
dev.off()

# par(mfrow=c(1,1))
# for(i in 1:nrow(sigs_cosmic)){
#   plot(as.numeric(cbind(sigs_cosmic[,c(1,5)], sigs_cosmic[,-c(1,5)])[i,]), type='l', main=i)
# }
# 
# as.numeric(cbind(sigs_cosmic[,c(1,5)], sigs_cosmic[,-c(1,5)])[39,]) ## high in SBS1 (but also quite high in another)
# as.numeric(cbind(sigs_cosmic[,c(1,5)], sigs_cosmic[,-c(1,5)])[35,]) ## high in SBS1 (but also quite high in another)
# 
# rownames(sigs_cosmic) <- gsub('>', '/', rownames(sigs_cosmic))
# 
# a <- read.table("../../data/restricted/pcawg/pcawg_restricted_snv/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.vcf.gz", sep='\t', header = T)
# a <- read.table("../../data/restricted/pcawg/pcawg_restricted_snv_counts/0009b464-b376-4fbc-8a56-da538269a02f", sep='\t', header = T)
# a <- read.table("../../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/0009b464-b376-4fbc-8a56-da538269a02f_cluster_assignments.txt.gz", sep='\t', header = T)
# a <- read.table("../../data/restricted/pcawg/", sep='\t', header = T)
# head(a)
# 
# a <- readRDS("../../data/roo/Biliary-AdenoCA_nucleotidesubstitution3_ROO.RDS")
# 
# plot(a[a$mutation == rownames(sigs_cosmic)[35],'VAF'])


##-----------------------------------------------------------------------------------------------------##

ct <- "Lymph-BNHL" 

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS13', 'SBS5', 'SBS56'),
                                                                                      c('SBS1', 'SBS40')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS13', 'SBS5', 'SBS56', 'SBS1', 'SBS40'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS56', 'SBS5'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS56', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS13', 'SBS5'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS13', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"), ncol=2),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Subset of sigs'),
                     nrow=1),
  nrow=3, rel_heights = c(2, 1, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Lymph-CLL"


##-----------------------------------------------------------------------------------------------------##

ct <- "Ovary-AdenoCA"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 10000)))
hypermutv2 <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 25000)))

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS3', 'SBS26'),
                                                                                      c('SBS41', 'SBS40')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS3', 'SBS26', 'SBS40', 'SBS41'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]

diagRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS3', 'SBS26'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS26'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS3', 'SBS41'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS41'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"), ncol=2),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS40', 'SBS26'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS40', 'SBS26'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS40', 'SBS41'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS40', 'SBS41'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"), ncol=2),
  give_plot_bleeding(names_sigs = c('SBS5', 'SBS3', 'SBS41'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS3', 'SBS41'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../"),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Amalgamation of sigs'),
                     nrow=1),
  nrow=5, rel_heights = c(1.5,0.5, 0.5, 0.5, 1)))
dev.off()


plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, vector_cats_to_logR(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)),
           sort_by_slope = T)

## no SBS40
new_sigs[[paste0(ct, '_1')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8",  "SBS13", "SBS18", "SBS26", "SBS35", "SBS39", "SBS41"))

diagDM_newsigs[[paste0(ct, '_1')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_1')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_1')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_1')]]$Y)), sort_by_slope = T)

## no SBS40, SBS41
new_sigs[[paste0(ct, '_2')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8",  "SBS13", "SBS18", "SBS26", "SBS35", "SBS39"))

diagDM_newsigs[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_2')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_2')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_2')]]$Y)), sort_by_slope = T)

## no SBS26, SBS40, SBS41
new_sigs[[paste0(ct, '_3')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8",  "SBS13", "SBS18", "SBS35", "SBS39"))

diagDM_newsigs[[paste0(ct, '_3')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_3')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_3')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_3')]]$Y)), sort_by_slope = T)

## no SBS26, SBS35, SBS40, SBS41
new_sigs[[paste0(ct, '_4')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8",  "SBS13", "SBS18", "SBS39"))

diagDM_newsigs[[paste0(ct, '_4')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_4')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_4')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_4')]]$Y)), sort_by_slope = T)

## no SBS26, SBS35, SBS39, SBS40, SBS41
new_sigs[[paste0(ct, '_5')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8",  "SBS13", "SBS18"))

diagDM_newsigs[[paste0(ct, '_5')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_5')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
fullDM_newsigs[[paste0(ct, '_5')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_5')]],
                                                      model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

plot_betas(diagDM_newsigs[[paste0(ct, '_5')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)), sort_by_slope = T)
plot_betas(fullDM_newsigs[[paste0(ct, '_5')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)), sort_by_slope = T)

## no SBS18, SBS26, SBS35, SBS39, SBS40, SBS41
new_sigs[[paste0(ct, '_6')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8",  "SBS13"))

diagDM_newsigs[[paste0(ct, '_6')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_6')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
fullDM_newsigs[[paste0(ct, '_6')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_6')]],
                                                      model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

## no SBS8, SBS18, SBS26, SBS35, SBS39, SBS40, SBS41
new_sigs[[paste0(ct, '_7')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS13"))

diagDM_newsigs[[paste0(ct, '_7')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_7')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
fullDM_newsigs[[paste0(ct, '_7')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_7')]],
                                                      model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

## different baseline
give_all_baselines <- function(dataset_arg){
  lapply(1:ncol(dataset_arg$Y), function(i){
    if(i == 1){
      .xx <- c((i+1):ncol(dataset_arg$Y) ,i )
    }else if(i == ncol(dataset_arg$Y)){
      .xx <- c(1:(i-1), i )
    }else{
      .xx <- c(1:(i-1), (i+1):ncol(dataset_arg$Y) ,i )
    }
    resort_columns(dataset_arg, .xx)
  })
}

all_datasets_ov_subset <- give_all_baselines(dataset_arg = new_sigs[[paste0(ct, '_7')]])
all_datasets_ov_subset_TMB <-lapply(all_datasets_ov_subset, function(i){
  wrapper_run_TMB(object = i, model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
  })
sapply(all_datasets_ov_subset_TMB, wald_TMB_wrapper) ## all consistent


# no sbs18
plot_betas(diagDM_newsigs[[paste0(ct, '_6')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_6')]]$Y)), sort_by_slope = T)
plot_betas(fullDM_newsigs[[paste0(ct, '_6')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_6')]]$Y)), sort_by_slope = T)

## no sbs8
plot_betas(diagDM_newsigs[[paste0(ct, '_7')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_7')]]$Y)), sort_by_slope = T)
plot_betas(fullDM_newsigs[[paste0(ct, '_7')]], names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_7')]]$Y)), sort_by_slope = T)

saveRDS(fullDM_newsigs[[paste0(ct, '_5')]],
        file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_subset5', '.RDS'))
saveRDS(diagDM_newsigs[[paste0(ct, '_5')]],
        file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_subset5', '.RDS'))


# pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_betas_subset.pdf"), onefile=FALSE, height = 3, width = 5)
tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_betas_subset.tex"), onefile=FALSE, height = 2, width = 5.5)
grid.arrange(plot_betas(diagDM_newsigs[[paste0(ct, '_5')]],
                        names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)),
                        sort_by_slope = T, title = 'diagREDMDL'),
             plot_betas(fullDM_newsigs[[paste0(ct, '_5')]],
                        names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)),
                        sort_by_slope = T, title = 'fullREDMDL'), nrow=1)
dev.off()

tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_betas_subset_fullREDM.tex"), onefile=FALSE, height = 1.8, width = 2.8)
plot_betas(fullDM_newsigs[[paste0(ct, '_5')]],
          names_cats =  vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)),
          sort_by_slope = T, title = NULL)
dev.off()

tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_betas_subset_allsigsdiagREDM.tex"), onefile=FALSE, height = 1.8, width = 4)
plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP,
           names_cats =  vector_cats_to_logR(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)),
           sort_by_slope = T, title = NULL)
dev.off()


compare_signaturefit_to_data(read_info_list[[ct]]$dataset_active_sigs,
                             read_info_list[[ct]]$dataset_nucleotidesubstitution3, sigs_cosmic0)
compare_signaturefit_to_data(new_sigs[[paste0(ct, paste0('_5'))]],
                             read_info_list[[ct]]$dataset_nucleotidesubstitution3, sigs_cosmic0)
compare_signaturefit_to_data(new_sigs[[paste0(ct, paste0('_6'))]],
                             read_info_list[[ct]]$dataset_nucleotidesubstitution3, sigs_cosmic0)
compare_signaturefit_to_data(new_sigs[[paste0(ct, paste0('_7'))]],
                             read_info_list[[ct]]$dataset_nucleotidesubstitution3, sigs_cosmic0)

give_barplot_from_obj(new_sigs[[paste0(ct, '_5')]], levels_signatures=sigs_cosmic, legend_on = T, legend_bottom = T)

tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_barplot_subset.tex"),
                 height = 1.8, width = 2.5)
give_barplot_from_obj(new_sigs[[paste0(ct, '_5')]], only_normalised = T, levels_signatures=sigs_cosmic, nrow_plot = 1,
                      title = '', title_facets = c(NA, NA, 'Cl', 'Scl'))
dev.off()

tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_barplot_allsigs.tex"),
                 height = 1.8, width = 2.5)
give_barplot_from_obj(read_info_list[[ct]]$dataset_active_sigs, only_normalised = T, levels_signatures=sigs_cosmic, nrow_plot = 1,
                      title = '', title_facets = c(NA, NA, 'Cl', 'Scl'))
dev.off()

tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_changes_in_beta.tex"),
                 height = 2.5, width = 4.5)
compare_betas_tmb(tmb_obj_1 = read_info_list[[ct]]$diagRE_DMDL_SP,
                  names_cats1 = colnames(read_info_list[[ct]]$dataset_active_sigs$Y),
                  tmb_obj_2 = diagDM_newsigs[[paste0(ct, '_5')]],
                  names_cats2 = colnames(new_sigs[[paste0(ct, '_5')]]$Y),
                  names_groups = c('Original signatures', 'Subset'), include_missing_as_inf = T)+
  guides(alpha='none')
dev.off()

compare_betas_tmb(tmb_obj_1 = fullDM_newsigs[[paste0(ct, '_5')]],
                  names_cats1 = colnames(new_sigs[[paste0(ct, '_5')]]$Y),
                  tmb_obj_2 = diagDM_newsigs[[paste0(ct, '_5')]],
                  names_cats2 = colnames(new_sigs[[paste0(ct, '_5')]]$Y),
                  names_groups = c('Original signatures', 'Subset'), include_missing_as_inf = T)+
  guides(alpha='none')

plot_lambdas(diagDM_newsigs[[paste0(ct, '_2')]])
plot_lambdas(diagDM_newsigs[[paste0(ct, '_3')]])
plot_lambdas(diagDM_newsigs[[paste0(ct, '_4')]])
plot_lambdas(fullDM_newsigs[[paste0(ct, '_5')]])

hypermut

dataset_extra[[paste0(ct, '_nucleotidesubstitution1nonhypermutated')]] <- give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_nucleotidesubstitution1, samples_to_remove = hypermutv2)
dim(dataset_extra[[paste0(ct, '_nucleotidesubstitution1nonhypermutated')]]$Y)
dim(read_info_list[[ct]]$dataset_nucleotidesubstitution1$Y)
fullDM_newsigs[[paste0(ct, '_nucleotidesubstitution1nonhypermutated')]] <- wrapper_run_TMB(object = dataset_extra[[paste0(ct, '_nucleotidesubstitution1nonhypermutated')]],
                                                                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

## the hypermutated sample is not a problem
grid.arrange(plot_betas(nucleotide1[[ct]], names_cats = names_trinucleotide, title = 'all samples'),
             plot_betas(fullDM_newsigs[[paste0(ct, '_nucleotidesubstitution1nonhypermutated')]], names_cats = names_trinucleotide, title = 'no hypermutated sample'))
give_barplot_from_obj(read_info_list[[ct]]$dataset_nucleotidesubstitution1, only_normalised = F, nrow_plot = 1, scale_color_manual_vec = nucleotide_colours)

tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_nucleotidebetas.tex"),
                 height = 2, width = 3)
plot_betas(nucleotide1[[ct]], names_cats = gsub(">", "$>$", names_trinucleotide), title = NULL)
dev.off()

plot_lambdas(nucleotide1[[ct]])


tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_nucleotidebarplot.tex"),
                 height = 2, width = 3)
give_barplot_from_obj(rename_Y_dollar(read_info_list[[ct]]$dataset_nucleotidesubstitution1), only_normalised = T, nrow_plot = 1, legend_on = F,
                      legend_bottom = T, scale_color_manual_vec = nucleotide_colours_dollar, title_facets = c(NA, NA, 'Clonal', 'Subclonal'))
dev.off()

for(ct in enough_samples){
  cat(ct, '\n')
  tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_sigma.tex"),
                   height = 3.5, width = 3.5)
  plot_covariance_mat(nucleotide1[[ct]], names_cats =  gsub(">", "$>$", names_trinucleotide), title = '',
                            lims=c(-1, 1))
  dev.off()
}
for(ct in enough_samples){
  cat(ct, '\n')
  tikzDevice::tikz(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_sigmanocluster.tex"),
                   height = 3.5, width = 3.5)
  plot_covariance_mat(nucleotide1[[ct]], names_cats =  gsub(">", "$>$", names_trinucleotide), title = '',
                      lims=c(-1, 1), arg_cluster_rows = F, arg_cluster_cols = F)
  dev.off()
}
ct <- "Ovary-AdenoCA"


## adding MMR signatures: SBS6 and SBS20, as in theory MMR is important in ov cancer, but were not included in PCAWG
new_sigs[[paste0(ct, '_newsigs1')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS13", "SBS18", "SBS6", "SBS20"))
diagDM_newsigs[[paste0(ct, '_newsigs1')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_newsigs1')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
fullDM_newsigs[[paste0(ct, '_newsigs1')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_newsigs1')]],
                                                      model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

give_barplot_from_obj(new_sigs[[paste0(ct, '_newsigs1')]])
plot_betas(diagDM_newsigs[[paste0(ct, '_newsigs1')]], names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_newsigs1')]]$Y)), sort_by_slope = T)

## only adding SBS6, not SBS20 (as there is none in ovarian cancer)
new_sigs[[paste0(ct, '_newsigs2')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                            subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5", 'SBS6', "SBS13", "SBS18"))

give_barplot_from_obj(new_sigs[[paste0(ct, '_newsigs2')]])
diagDM_newsigs[[paste0(ct, '_newsigs2')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_newsigs2')]],
                                                             model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
### which good, below!
plot_betas(diagDM_newsigs[[paste0(ct, '_newsigs2')]], names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_newsigs2')]]$Y)), sort_by_slope = T)

## adding SBS8
new_sigs[[paste0(ct, '_newsigs3')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                            subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5", "SBS8", "SBS13", "SBS18", "SBS6", "SBS20"))
diagDM_newsigs[[paste0(ct, '_newsigs3')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_newsigs3')]],
                                                             model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_newsigs3')]], sort_by_slope = T,
           names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_newsigs3')]]$Y)))

## adding SBS8, no 20
new_sigs[[paste0(ct, '_newsigs4')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                            subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5", "SBS8", "SBS13", "SBS18", "SBS6"))
diagDM_newsigs[[paste0(ct, '_newsigs4')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_newsigs4')]],
                                                             model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
fullDM_newsigs[[paste0(ct, '_newsigs4')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_newsigs4')]],
                                                             model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

plot_betas(diagDM_newsigs[[paste0(ct, '_newsigs4')]], sort_by_slope = T,
           names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_newsigs4')]]$Y)))
plot_betas(fullDM_newsigs[[paste0(ct, '_newsigs4')]], sort_by_slope = T,
           names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_newsigs4')]]$Y)))

plot_covariance_mat(fullDM_newsigs[[paste0(ct, '_newsigs4')]], names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_newsigs4')]]$Y)))

plot_covariance_mat(fullDM_newsigs[[paste0(ct, '_5')]],
                    names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5')]]$Y)))
fullDM_newsigs[[paste0(ct, '_5')]]

new_sigs[[paste0(ct, '_5wSBS36')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                     subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8",  "SBS13", "SBS36", "SBS18"))

diagDM_newsigs[[paste0(ct, '_5wSBS36')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_5wSBS36')]],
                                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_5wSBS36')]], sort_by_slope = T,
           names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5wSBS36')]]$Y)))


new_sigs[[paste0(ct, '_5wSBS36noSBS13')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                           subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8", "SBS36", "SBS18"))

diagDM_newsigs[[paste0(ct, '_5wSBS36noSBS13')]] <- wrapper_run_TMB(object = new_sigs[[paste0(ct, '_5wSBS36noSBS13')]],
                                                            model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
plot_betas(diagDM_newsigs[[paste0(ct, '_5wSBS36noSBS13')]], sort_by_slope = T,
           names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5wSBS36noSBS13')]]$Y)))

grid.arrange(plot_betas(diagDM_newsigs[[paste0(ct, '_5wSBS36')]], sort_by_slope = T,
                        names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5wSBS36')]]$Y))),
plot_betas(diagDM_newsigs[[paste0(ct, '_5wSBS36noSBS13')]], sort_by_slope = T,
           names_cats = vector_cats_to_logR(colnames(new_sigs[[paste0(ct, '_5wSBS36noSBS13')]]$Y)))
)
wald_TMB_wrapper(diagDM_newsigs[[paste0(ct, '_5wSBS36')]])
wald_TMB_wrapper(diagDM_newsigs[[paste0(ct, '_5wSBS36noSBS13')]])


## https://media.springernature.com/lw685/springer-static/image/art%3A10.1186%2Fs12967-022-03259-0/MediaObjects/12967_2022_3259_Fig3_HTML.png?as=webp
new_sigs[[paste0(ct, '_trying_to_stratify')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                                  subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8", "SBS13", "SBS24", "SBS30", "SBS32", "SBS36", "SBS35", "SBS39", "SBS41"))
createBarplot(normalise_rw(make_rownames_unique(new_sigs[[paste0(ct, '_trying_to_stratify')]]$Y)), order_labels = rownames(make_rownames_unique(new_sigs[[paste0(ct, '_trying_to_stratify')]]$Y))[order(normalise_rw(make_rownames_unique(new_sigs[[paste0(ct, '_trying_to_stratify')]]$Y))[,'SBS3'])])


rownames(new_sigs[[paste0(ct, '_5')]]$Y)[duplicated(rownames(new_sigs[[paste0(ct, '_5')]]$Y))] <- paste0(rownames(new_sigs[[paste0(ct, '_5')]]$Y)[duplicated(rownames(new_sigs[[paste0(ct, '_5')]]$Y))], '_2')
rownames(read_info_list[[ct]]$dataset_active_sigs$Y)[duplicated(rownames(read_info_list[[ct]]$dataset_active_sigs$Y))] <- paste0(rownames(read_info_list[[ct]]$dataset_active_sigs$Y)[duplicated(rownames(read_info_list[[ct]]$dataset_active_sigs$Y))], '_2')

hclust_ov <- hclust(dist(as(compositions::clr(normalise_rw(new_sigs[[paste0(ct, '_5')]]$Y)), 'matrix')))
hclust_ov$labels_old <- hclust_ov$labels
hclust_ov$labels = ifelse(grepl("_2", hclust_ov$labels), '*', '.')
plot(hclust_ov)

## with COSMIC v2
## here SBS3 does not show a clear stratification at all either
new_sigs[[paste0(ct, '_trying_to_stratifyv2')]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                                                      subset_signatures = c("SBS1",  "SBS2",  "SBS3",  "SBS5",  "SBS8", "SBS13", "SBS24", "SBS30"), signature_version = 'v2')
createBarplot(normalise_rw(make_rownames_unique(new_sigs[[paste0(ct, '_trying_to_stratifyv2')]]$Y)), order_labels = rownames(make_rownames_unique(new_sigs[[paste0(ct, '_trying_to_stratifyv2')]]$Y))[order(normalise_rw(make_rownames_unique(new_sigs[[paste0(ct, '_trying_to_stratifyv2')]]$Y))[,'SBS3'])])



createBarplot(normalise_rw(new_sigs[[paste0(ct, '_5')]]$Y[hclust_ov$labels_old,]))
createBarplot(normalise_rw(new_sigs[[paste0(ct, '_5')]]$Y), order_labels = rownames(new_sigs[[paste0(ct, '_5')]]$Y)[order(normalise_rw(new_sigs[[paste0(ct, '_5')]]$Y)[,'SBS3'])])

createBarplot(normalise_rw(read_info_list[[ct]]$dataset_active_sigs$Y), order_labels = rownames(read_info_list[[ct]]$dataset_active_sigs$Y)[order(normalise_rw(read_info_list[[ct]]$dataset_active_sigs$Y)[,'SBS3'])])
createBarplot(normalise_rw(read_info_list[[ct]]$dataset_active_sigs_MSE$Y), order_labels = rownames(read_info_list[[ct]]$dataset_active_sigs_MSE$Y)[order(normalise_rw(read_info_list[[ct]]$dataset_active_sigs_MSE$Y)[,'SBS3'])])

new_sigs_5_ov_stratify <- split_matrix_in_half(normalise_rw(new_sigs[[paste0(ct, '_5')]]$Y))
pheatmap::pheatmap(new_sigs_5_ov_stratify[[2]]-new_sigs_5_ov_stratify[[1]])
pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_new_sigs_5_ov_stratify.pdf"))
print(pheatmap::pheatmap(as(compositions::clr(new_sigs_5_ov_stratify[[2]]), 'matrix')-as(compositions::clr(new_sigs_5_ov_stratify[[1]]), 'matrix')))
dev.off()

singleRE_DM[[ct]] <- try({wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "singleRE_DM", use_nlminb=T, smart_init_vals=F)})
grid.arrange(
plot_betas(diagRE_DMDL_extra[[ct]]),
plot_betas(singleRE_DM[[ct]]))

## very interesting (two groups)
par(mfrow=c(1,1))
plot(density(  singleRE_DM[[ct]]$par.random))
plot(density(  diagRE_DMDL_extra[[ct]]$par.random))
## no such clustering here
par(mfrow=c(2,5))
for(i in 1:(length(python_like_select_name(diagRE_DMDL_extra[[ct]]$par.fixed, 'beta'))/2)){
  plot(density(matrix( diagRE_DMDL_extra[[ct]]$par.random, ncol=length(python_like_select_name(diagRE_DMDL_extra[[ct]]$par.fixed, 'beta'))/2)[,i]))
}

pairs(matrix( diagRE_DMDL_extra[[ct]]$par.random, ncol=length(python_like_select_name(diagRE_DMDL_extra[[ct]]$par.fixed, 'beta'))/2))

## PCA of random effects
plot(prcomp(matrix( diagRE_DMDL_extra[[ct]]$par.random, ncol=length(python_like_select_name(diagRE_DMDL_extra[[ct]]$par.fixed, 'beta'))/2))$x[,1:2])


##-----------------------------------------------------------------------------------------------------##

ct <- "Panc-AdenoCA"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 21000)))

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS3', 'SBS5', 'SBS17a', 'SBS20')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS3', 'SBS5', 'SBS17a', 'SBS20'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

fullRE_DMDL_extra[[ct]]
saveRDS(object = fullRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))
# fullRE_DMDL_extra[[ct]] = readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))

amalgamated_extra[[paste0(ct, '_2')]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS3', 'SBS5')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS3', 'SBS5'))]))))
diagRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_2')]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = diagRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_nonhypermutated', '.RDS'))

amalgamated_extra[[paste0(ct, '_3')]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                                            list_groupings = c(list(c('SBS3', 'SBS5'),
                                                                                                    c('SBS8', 'SBS13')),
                                                                                               as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS3', 'SBS5', 'SBS8', 'SBS13'))]))))
diagRE_DMDL_extra[[paste0(ct, '_3')]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_3')]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = diagRE_DMDL_extra[[paste0(ct, '_3')]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, '_3diagRE_DMDL_nonhypermutated', '.RDS'))


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Subset of sigs'),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS3', 'SBS40'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  give_plot_bleeding(names_sigs = c('SBS3', 'SBS5'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  give_plot_bleeding(names_sigs = c('SBS17a', 'SBS5'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS17a', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  give_plot_bleeding(names_sigs = c('SBS17a', 'SBS40'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS17a', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  rel_widths = c(1,1,1, 1), ncol=4),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS20', 'SBS5'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS20', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  give_plot_bleeding(names_sigs = c('SBS20', 'SBS40'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS20', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
                     give_plot_bleeding_wrapper( c('SBS51', 'SBS28')), ncol=3),
  nrow=4, rel_heights = c(0.5,.8,0.4, 0.4 )))
dev.off()

mean(normalise_rw(read_info_list[[ct]]$dataset_active_sigs$Y)[,'SBS3'])
mean(normalise_rw(read_info_list[["Panc-AdenoCA"]]$dataset_active_sigs$Y)[,'SBS3'])

mean(normalise_rw(read_info_list[[ct]]$dataset_active_sigs$Y)[,'SBS13'])
mean(normalise_rw(read_info_list[["Panc-AdenoCA"]]$dataset_active_sigs$Y)[,'SBS13'])


all_bleeding[(all_bleeding$ct == "Panc-AdenoCA") & (all_bleeding$sig == 'SBS3'),]

head(new_sigs[[paste0(ct, '_5')]]$Y)

normalised_by_s1s5_ovary <- sweep(new_sigs[[paste0(ct, '_5')]]$Y, 1, new_sigs[[paste0(ct, '_5')]]$Y[,'SBS1']+new_sigs[[paste0(ct, '_5')]]$Y[,'SBS5'], '/')
normalised_by_s1s5_ovary <- split_matrix_in_half(normalised_by_s1s5_ovary)
all(rownames(normalised_by_s1s5_ovary[[1]]) == rownames(normalised_by_s1s5_ovary[[2]])) ## should all be true

inf_to_large <- function(m, num_high=4){
  m[(is.infinite(m) & m>0)] <- num_high
  m[(is.infinite(m) & m<0)] <- -num_high
  m
}

ratios_normalised_by_s1s5_ovary <- normalised_by_s1s5_ovary[[2]]/normalised_by_s1s5_ovary[[1]]

pheatmap::pheatmap(inf_to_large(log(normalised_by_s1s5_ovary[[2]]/normalised_by_s1s5_ovary[[1]])), cluster_rows = F, cluster_cols = F)

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_changes_in_ALR_patient_stratification.pdf"), onefile=FALSE, height = 12)
print(pheatmap::pheatmap(inf_to_large(log(normalised_by_s1s5_ovary[[2]]/normalised_by_s1s5_ovary[[1]]))))
dev.off()

plot(density(remove_na(ratios_normalised_by_s1s5_ovary[,'SBS2']))); abline(v = 1, col='blue', lty='dashed')
plot(density(remove_na(ratios_normalised_by_s1s5_ovary[,'SBS13']))); abline(v = 1, col='blue', lty='dashed')

pheatmap::pheatmap(apply(ratios_normalised_by_s1s5_ovary > 0, 2, as.numeric))
pheatmap::pheatmap(apply(ratios_normalised_by_s1s5_ovary > 0, 2, as.numeric))

table(ratios_normalised_by_s1s5_ovary[,'SBS13'] > ratios_normalised_by_s1s5_ovary[,'SBS2'])
table(ratios_normalised_by_s1s5_ovary[,'SBS13'] > ratios_normalised_by_s1s5_ovary[,'SBS8'])
table(ratios_normalised_by_s1s5_ovary[,'SBS13'] > ratios_normalised_by_s1s5_ovary[,'SBS3'])

## does SBS3 only increase in a fraction of samples?
plot(density(remove_na(ratios_normalised_by_s1s5_ovary[,'SBS3']))); abline(v = 1, col='blue', lty='dashed')

##-----------------------------------------------------------------------------------------------------##

ct <- "Panc-Endocrine"


pancreas <- list(PancEndocrine=betas_df_to_softmax(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP, return_df = T) %>% dplyr::filter(type_beta == 'Slope')),
                 PancAdenoCA = betas_df_to_softmax(plot_betas(read_info_list[["Panc-AdenoCA"]]$diagRE_DMDL_SP, read_info_list[["Panc-AdenoCA"]]$logR_notsorted_SP, return_df = T) %>% dplyr::filter(type_beta == 'Slope')))

pancreas$PancAdenoCA <- pancreas$PancAdenoCA[match(names(pancreas$PancEndocrine), names((pancreas$PancAdenoCA)))]

pancreas2 <- list(PancEndocrine=betas_df_to_softmax(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP, return_df = T) %>% dplyr::filter(type_beta == 'Slope')),
                 PancAdenoCA = betas_df_to_softmax(plot_betas(diagRE_DMDL_extra[["Panc-AdenoCA"]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0("Panc-AdenoCA", '_2')]]$Y)), return_df = T) %>% dplyr::filter(type_beta == 'Slope')))
names(pancreas2$PancAdenoCA)[names(pancreas2$PancAdenoCA) == '3+'] <- '5'

pancreas2$PancAdenoCA <- pancreas2$PancAdenoCA[match(names(pancreas2$PancEndocrine), names((pancreas2$PancAdenoCA)))]

pancreas3 <- list(PancEndocrine=betas_df_to_softmax(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP, return_df = T) %>% dplyr::filter(type_beta == 'Slope')),
                  PancAdenoCA = betas_df_to_softmax(plot_betas(diagRE_DMDL_extra[[paste0("Panc-AdenoCA", "_3")]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0("Panc-AdenoCA", "_3")]]$Y)), return_df = T) %>% dplyr::filter(type_beta == 'Slope')))
names(pancreas3$PancAdenoCA)[names(pancreas3$PancAdenoCA) == '3+'] <- '5'
names(pancreas3$PancAdenoCA)[names(pancreas3$PancAdenoCA) == '8+'] <- '13'

pancreas3$PancAdenoCA <- pancreas3$PancAdenoCA[match(names(pancreas3$PancEndocrine), names((pancreas3$PancAdenoCA)))]

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS5', 'SBS39')),
  give_plot_bleeding_wrapper( c('SBS5', 'SBS26')), ncol=2),
  cowplot::plot_grid(ggplot(cbind.data.frame(sig=names(pancreas$PancEndocrine),
                                             betaslopeest=do.call('cbind', pancreas),
                                             abundance=(colSums(read_info_list[[ct]]$dataset_active_sigs$Y))),
         aes(x=betaslopeest.PancEndocrine, y=betaslopeest.PancAdenoCA, label=sig))+
    geom_point(aes(size=abundance))+geom_label_repel()+theme_bw()+geom_smooth(method = lm)+
      scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+guides(size=F),
  ggplot(cbind.data.frame(sig=names(pancreas2$PancEndocrine), betaslopeest=do.call('cbind', pancreas2),
                          abundance=(colSums(read_info_list[[ct]]$dataset_active_sigs$Y))),
         aes(x=betaslopeest.PancEndocrine, y=betaslopeest.PancAdenoCA, label=sig))+
    geom_point(aes(size=abundance))+geom_label_repel(min.segment.length = 0.2, , max.overlaps = 100)+theme_bw()+geom_smooth(method = lm)+
    scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+guides(size=F),
  ggplot(cbind.data.frame(sig=names(pancreas3$PancEndocrine), betaslopeest=do.call('cbind', pancreas3),
                          abundance=(colSums(read_info_list[[ct]]$dataset_active_sigs$Y))),
         aes(x=betaslopeest.PancEndocrine, y=betaslopeest.PancAdenoCA, label=sig))+
    geom_point(aes(size=abundance))+geom_label_repel(min.segment.length = 0.2, max.overlaps = 100)+theme_bw()+geom_smooth(method = lm)+
    scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+guides(size=F), ncol=3),
  nrow=4, rel_heights = c(0.5,0.8, 0.5, 1)))
dev.off()


##-----------------------------------------------------------------------------------------------------##

ct <- "Prost-AdenoCA"

hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 8000)))
hypermut

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS40', 'SBS33', 'SBS37')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS40', 'SBS33', 'SBS37'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = fullRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
                                title='fullRE_DMDL Amalgamation of sigs', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  nrow=4, rel_heights = c(0.5,0.5, 1, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Skin-Melanoma.cutaneous"


amalgamated_extra[[ct]] <- give_amalgamated_exposures_TMBobj( read_info_list[[ct]]$dataset_active_sigs,
                                                              list_groupings = c(list(c('SBS5', 'SBS2', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d'),
                                                                                      c('SBS17a', 'SBS17b')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS5', 'SBS2', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS17a', 'SBS17b'))])))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 15)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40', 'SBS58')),
  # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40')),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a')),
                     give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7b')), ncol=2),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7c')),
                     give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7d')), ncol=2),
  give_plot_bleeding_wrapper( c('SBS1', 'SBS7c')),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
                                title='fullRE_DMDL Amalgamation of sigs', sort_by_slope = T),
                     nrow=1),
  nrow=7, rel_heights = c(0.5,1, .6, .6, 0.5, 0.5, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##
ct <- "Stomach-AdenoCA"

amalgamated_extra[[paste0(ct)]] <- (give_amalgamated_exposures_TMBobj(read_info_list[[ct]]$dataset_active_sigs,
                                                                            list_groupings = c(list(c('SBS20', 'SBS43', 'SBS44'),
                                                                                                    c('SBS26', 'SBS28', 'SBS9', 'SBS5')),
                                          as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS20', 'SBS43', 'SBS44', 'SBS26', 'SBS28', 'SBS5', 'SBS9'))]))))
diagRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct)]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = diagRE_DMDL_extra[[paste0(ct)]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, '_diagRE_DMDL_amalgamation', '.RDS'))


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  cowplot::plot_grid(plot_betas(diagRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct)]]$Y)),
                                title='Sorted sigs diagRE_DMDL_SP amalgamation', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40', 'SBS58')),
  # # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40')),
  # cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a')),
  #                    give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7b')), ncol=2),
  # cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7c')),
  #                    give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7d')), ncol=2),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS3', 'SBS9', 'SBS20', 'SBS26', 'SBS28', 'SBS43', 'SBS44', 'SBS51')),
                     give_plot_bleeding_wrapper( c('SBS3', 'SBS9', 'SBS20', 'SBS26', 'SBS28', 'SBS43', 'SBS5', 'SBS51')),nrow=1),
  give_plot_bleeding_wrapper(c("SBS1",   "SBS2",   "SBS3",   "SBS5",   "SBS9",   "SBS13",  "SBS15",  "SBS17a", "SBS17b", "SBS18",
                               "SBS20",  "SBS21",  "SBS26",  "SBS28",  "SBS40" , "SBS41",  "SBS43",  "SBS44",  "SBS51",  "SBS58")),
  # cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
  #                               title='fullRE_DMDL Amalgamation of sigs', sort_by_slope = T),
  #                    nrow=1),
  nrow=5, rel_heights = c(0.5,0.5, 0.8, 0.5, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Thy-AdenoCA"

hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 3000)))
hypermut

amalgamated_extra[[ct]] <- give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS2', 'SBS13', 'SBS58')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS2', 'SBS13', 'SBS58'))])))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
saveRDS(object = fullRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))



amalgamated_extra[[paste0(ct, '_2')]] <- give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                             list_groupings = c(list(c('SBS2', 'SBS13', 'SBS58', 'SBS40')),
                                                                                as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS2', 'SBS13', 'SBS58', 'SBS40'))])))
fullRE_DMDL_extra[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_2')]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[paste0(ct, '_2')]]

c

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 15)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40', 'SBS58')),
  # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40')),
  # cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a')),
  #                    give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7b')), ncol=2),
  # cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7c')),
  #                    give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7d')), ncol=2),
  # give_plot_bleeding_wrapper( c('SBS1', 'SBS7c')),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
                                title='fullRE_DMDL Amalgamation of sigs', sort_by_slope = T),
                     plot_betas(fullRE_DMDL_extra[[paste0(ct, '_2')]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_2')]]$Y)),
                                title='fullRE_DMDL Amalgamation 2 of sigs', sort_by_slope = T),
                     nrow=1),
  cowplot::plot_grid(plot_betas(diagRE_DMDL_extra[[paste0(ct, '_3')]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_3')]]$Y)),
                                title='diagRE_DMDL Amalgamation 3 of sigs', sort_by_slope = T),
                     plot_betas(fullRE_DMDL_extra[[paste0(ct, '_3')]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_3')]]$Y)),
                                title='fullRE_DMDL Amalgamation 3 of sigs', sort_by_slope = T),
                     nrow=1),
  nrow=4, rel_heights = c(0.5,1, .6, .6)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Uterus-AdenoCA"

hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_active_sigs$Y)) > 10000)))
hypermut

amalgamated_extra[[paste0(ct, '_nohypermut')]] <- give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut)

dim(read_info_list[[ct]]$dataset_active_sigs$Y)
dim(amalgamated_extra[["Uterus-AdenoCA_nohypermut"]]$Y)
diagRE_DMDL_extra[[paste0(ct)]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_nohypermut')]],
                                                         model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_extra[[paste0(ct)]]


amalgamated_extra[[paste0(ct, '_2')]] <- give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_active_sigs, samples_to_remove = hypermut),
                                                                           list_groupings = c(list(c('SBS5', 'SBS28', 'SBS44', 'SBS26'),
                                                                                                   c('SBS40', 'SBS3', 'SBS6', 'SBS15', 'SBS14'),
                                                                                                   c('SBS10a', 'SBS10b'),
                                                                                                   c('SBS2', 'SBS13')),
 as.list(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_active_sigs$Y) %in% c('SBS5', 'SBS28', 'SBS44', 'SBS26', 'SBS6', 'SBS14', 'SBS15', 'SBS10a', 'SBS10b', 'SBS40', 'SBS3', 'SBS2', 'SBS13'))])))
colnames(amalgamated_extra[[paste0(ct, '_2')]]$Y)
fullRE_DMDL_extra[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_2')]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[paste0(ct, '_2')]]

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 15)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  cowplot::plot_grid(plot_betas(diagRE_DMDL_extra[[paste0(ct)]],
                                vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_nohypermut')]]$Y)),
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     plot_betas(fullRE_DMDL_extra[[paste0(ct, '_2')]],
                                vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_2')]]$Y)),
                                title='Amalgamation fullRE_DMDL_extra', sort_by_slope = T),
                                        nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS28', 'SBS5', 'SBS40')),
                     give_plot_bleeding_wrapper( c('SBS3', 'SBS5', 'SBS40')), nrow=1),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS6', 'SBS5', 'SBS40')),
                     give_plot_bleeding_wrapper( c('SBS15', 'SBS6')), nrow=1),
  nrow=5, rel_heights = c(0.5,0.5, 1, .5, .5)))
dev.off()

##-----------------------------------------------------------------------------------------------------##
sbs2_sbs8 <- lapply(read_info_list, function(i) try(normalise_rw(i$dataset_active_sigs$Y)[,c('SBS2', 'SBS8')]))
sbs13_sbs8 <- lapply(read_info_list, function(i) try(normalise_rw(i$dataset_active_sigs$Y)[,c('SBS13', 'SBS8')]))
head(melt(sapply(sbs2_sbs8[sapply(sbs2_sbs8, typeof) == "double"], function(i) i[,1]/i[,2])))
par(mfrow=c(1,1))
df_sbs2_sbs8 <- data.frame(do.call('rbind', sbs2_sbs8[sapply(sbs2_sbs8, typeof) == "double"]))
df_sbs13_sbs8 <- data.frame(do.call('rbind', sbs13_sbs8[sapply(sbs13_sbs8, typeof) == "double"]))
df_sbs2_sbs8$ct <- rep(names(sbs2_sbs8[sapply(sbs2_sbs8, typeof) == "double"]), sapply(sbs2_sbs8[sapply(sbs2_sbs8, typeof) == "double"], nrow))
df_sbs2_sbs8$clonal = unlist(sapply(sapply(table(df_sbs2_sbs8$ct), function(i) i/2), function(j) rep(c("Clonal", 'Subclonal'), each=j)))
df_sbs13_sbs8$ct <- rep(names(sbs13_sbs8[sapply(sbs13_sbs8, typeof) == "double"]), sapply(sbs13_sbs8[sapply(sbs13_sbs8, typeof) == "double"], nrow))
df_sbs13_sbs8$clonal = unlist(sapply(sapply(table(df_sbs13_sbs8$ct), function(i) i/2), function(j) rep(c("Clonal", 'Subclonal'), each=j)))
cowplot::plot_grid(ggplot(df_sbs2_sbs8, aes(x=SBS2, y=SBS8, col=gsub("\\..*", "", ct)))+
                     geom_point(col='black', size=3, alpha=0.2)+
                     geom_point()+#scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
                     theme_bw()+
    facet_wrap(.~ct+clonal, nrow=2, scales = "free")+
    scale_color_manual(values = pcawg_palette)+
    geom_smooth(method = lm, col='black')+guides(col=F),
  ggplot(df_sbs13_sbs8, aes(x=SBS13, y=SBS8, col=gsub("\\..*", "", ct)))+
    geom_point(col='black', size=3, alpha=0.2)+
    geom_point()+#scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
    theme_bw()+
    facet_wrap(.~ct+clonal, nrow=2, scales = "free")+
  scale_color_manual(values = pcawg_palette)+
    geom_smooth(method = lm, col='black')+guides(col=F), nrow=2)
ggsave("../../results/results_TMB/pcawg/SBS8_APOBEC_correlation_raw.pdf", height = 5, width = 10)
##-----------------------------------------------------------------------------------------------------##
python_like_select(list.files("../../results/results_TMB/pcawg/reports_per_cancer_type"), ct)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Analysing the number of mutations in each group 
#'  there is an increase in variance from clonal to subclonal in the majority of cancer types.
#'  Can it be because there are fewer mutations in clonal than in subclonal, for most cancer types?
#'  Answer: no. whether the overdispersion is hihger in the first or second group does not depend on
#'  whether there are more mutations in one group or the other, since there are three ct (of the 23) 
#'  with more mutations in the second group than in the first: CNS-GBM, CNS-PiloAstro, ColoRect-AdenoCA,
#'   and of those only piloastro has a higher dispersion in the first group.

idx_ct <- 1

num_muts_groups <- lapply(read_info_list, function(read_info_list_ct)sapply(split_matrix_in_half(read_info_list_ct$dataset_active_sigs$Y), sum))
sapply(num_muts_groups, function(i) i[1] > i[2])
xtable::xtable((cbind(apply(do.call('rbind',num_muts_groups), 2, as.integer), sapply(num_muts_groups, function(i) i[1] > i[2]))))

num_muts_groups <- melt(num_muts_groups)
table(sapply(num_muts_groups, function(i) i[1] > i[2])) num_muts_groups$group <- c('Clonal', 'subclonal')

ggplot(num_muts_groups, aes(x=L1, y=value, group=group, col=group))+geom_bar(stat = "identity")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  

## IN PAPER CURRENTLY IT IS THE NONEXO-SP THAT IS BEING USED, BUT IS THAT THE CASE FOR BETAS, TOO? (fig 6, 7, 8)
# KRCC_DMDL <- wrapper_run_TMB(object = read_info_list$`Kidney-RCC.papillary`$dataset_active_sigs,
#                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F, return_opt=T)
# KRCC_DMDLworked <- KRCC_DMDL
# system("mv ../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_Kidney-RCC.papillary_signaturesPCAWG.RDS ../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_Kidney-RCC.papillary_signaturesPCAWG_prev.RDS")
# saveRDS(KRCC_DMDL$return_report, paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_", 
#                           "Kidney-RCC.papillary", "_signaturesPCAWG.RDS"))
# saveRDS(KRCC_DMDL$return_report, paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_", 
#                                         "Kidney-RCC.papillary", "_signaturesPCAWG_whole_nlminb.RDS"))

diagRE_DMDL <- lapply(read_info_list, function(i) i$diagRE_DMDL_SP)
diagRE_DMDL$`Kidney-RCC.papillary`
sapply(diagRE_DMDL, function(i) i$par.fixed)
ovrdisp <- do.call('rbind.data.frame', lapply(1:length(diagRE_DMDL), try(function(idx){
  if(diagRE_DMDL[[idx]]$pdHess){
    cbind.data.frame( plot_lambdas(diagRE_DMDL[[idx]], return_df=T, plot=F), ct=names(diagRE_DMDL)[idx])
  }else{
    c(NA, NA)
  }
})))
ovrdisp[ovrdisp$name == 'Lambda 1','name'] = 'Clonal'
ovrdisp[ovrdisp$name == 'Lambda 2','name'] = 'Subclonal'
# ovrdisp$differentially_abundant = ifelse(ovrdisp$ct %in% names(differential_precision[(differential_precision <= 0.05)]), yes = '*', no = '')
differential_precision_2 <- p.adjust(sapply(diagRE_DMDL, ttest_TMB_wrapper_overdisp), method = 'fdr')
names(differential_precision_2) <- names(diagRE_DMDL)

ovrdisp$differential_precision_2 = ifelse(ovrdisp$ct %in% names(differential_precision_2[(differential_precision_2 <= 0.05)]), yes = '*', no = '')
ggplot(ovrdisp, aes(x=ct,  y=`Estimate`, group=name, col=name))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=`Estimate`-`Std..Error`,
                    ymax=`Estimate`+`Std..Error`), width=.1, position=position_dodge(width=0.5))+
  theme_bw()+
  geom_text(aes(y=Inf, label=differential_precision_2, vjust=1.8), col='black')+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='', y='Log lambda', col='Group')+theme(legend.position = "bottom")+
  theme(        legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(-10,-10,-10,-10),
                plot.margin = unit(c(1,1,1,1), "cm"))+
  scale_color_manual(values=c('#3b4d61', '#ef9d10'))

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## How many clusters do the files have?
list_files_pcawg_subclonal <- list.files("../../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/", full.names = T)
list_files_pcawg_subclonal <- list_files_pcawg_subclonal[grep('cluster_assignments', list_files_pcawg_subclonal)]

num_clusters <- sapply(list_files_pcawg_subclonal, function(i) python_like_select(colnames(read.table(i, head=T)), 'cluster'))
saveRDS(num_clusters, "../../data/pcawg_robjects_cache/num_clusters_samples.RDS")

clusters_probs <- sapply(list_files_pcawg_subclonal, function(i) python_like_select_colnames(read.table(i, head=T), 'cluster'))
saveRDS(clusters_probs, "../../data/pcawg_robjects_cache/clusters_probs_samples.RDS")

num_clusters_df <- data.frame(num_clusters=sapply(num_clusters, length),
                              sample=gsub("_cluster_assignments.txt.gz", "", basename(names(num_clusters))))
rownames(num_clusters_df) <- num_clusters_df$sample
num_clusters_df$ct <- pcawg_meta$Project[match(num_clusters_df$sample, gsub("\\..*", "", pcawg_meta$File.Name))]

ggplot(num_clusters_df, aes(x=ct, y=num_clusters))+geom_boxplot()+geom_jitter()+theme_bw()+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
ggsave("../../results/exploratory/num_clusters_samples_per_cancer_type.pdf", width = 8, height = 3)

num_clusters_df$cthistological = pcawg_meta2$histology_detailed[match(num_clusters_df$sample,
                                                                      pcawg_meta2$samplename)]
ggplot(num_clusters_df, aes(x=cthistological, y=num_clusters))+geom_boxplot()+geom_jitter()+theme_bw()+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
ggsave("../../results/exploratory/num_clusters_samples_per_cancer_type_histological.pdf", width = 8, height = 3)

dim(num_clusters_df[is.na(num_clusters_df$cthistological),])

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
xtable::xtable(cbind(num_active=sapply(enough_samples, function(ct) ncol(read_info_list[[ct]]$dataset_active_sigs$Y)),
      num_nonexo=sapply(enough_samples, function(ct) ncol(read_info_list[[ct]]$dataset_nonexo$Y))))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Mean-variance relationship
meanvar <- lapply(enough_samples, function(ct){
  a <- data.frame(t(apply(read_info_list[[ct]]$dataset_active_sigs$Y, 2, function(i){
    cbind(mean=mean(i), var=var(i))
  })))
  colnames(a) <- c('mean', 'var')
  a
})
names(meanvar) <- enough_samples

meanvar_normalised <- lapply(enough_samples, function(ct){
  a <- data.frame(t(apply(normalise_rw(read_info_list[[ct]]$dataset_active_sigs$Y), 2, function(i){
    cbind(mean=mean(i), var=var(i))
  })))
  colnames(a) <- c('mean', 'var')
  a
})
names(meanvar_normalised) <- enough_samples

ggplot(melt(meanvar, id.vars=c('mean', 'var')), aes(x=mean, y=var))+
  geom_point()+facet_wrap(.~L1, scales = "free", ncol=4)+theme_bw()
ggsave("../../results/exploratory/mean_var_correlation_pcawg.pdf", width = 8)
ggplot(melt(meanvar_normalised, id.vars=c('mean', 'var')), aes(x=mean, y=var))+
  geom_point()+facet_wrap(.~L1, scales = "free", ncol=4)+theme_bw()
ggsave("../../results/exploratory/mean_var_correlation_pcawg_normalised.pdf", width = 8)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Stratification of samples

require(pheatmap)

sapply(enough_samples, function(ct){
  cat(ct, '\n')
  xxxx <- split_rows(read_info_list[[ct]]$dataset_active_sigs$Y, f=duplicated(rownames(read_info_list[[ct]]$dataset_active_sigs$Y)))
  stopifnot(rownames(xxxx[[1]]) == rownames(xxxx[[2]]))
  
  ph1 <- pheatmap::pheatmap(normalise_rw(xxxx[[1]]), show_rownames = F)
  ph2 <- pheatmap::pheatmap(normalise_rw(xxxx[[2]])[ph1$tree_row$order,], show_rownames = F, cluster_rows = F)
  
  pdf(paste0("../../results/exploratory/stratification_pcawg_samples/stratification_", gsub("[.]", "_", ct), ".pdf"),
      width = 6, height=3)
  print(cowplot::plot_grid(ph1[[4]], ph2[[4]], rel_widths = c(0.5,0.44), labels=ct))
  dev.off()
})

 ##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## LRT test instead of Wald test
diagDM_LRT_reduced <- lapply(enough_samples, function(ct){
  cat(ct, '\n')
  x <- wrapper_run_TMB(object = read_info_list[[ct]]$dataset_active_sigs,
                       model = "diagRE_DMreducedmodel", use_nlminb=T, smart_init_vals=F, return_opt=T)
  saveRDS(x, file = paste0("../../results/results_TMB/pcawgLRT/diagREDMreducedmodel_", ct, '.RDS'))
  x
})
names(diagDM_LRT_reduced) <- enough_samples

diagDM_LRT_full <- lapply(enough_samples, function(ct){
  cat(ct, '\n')
  x <- wrapper_run_TMB(object = read_info_list[[ct]]$dataset_active_sigs,
                       model = "diagRE_DM", use_nlminb=T, smart_init_vals=F, return_opt=T)
  saveRDS(x, file = paste0("../../results/results_TMB/pcawgLRT/diagREDMfullmodel_", ct, '.RDS'))
  x
})
names(diagDM_LRT_full) <- enough_samples

diagDM_LRT_reduced[[1]]$opt$objective
diagDM_LRT_full[[1]]$opt$objective

diagDM_LRT_reduced[[2]]$opt$objective
diagDM_LRT_full[[2]]$opt$objective

2*(diagDM_LRT_full[[1]]$opt$objective-diagDM_LRT_reduced[[1]]$opt$objective)

pvals_LRT <- sapply(enough_samples, function(ct){
  pchisq(q = 2*(diagDM_LRT_reduced[[ct]]$opt$objective-diagDM_LRT_full[[ct]]$opt$objective),
       df = ncol(read_info_list[[ct]]$dataset_active_sigs$Y)-1, lower.tail = F)
})

pvals_LRT

pvals_wald <- sapply(diagDM_LRT_full, function(i) wald_TMB_wrapper(i$return_report) )

pdf("../../results/results_TMB/LRT_Wald_pval_comparison.pdf", height = 4, width = 3.5)
plot(log(pvals_wald), log(pvals_LRT))
abline(coef = c(0,1))
dev.off()

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Comparing two extractions of signatures

read_info_list[[1]]$dataset_all_sigs
read_info_list[[1]]$dataset_active_sigs

MSE_ROO2_files <- list.files("../../data/roo/", full.names = T)[grepl('signaturesMSE_ROO2', list.files("../../data/roo/"))]
MSE_ROO2 <- sapply(MSE_ROO2_files, readRDS)
names(MSE_ROO2) <- gsub("_signaturesMSE_ROO2.RDS", "", basename(MSE_ROO2_files))

MSE_ROO2$`Bladder-TCC`@count_matrices_all

dim(MSE_ROO2$`Bone-Osteosarc`@count_matrices_active[[1]])

dim(read_info_list$`Bone-Osteosarc`$dataset_active_sigs$Y)

ct_it <- 'Bone-Osteosarc'

pdf("../../results/exploratory/mutsigexplorer_signature_exposures_comparison.pdf")
for(ct_it in names(MSE_ROO2)){
  
  if(ct_it %in% names(read_info_list)){
    mat1 <- do.call('rbind', MSE_ROO2[[ct_it]]@count_matrices_all)
    mat2 <- read_info_list[[ct_it]]$dataset_all_sigs$Y
    
    not_same_sigs <- (all(dim(mat1) == dim(mat2)))
    
    colnames(mat2)[!(colnames(mat2) %in% colnames(mat1))]
    colnames(mat1)[!(colnames(mat1) %in% colnames(mat2))]
    remove_na <- function(i) i[!is.na(i)]
    mat2 <- mat2[,remove_na(match(colnames(mat1), colnames(mat2)))]
    mat1 <- mat1[,remove_na(match(colnames(mat2), colnames(mat1)))]
    
    
    if(! not_same_sigs){
      if(nrow(mat1) != nrow(mat2)){
       print(ggplot())
        next
      }
      title=paste0('All sigs\n', ct_it, '\n Not the same number of sigs')
    }else{
      title <- paste0('All sigs\n', ct_it)
    }
  
    df_MSE_comparison <- data.frame(MSE=as.vector(mat1), QP=as.vector(mat2),
                                    sig=rep(colnames(mat1), each=nrow(mat1)))
    print(ggplot(df_MSE_comparison, aes(x=MSE, y=QP, col=sig))+
      geom_abline(intercept = 0, slope = 1, lty='dashed')+
      geom_point()+theme_bw()+ggtitle(title)+
      guides(col=F))
    
    
    mat1 <- do.call('rbind', MSE_ROO2[[ct_it]]@count_matrices_active)
    mat2 <- read_info_list[[ct_it]]$dataset_active_sigs$Y
    not_same_sigs <- (all(dim(mat1) == dim(mat2)))
    
    colnames(mat2)[!(colnames(mat2) %in% colnames(mat1))]
    colnames(mat1)[!(colnames(mat1) %in% colnames(mat2))]
    remove_na <- function(i) i[!is.na(i)]
    mat2 <- mat2[,remove_na(match(colnames(mat1), colnames(mat2)))]
    mat1 <- mat1[,remove_na(match(colnames(mat2), colnames(mat1)))]
    
    if(not_same_sigs){
      if(nrow(mat1) != nrow(mat2)){
        print(ggplot())
        next
      }
      title=paste0('Active sigs\n', ct_it, '\n Not the same number of sigs')
    }else{
      title <- paste0('Active sigs\n', ct_it)
    }
    df_MSE_comparison <- data.frame(MSE=as.vector(mat1), QP=as.vector(mat2),
                                    sig=rep(colnames(mat1), each=nrow(mat1)))
    print(ggplot(df_MSE_comparison, aes(x=MSE, y=QP, col=sig))+
      geom_abline(intercept = 0, slope = 1, lty='dashed')+
      geom_point()+theme_bw()+ggtitle(title))
  }
}
dev.off()
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##

nucleotide1 <- sapply(read_info_list, `[`, 'fullREDM_nucleotide1')
names_trinucleotide <- vector_cats_to_logR(colnames(read_info_list[[1]]$dataset_nucleotidesubstitution1$Y))
nucleotide1 <- sapply(read_info_list, `[`, 'fullREDM_nucleotide1')
names(nucleotide1) <- names(read_info_list)

system("open ../../results/results_TMB/pcawg/reports_per_cancer_type/")
pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/nucleotide1_betas.pdf"),
    height = 3, width = 3)
for(ct_it in  names(read_info_list)){
  print(plot_betas(nucleotide1[[ct_it]], names_cats = names_trinucleotide, title = ct_it))
}
dev.off()

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/nucleotide1_sigma.pdf"),
    height = 3, width = 3)
for(ct_it in  names(read_info_list)){
  print(plot_covariance_mat(nucleotide1[[ct_it]], names_cats = names_trinucleotide, title = ct_it,
                            lims=c(-1, 1)))
}
dev.off()
pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/nucleotide1_sigmanoclust.pdf"),
    height = 3, width = 3)
for(ct_it in  names(read_info_list)){
  print(plot_covariance_mat(nucleotide1[[ct_it]], names_cats = names_trinucleotide, title = ct_it,
                            arg_cluster_rows = F, arg_cluster_cols = F, lims=c(-1, 1)))
}
dev.off()

betas_nucleotides <- lapply(nucleotide1, function(i) plot_betas(i, return_df = T))
betas_nucleotides <- lapply(betas_nucleotides, function(i){
  i$LogR <- names_trinucleotide[i$LogR]
  # rownames(i) <- make.names(i$LogR, unique = T)
  i
})

plot_betas(nucleotide1$`Ovary-AdenoCA`, names_cats = names_trinucleotide)
give_barplot_from_obj(read_info_list$`Ovary-AdenoCA`$dataset_nucleotidesubstitution1, only_normalised = F, nrow_plot = 1)
give_barplot_from_obj(read_info_list$`Ovary-AdenoCA`$dataset_nucleotidesubstitution1, only_normalised = T, nrow_plot = 1)

betas_nucleotides_slopes <- do.call('cbind', lapply(betas_nucleotides, function(i) i%>% filter(type_beta == 'Slope' ) %>% select(Estimate)))
rownames(betas_nucleotides_slopes) <- names_trinucleotide
betas_nucleotides_intercepts <- do.call('cbind', lapply(betas_nucleotides, function(i) i%>% filter(type_beta == 'Intercept' ) %>% select(Estimate)))
rownames(betas_nucleotides_intercepts) <- names_trinucleotide
colnames(betas_nucleotides_intercepts) <- colnames(betas_nucleotides_slopes) <- names(read_info_list)

pheatmap::pheatmap(betas_nucleotides_slopes)
pheatmap::pheatmap(betas_nucleotides_intercepts)
# ggplot(melt(betas_nucleotides) %>% filter(type_beta == 'Slope', ), aes(x=L1, y=value))+
#   geom_point()

betas_nucleotides_slopes_cors <- outer(1:nrow(betas_nucleotides_slopes),
                                       1:nrow(betas_nucleotides_slopes), Vectorize(function(i,j){
  cor(x = unlist(betas_nucleotides_slopes[i,]), y = unlist(betas_nucleotides_slopes[j,]))
}))


library(ComplexHeatmap)
rownames(betas_nucleotides_slopes_cors) <- colnames(betas_nucleotides_slopes_cors) <- gsub(">", "$>$", rownames(betas_nucleotides_slopes))
system("open ../../results/results_TMB/pcawg/reports_per_cancer_type/")
tikzDevice::tikz("../../results/results_TMB/pcawg/reports_per_cancer_type/cors_trinucleotide.tex", height = 3, width = 2)
ht <- ComplexHeatmap::Heatmap(betas_nucleotides_slopes_cors, name = 'Correlation', cluster_rows = F,
                              heatmap_legend_param = list(direction = "horizontal"), 
                              column_names_gp = grid::gpar(fontsize = 8),
                              row_names_gp = grid::gpar(fontsize = 8))
draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom")
# pheatmap::pheatmap(betas_nucleotides_slopes_cors, cex=0.9, cluster_rows = F)
dev.off()

rownames(betas_nucleotides_slopes) <- gsub(">", "$>$", rownames(betas_nucleotides_slopes))
tikzDevice::tikz("../../results/results_TMB/pcawg/reports_per_cancer_type/cors_trinucleotide2.tex", height = 3, width = 4)
ggplot(melt(as(betas_nucleotides_slopes, 'matrix')), aes(x=Var2, col=Var1, y=value))+geom_point()+
  geom_hline(yintercept = 0, lty='dashed')+theme_bw()+geom_line(aes(group=Var1))+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank())+
  labs(y='Beta slope')+guides(col=guide_legend(nrow=2,byrow=TRUE))+
  scale_color_manual(values = nucleotide_colours_logR)
dev.off()

tikzDevice::tikz("../../results/results_TMB/pcawg/reports_per_cancer_type/cors_trinucleotide3.tex", height = 3, width = 5)
ggplot(melt(as(betas_nucleotides_slopes, 'matrix')), aes(x=Var2, col=Var1, y=value))+geom_point()+
  geom_hline(yintercept = 0, lty='dashed')+theme_bw()+geom_line(aes(group=Var1))+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank())+
  labs(y='Beta slope')+guides(col=guide_legend(nrow=1,byrow=TRUE))+
  scale_color_manual(values = nucleotide_colours_logR)
dev.off()

tikzDevice::tikz("../../results/results_TMB/pcawg/reports_per_cancer_type/cors_trinucleotide3sorted.tex", height = 3, width = 5)
ggplot(melt(as(betas_nucleotides_slopes, 'matrix')),
       aes(x=factor(Var2,levels=names(sort(colMeans(betas_nucleotides_slopes)))),
           col=Var1, y=value))+geom_point()+
  geom_hline(yintercept = 0, lty='dashed')+theme_bw()+geom_line(aes(group=Var1))+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank())+
  labs(y='Beta slope')+guides(col=guide_legend(nrow=1,byrow=TRUE))+
  scale_color_manual(values = nucleotide_colours_logR)
dev.off()

tikzDevice::tikz("../../results/results_TMB/pcawg/reports_per_cancer_type/cors_trinucleotide3sorted_v2.tex", height = 3, width = 5)
ggplot(melt(as(betas_nucleotides_slopes, 'matrix')),
       aes(x=factor(Var2,levels=names(sort(colMeans(betas_nucleotides_slopes)))),
           col=Var1, y=value))+geom_point()+
  geom_hline(yintercept = 0, lty='dashed')+theme_bw()+geom_line(aes(group=Var1))+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank())+
  # labs(y=("$\\widehat{\\betab}_1$"))+
  labs(y=("$\\hat{\\beta}_1$"))+
  guides(col=guide_legend(nrow=1,byrow=TRUE))+
  scale_color_manual(values = nucleotide_colours_logR)
dev.off()


library(GGally)
ggpairs(data.frame(t(betas_nucleotides_slopes)))

## general trends
lapply(read_info_list, function(i) normalise_cl(sapply(split_matrix_in_half(i$dataset_nucleotidesubstitution1$Y), colSums)))

pvals_nucleotide1 <- sapply(nucleotide1, wald_TMB_wrapper, fail_non_converged = F)
pvals_nucleotide1 <- sapply(nucleotide1, wald_TMB_wrapper, fail_non_converged = T)
pvals_nucleotide1[p.adjust(pvals_nucleotide1) > 0.05]
pvals_nucleotide1[is.na(pvals_nucleotide1)]

new_list_for_pert <- function(lst, nme){
  lst2 <- replicate(nme, list())
  lst2[nme] <- lst
  lst2
}
pert_nucleotide1 <- lapply(1:length(nucleotide1), function(i) give_min_pert(idx_sp = i, list_runs = nucleotide1,
                                                        logR_names_vec = new_list_for_pert(list(names_trinucleotide), i)))
# give_min_pert(idx_sp = 1, list_runs = nucleotide1, logR_names_vec = list(colnames(read_info_list[[1]]$dataset_nucleotidesubstitution1$Y)))
pert_nucleotide1 <- do.call('rbind', sapply(pert_nucleotide1, `[`, 'betas_perturbed'))
rownames(pert_nucleotide1) <- names(nucleotide1)
pert_nucleotide1
#dataset_nucleotidesubstitution1 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution1", path_to_data = "../../data/"),
##-----------------------------------------------------------------------------------------------------##

## Differential precision in data
pdf("../../results/results_TMB/pcawg/reports_per_cancer_type/ternary_plots_highestsigs.pdf", height = 3, width = 7)
for(ct in enough_samples){
  ct1 <- read_info_list[[ct]]$dataset_active_sigs$Y[read_info_list[[ct]]$dataset_active_sigs$x[,2] == 0,]
  ct2 <- read_info_list[[ct]]$dataset_active_sigs$Y[read_info_list[[ct]]$dataset_active_sigs$x[,2] == 1,]
  sigs_select <- names(sort(colSums(ct1+ct2), decreasing = T)[1:3])
  par(mfrow=c(1,2), mar=c(0,0,0,0))
  give_ternary(probs =  rm_na_rows(normalise_rw(ct1[,sigs_select])), add_par = F, opacity = 1, main=paste0('Clonal\n', ct, '\nLarge sigs'),
               col='black', pch=19, cex=.5)
  give_ternary(probs =  rm_na_rows(normalise_rw(ct2[,sigs_select])), add_par = F, opacity = 1, main=paste0('Subclonal\n', ct, '\nLarge sigs'),
               col='black', pch=19, cex=.5)
  sigs_select <- names(sort(colSums(ct1+ct2), decreasing = F)[1:3])
  give_ternary(probs =  rm_na_rows(normalise_rw(ct1[,sigs_select])), add_par = F, opacity = 1, main=paste0('Clonal\n', ct, '\nSmall sigs'),
               col='black', pch=19, cex=.5)
  give_ternary(probs =  rm_na_rows(normalise_rw(ct2[,sigs_select])), add_par = F, opacity = 1, main=paste0('Subclonal\n', ct, '\nSmall sigs'),
               col='black', pch=19, cex=.5)
}
dev.off()

##-----------------------------------------------------------------------------------------------------##
## Only using tracksig exposures
tracksig_exposures <- list()
fullDM_tracksig <- list()
diagDM_tracksig <- list()

tracksig_sigs <- list(Bone_Osteosarc = c('SBS2+SBS13', 'SBS3', 'SBS5', 'SBS40'),
                      Breast_AdenoCA = c('SBS1', 'SBS2+SBS13', 'SBS3', 'SBS5', 'SBS18'),
                      CNS_GBM = c('SBS1', 'SBS5', 'SBS40'),
                      CNS_Medullo =c('SBS1', 'SBS5', 'SBS18', 'SBS39', 'SBS40'),        
                      #CNS-PiloAstro
                      ColoRect_AdenoCA=c('SBS1', 'SBS5', 'SBS18', 'SBS40', 'SBS44'),
                      Eso_AdenoCA = c('SBS1', 'SBS2+SBS13', 'SBS3', 'SBS5', 'SBS17a+SBS17b', 'SBS18', 'SBS40'),
                      Head_SCC = c('SBS1', 'SBS2+SBS13', 'SBS5', 'SBS40'),
                      Kidney_ChRCC = c('SBS1', 'SBS2+SBS13', 'SBS5', 'SBS40'),
                      Kidney_RCC.clearcell = c('SBS1', 'SBS5', 'SBS40', 'SBS41'),
                      Kidney_RCC.papillary = c('SBS1', 'SBS5', 'SBS40', 'SBS41'),
                      Liver_HCC = c('SBS4', 'SBS5', 'SBS12', 'SBS16', 'SBS29', 'SBS40'),
                      Lung_SCC = c('SBS2+SBS13', 'SBS4', 'SBS5'),
                      Lymph_BNHL = c('SBS5', 'SBS9', 'SBS17a+SBS17b', 'SBS40'),
                      Lymph_CLL = c('SBS1', 'SBS5', 'SBS9', 'SBS40'),
                      Ovary_AdenoCA = c('SBS1', 'SBS2+SBS13', 'SBS3', 'SBS5', 'SBS40'), 
                      Panc_AdenoCA = c('SBS1', 'SBS2+SBS13', 'SBS3', 'SBS5', 'SBS17a+SBS17b', 'SBS18', 'SBS40'),
                      Panc_Endocrine = c('SBS1', 'SBS3', 'SBS5', 'SBS8', 'SBS30', 'SBS36'),
                      Prost_AdenoCA = c('SBS1', 'SBS3', 'SBS5', 'SBS18', 'SBS40'),
                      Skin_Melanoma.cutaneous=c('SBS5', 'SBS7a+SBS7b+SBS7c+SBS7d'),
                      Stomach_AdenoCA = c('SBS1', 'SBS5', 'SBS17a+SBS17b', 'SBS18', 'SBS40'),
                      Thy_AdenoCA= c('SBS1', 'SBS2+SBS13', 'SBS5', 'SBS40'),
                      Uterus_AdenoCA = c('SBS1', 'SBS2+SBS13', 'SBS5', 'SBS40', 'SBS44'))

vector_to_amalgamate_list <- function(i){
 sapply(i, function(j){
   if(grepl('[+]', j)){
     strsplit(j, '[+]')[[1]]
   }else{
     j
   }
 }) 
}

enough_samples <- enough_samples[enough_samples != 'CNS-PiloAstro']
for(ct in enough_samples){
  tracksig_exposures[[ct]] <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                         subset_signatures = unlist(sapply(tracksig_sigs[[gsub('-', '_', ct)]], function(i) strsplit(i, '[+]'))))
  tracksig_exposures[[ct]] <- give_amalgamated_exposures_TMBobj(tracksig_exposures[[ct]],
                                   vector_to_amalgamate_list(tracksig_sigs[[gsub('-', '_', ct)]]))
  colnames(tracksig_exposures[[ct]]$Y) <- tracksig_sigs[[gsub('-', '_', ct)]]
  
  # amalgamate if needed
  fullDM_tracksig[[ct]] <- wrapper_run_TMB(object = tracksig_exposures[[ct]],
                                          model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
  saveRDS(fullDM_tracksig[[ct]],
          paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_tracksig', '.RDS'))

  diagDM_tracksig[[ct]] <- wrapper_run_TMB(object = tracksig_exposures[[ct]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
  saveRDS(diagDM_tracksig[[ct]],
          paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_tracksig', '.RDS'))

}

ff <- list.files('../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/')
ff <- ff[grepl('tracksig', ff)]
ff

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/tracksig_betas.pdf"), height = 3)
for(ct in enough_samples){
  plot_betas(fullDM_tracksig[[ct]], names_cats = vector_cats_to_logR(colnames(tracksig_exposures[[ct]]$Y)),
             title = ct)
}
dev.off()

effect_size3_SP_tracksig <- sapply(enough_samples, function(ct){
  .xx <- tracksig_exposures[[ct]]
  try(give_totalperturbation_TMBobj_sigaverage(.xx, addone=T))
})

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/tracksig_exposures.pdf"), height = 3)
for(ct in enough_samples){
  give_barplot_from_obj(tracksig_exposures[[ct]], legend_on = T, legend_bottom = T, plot=F, title =ct, only_normalised = T, nrow_plot = 1,
                        levels_signatures=sigs_cosmic, arg_title='')
}
dev.off()

########################################################################
source("../3_analysis/helper/pcawg.colour.palette.R")

tracksig_pvls <- sapply(fullDM_tracksig, wald_TMB_wrapper)
tracksig_pvlsdiag <- sapply(diagDM_tracksig, wald_TMB_wrapper)
tracksig_pvls['Thy-AdenoCA'] <- tracksig_pvlsdiag['Thy-AdenoCA']
tracksig_pvls['CNS-Medullo'] <- tracksig_pvlsdiag['CNS-Medullo']
tracksig_pvls <- p.adjust(tracksig_pvls)

tracksig = read.csv("../../data/restricted/tracksig/changepoints_stats_tracksig.csv", stringsAsFactors = FALSE)
tracksig = tracksig %>% group_by(type) %>%
  dplyr::summarize(count = n(),bool_changepoints=sum(n_changepoints > 0)) %>%
  mutate(tracksig_frac= bool_changepoints/count)
tracksig = cbind.data.frame(tracksig_pvls,
                            tracksig[match(names(tracksig_pvls), tracksig$type),],
                            effect_size3_SP=effect_size3_SP_tracksig[match(names(tracksig_pvls), names(effect_size3_SP_tracksig))])
tracksig$ct = rownames(tracksig)
tracksig$minpvals = -log2(tracksig$tracksig_pvls)

pcawg_palette <- pcawg.colour.palette(gsub("\\..*", "", tracksig$ct), scheme = "tumour.subtype")
names(pcawg_palette) <- tracksig$ct

tikzDevice::tikz("../../code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/tracksig_comparison_tracksigsigs_pval.tex", height = 3, width = 3)
ggplot(tracksig, aes(x=-log2(tracksig_pvls), y=tracksig_frac, label=ct, size=count))+
  geom_point(aes(col=ct))+geom_label_repel(size=3, col='black', max.overlaps = 2)+
  labs(x='-log2 p-value', y='Fraction of TrackSig samples\nwith some changepoint')+
  scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  theme_bw()+theme(legend.position = "bottom")+guides(col='none')+labs(size='N. obs')
dev.off()

tikzDevice::tikz("../../code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/tracksig_comparison_tracksigsigs_effectsize.tex", height = 3, width = 3)
ggplot(tracksig, aes(x=effect_size3_SP, y=tracksig_frac, label=ct))+
  geom_point(aes(size=minpvals,  col=ct))+geom_label_repel(size=3, max.overlaps = 4)+
  labs(x='Effect size', y='Fraction of TrackSig samples\nwith some changepoint', col="")+theme_bw()+
  scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  # guides(col=guide_legend(ncol=4), size=FALSE)+
  guides(col='none')+labs(size='N. obs')
dev.off()

ggplot(tracksig, aes(x=effect_size3_SP, y=-log2(tracksig_pvls), label=ct, col=ct))+
  geom_point(aes(size=minpvals))+geom_label_repel(max.overlaps = 5)+
  labs(x='Effect size', y='Fraction of TrackSig samples with some changepoint', col="")+theme_bw()+
  scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  guides(col=guide_legend(ncol=4), size=FALSE)

ggplot(tracksig, aes(x=-log2(tracksig_pvls), y=tracksig_frac, label=ct, size=count))+geom_point(aes(col=ct))+
  geom_label_repel(max.overlaps = 4, size=3, col='black')+
  labs(x='-log2 p-value of\n diagRE DMDL', y='Fraction of TrackSig samples\n with some changepoint')+
  scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  # scale_x_continuous(trans = "log2")+
  theme_bw()+theme(legend.position = "bottom")+guides(col=FALSE)+labs(size='N. obs')

ggplot(tracksig, aes(x=-log2(tracksig_pvls), y=tracksig_frac, label=ct, size=count))+geom_point(aes(col=ct))+
  geom_label_repel(size=3, col='black', max.overlaps = 4)+
  labs(x='-log2 p-value of\n diagRE DMDL (nonexo)', y='Fraction of TrackSig samples\n with some changepoint')+
  scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  theme_bw()+theme(legend.position = "bottom")+guides(col=FALSE)+labs(size='N. obs')

##--------------------------------------------------------------------------------------##
### Small signatures: how small are they?

## SBS17a, SBS17b in bone osteosarcoma
ct

colnames(read_info_list[[ct]]$dataset_active_sigs$Y)

## in both cases there are some signatures
split_matrix_in_half(read_info_list[[ct]]$dataset_active_sigs$Y[,c('SBS17a', 'SBS17b')])


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/normalised_exposures_per_ct.pdf"), height = 3)
for(ct in enough_samples){
  print(do.call('grid.arrange', c(grobs=lapply(colnames(read_info_list[[ct]]$dataset_active_sigs$Y), function(sig_it){
  give_boxplot_normalised_for_sig(read_info_list[[ct]]$dataset_active_sigs, sig = sig_it)
  }), main='dsad')))
}
dev.off()

give_boxplot_normalised_for_sig(read_info_list[[ct]]$dataset_active_sigs, sig = 'SBS17b')

do.call('grid.arrange', c(grobs=lapply(colnames(read_info_list[[ct]]$dataset_nucleotidesubstitution1$Y), function(sig_it){
  give_boxplot_normalised_for_sig(read_info_list[[ct]]$dataset_nucleotidesubstitution1, sig = sig_it)
})))

##--------------------------------------------------------------------------------------##
pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/betas_all_sorted.pdf"), height = 3)
for(ct in enough_samples){
  plot_betas(TMB_obj = read_info_list[[ct]]$diagRE_DMDL_SP, sort_by_slope = T, title = ct,
             vector_cats_to_logR(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)))
}
dev.off()

all_diagREDMDL_betas <- lapply(enough_samples, function(ct){
  plot_betas(TMB_obj = read_info_list[[ct]]$diagRE_DMDL_SP,
             names_cats= vector_cats_to_logR(colnames(read_info_list[[ct]]$dataset_active_sigs$Y)),
             return_df=T, plot=F, only_slope = T, line_zero=F, add_confint = T)
})
all_diagREDMDL_betas <- lapply(all_diagREDMDL_betas, function(i) cbind(i, numerator_LogR = gsub("/.*", "", i$LogR)))

all_diagREDMDL_betas_softmax <- lapply(all_diagREDMDL_betas,
                                       function(i) cbind.data.frame(Estimate=softmax(c(i$Estimate[i$type_beta == 'Slope'], 0)),
                                                        sig=c(i$numerator_LogR[i$type_beta == 'Slope'],
                                                              strsplit(i$LogR[1], '/')[[1]][2])))

names_sigs_unique <- gtools::mixedsort(unique(do.call('rbind', all_diagREDMDL_betas_softmax)$sig))
increases_matrix <- outer(names_sigs_unique, names_sigs_unique, Vectorize(function(sig_it1, sig_it2){
  mean(as.numeric(sapply(all_diagREDMDL_betas_softmax, function(j) j[j$sig == sig_it1,'Estimate'] > j[j$sig == sig_it2,'Estimate'] )), 
      na.rm = T)
}))
num_ct_in_common_increases_matrix <- outer(names_sigs_unique, names_sigs_unique, Vectorize(function(sig_it1, sig_it2){
  sum(!is.na(as.numeric(sapply(all_diagREDMDL_betas_softmax, function(j) j[j$sig == sig_it1,'Estimate'] > j[j$sig == sig_it2,'Estimate'] ))))
}))
colnames(increases_matrix) <- rownames(increases_matrix) <- names_sigs_unique
pheatmap_increases_matrix <- pheatmap::pheatmap(increases_matrix)
pheatmap_increases_matrix
# increases_matrix[upper.tri(increases_matrix)] <- NaN
increases_matrix_melt <- melt(increases_matrix)
increases_matrix_melt$size= melt(num_ct_in_common_increases_matrix)$value
increases_matrix_melt$Var1 = factor(increases_matrix_melt$Var1, levels=pheatmap_increases_matrix$tree_row$labels[pheatmap_increases_matrix$tree_row$order])
increases_matrix_melt$Var2 = factor(increases_matrix_melt$Var2, levels=pheatmap_increases_matrix$tree_col$labels[pheatmap_increases_matrix$tree_row$order])
# increases_matrix_melt <- increases_matrix_melt[as.numeric(increases_matrix_melt$Var1) > as.numeric(increases_matrix_melt$Var2),]
increases_matrix_melt <- increases_matrix_melt[!is.na(increases_matrix_melt$value),]

increases_matrix_melt[(increases_matrix_melt$Var1 == '7c') & (increases_matrix_melt$Var2 == '38'),]
increases_matrix_melt[(increases_matrix_melt$Var2 == '7c') & (increases_matrix_melt$Var1 == '38'),]

increases_matrix_melt$Var1 <- factor(paste0('SBS', increases_matrix_melt$Var1), levels=paste0('SBS', levels(increases_matrix_melt$Var1)))
increases_matrix_melt$Var2 <- factor(paste0('SBS', increases_matrix_melt$Var2), levels=paste0('SBS', levels(increases_matrix_melt$Var2)))
ggplot(increases_matrix_melt,
       aes(x=Var1, y=Var2, col=value,size=size))+geom_point()+theme_bw()+ # shape=(value>0.5), 
  scale_color_viridis() + theme(legend.position = "bottom", legend.box="vertical")+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='',y='', size='Number of samples in common', col='Fraction of higher coefficients')
ggsave("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/comparison_order_coefficients_heatmap.pdf",
       height = 8, width = 7)
saveRDS(increases_matrix_melt, "../../data/pcawg_robjects_cache/increases_matrix_melt.RDS")

cowplot::plot_grid(cowplot::plot_grid(plot.new(), ggplot(increases_matrix_melt,
                          aes(x=Var1, y=Var2, col=value,size=size))+geom_point()+theme_bw()+ # shape=(value>0.5), 
                     scale_color_viridis() + theme(legend.position = "left", legend.box="vertical")+
                     scale_y_discrete(position = "right")+
                     theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
                     labs(x='',y='', size='Number of signatures in common', col='Fraction of higher coefficients'),
                     nrow=2,
                     rel_heights = c(0.03, 0.9)),
                   ggdendro::ggdendrogram(pheatmap_increases_matrix$tree_row, rotate=T), rel_widths = c(3,1))
ggsave("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/comparison_order_coefficients_heatmap_with_tree.pdf",
       height = 8, width = 12)

cowplot::plot_grid(cowplot::plot_grid(plot.new(), ggplot(increases_matrix_melt,
                                                         aes(x=Var1, y=Var2, col=value,size=size))+geom_point()+theme_bw()+ # shape=(value>0.5), 
                                        scale_color_viridis() + #theme(legend.position = "left", legend.box="vertical")+
                                        scale_y_discrete(position = "right")+
                                        theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1), legend.position = "bottom")+
                                        labs(x='',y='', size='Number of signatures in common', col='Fraction of higher coefficients'),
                                      nrow=2,
                                      rel_heights = c(0.03, 0.9)),
                   plot_grid(plot.new(),ggdendro::ggdendrogram(pheatmap_increases_matrix$tree_row, rotate=T), plot.new(),
                             rel_heights = c(0.005, 0.9, 0.1), ncol=1), rel_widths = c(3,1))
ggsave("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/comparison_order_coefficients_heatmap_with_tree_bottomlegend.pdf",
       height = 8, width = 12)

pairs(t(increases_matrix[c('1', '5', '40'),]))

all_diagREDMDL_betas_softmax_allsigs <- data.frame(t(sapply(all_diagREDMDL_betas_softmax, function(i) i[match(names_sigs_unique, i$sig,),'Estimate'])), ct=enough_samples)
colnames(all_diagREDMDL_betas_softmax_allsigs) <- c(paste0('SBS', names_sigs_unique), 'ct')

all_diagREDMDL_betas_softmax_1_5_40 <- t(sapply(all_diagREDMDL_betas_softmax, function(i) i[match(names_sigs_unique, i$sig,),'Estimate']))
colnames(all_diagREDMDL_betas_softmax_1_5_40) <- paste0('SBS', names_sigs_unique)
all_diagREDMDL_betas_softmax_1_5_40 <- data.frame(all_diagREDMDL_betas_softmax_1_5_40, ct=enough_samples)
pairs(all_diagREDMDL_betas_softmax_1_5_40)

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs1sbs5_and_apobec_softmax_anno.tex",
                 height=2.5, width=2.5)
ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes(x=SBS1, y=SBS5, col=ct, label=substr(ct, 1, 3), group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+ geom_point(color='black', size=2)+ geom_point()+ geom_label_repel()+
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", col='black', size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.03, vjust=2, hjust=-.4)+
  labs(x='Beta slope of SBS1', y='Beta slope of SBS5')+guides(col='none')
dev.off()

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs1sbs5_and_apobec_softmax_noanno.tex",
                 height=2.2, width=2.2)
ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes(x=SBS1, y=SBS5, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.38, vjust=1, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS1', y='Beta slope of SBS5')+guides(col='none')
dev.off()

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs13sbs8_noanno.tex",
                 height=2.2, width=2.2)
ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes(x=SBS13, y=SBS8, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.25, vjust=1, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS13', y='Beta slope of SBS8')+guides(col='none')
dev.off()

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs2sbs8_noanno.tex",
                 height=2.2, width=2.2)
ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes(x=SBS2, y=SBS8, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.25, vjust=1, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS2', y='Beta slope of SBS8')+guides(col='none')
dev.off()

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs2sbs18_noanno.tex",
                 height=2.2, width=2.2)
ggplot(all_diagREDMDL_betas_softmax_1_5_40[!is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS2) | is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS18),],
       aes(x=SBS2, y=SBS18, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.1, vjust=1, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS2', y='Beta slope of SBS18')+guides(col='none')
dev.off()

ggplot(all_diagREDMDL_betas_softmax_1_5_40[!is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS3) | is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS8),],
       aes(x=SBS3, y=SBS8, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.1, vjust=1, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS3', y='Beta slope of SBS8')+guides(col='none')

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs40sbs5_noanno.tex",
                 height=2.2, width=2.2)
ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes(x=SBS5, y=SBS40, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.45, vjust=1, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS5', y='Beta slope of SBS40')+guides(col='none')
dev.off()

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs3sbs18_noanno.tex",
                 height=2.2, width=2.2)
ggplot(all_diagREDMDL_betas_softmax_1_5_40[!is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS3) | is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS18),],
       aes(x=SBS3, y=SBS18, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.1, vjust=0.2, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS3', y='Beta slope of SBS18')+guides(col='none')
dev.off()

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs3sbs2_noanno.tex",
                 height=2.2, width=2.2)
ggplot(all_diagREDMDL_betas_softmax_1_5_40[!is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS3) | is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS2),],
       aes(x=SBS3, y=SBS2, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.1, vjust=0.2, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS3', y='Beta slope of SBS2')+guides(col='none')
dev.off()

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs2sbs13_noanno.tex",
                 height=2.2, width=2.2)
ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes(x=SBS2, y=SBS13, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.25, vjust=0.3, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS2', y='Beta slope of SBS13')+guides(col='none')
dev.off()

tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs17asbs17b_noanno.tex",
                 height=2.2, width=2.2)
ggplot(all_diagREDMDL_betas_softmax_1_5_40[!is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS17a) | is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS17b),],
       aes(x=SBS17a, y=SBS17b, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.15, vjust=0.2, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS17a', y='Beta slope of SBS17b')+guides(col='none')
dev.off()

ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes(x=SBS2, y=SBS8, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", col='black', size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.03, vjust=2, hjust=-.4)+
  # labs(x='Beta slope of SBS1', y='Beta slope of SBS5')+
  guides(col='none')

# tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs1sbs5_and_apobec_softmax_anno_SBS40.tex",
#                  height=2.5, width=2.5)
# ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes(x=SBS1, y=SBS40, col=ct, label=substr(ct, 1, 3)))+
#   geom_abline(slope = 1, intercept = 0, lty='dashed')+ geom_point()+ geom_label_repel()+
#   theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
#   geom_smooth(method = "lm")+
#   ggpubr::stat_cor(method = "pearson", label.x = 1.1, label.y = Inf, vjust=2, hjust=1)+
#   labs(x='Beta slope of SBS1', y='Beta slope of SBS40')+guides(col='none')
# dev.off()
# 
# tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_sbs1sbs5_and_apobec_softmax_anno_SBS40_2.tex",
#                  height=2.5, width=2.5)
# ggplot(all_diagREDMDL_betas_softmax_1_5_40, aes(x=SBS5, y=SBS40, col=ct, label=substr(ct, 1, 3)))+
#   geom_abline(slope = 1, intercept = 0, lty='dashed')+ geom_point()+ geom_label_repel()+
#   theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
#   geom_smooth(method = "lm")+
#   ggpubr::stat_cor(method = "pearson", label.x = 1.1, label.y = Inf, vjust=2, hjust=1)+
#   labs(x='Beta slope of SBS5', y='Beta slope of SBS40')+guides(col='none')
# dev.off()
# 
# tikzDevice::tikz("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_SBS17.tex",
#                  height=2.5, width=2.5)
# ggplot(all_diagREDMDL_betas_softmax_allsigs, aes(x=SBS17a, y=SBS17b, col=ct, label=substr(ct, 1, 3)))+
#   geom_abline(slope = 1, intercept = 0, lty='dashed')+ geom_point()+ geom_label_repel()+
#   theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
#   geom_smooth(method = "lm")+
#   ggpubr::stat_cor(method = "pearson", label.x = 1.1, label.y = Inf, vjust=2, hjust=1)+
#   labs(x='Beta slope of SBS17a', y='Beta slope of SBS17b')+guides(col='none')
# dev.off()

ggplot(all_diagREDMDL_betas_softmax_allsigs, aes(x=SBS8, y=SBS13, col=ct, label=substr(ct, 1, 3)))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+ geom_point()+ geom_label_repel()+
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm")+
  ggpubr::stat_cor(method = "pearson", label.x = 1.1, label.y = Inf, vjust=2, hjust=1)+
  # labs(x='Beta slope of SBS17a', y='Beta slope of SBS17b')+
  guides(col='none')

ggplot(all_diagREDMDL_betas_softmax_allsigs, aes(x=SBS8, y=SBS26, col=ct, label=substr(ct, 1, 3)))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+ geom_point()+ geom_label_repel()+
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm")+
  ggpubr::stat_cor(method = "pearson", label.x = 1.1, label.y = Inf, vjust=2, hjust=1)+
  # labs(x='Beta slope of SBS17a', y='Beta slope of SBS17b')+
  guides(col='none')

ggplot(all_diagREDMDL_betas_softmax_allsigs, aes(x=SBS8, y=SBS18, col=ct, label=substr(ct, 1, 3)))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+ geom_point()+ geom_label_repel()+
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm")+
  ggpubr::stat_cor(method = "pearson", label.x = 1.1, label.y = Inf, vjust=2, hjust=1)+
  # labs(x='Beta slope of SBS17a', y='Beta slope of SBS17b')+
  guides(col='none')

ggplot(all_diagREDMDL_betas_softmax_allsigs, aes(x=SBS3, y=SBS13, col=ct, label=substr(ct, 1, 3)))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+ geom_point()+ geom_label_repel()+
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm")+
  ggpubr::stat_cor(method = "pearson", label.x = 1.1, label.y = Inf, vjust=2, hjust=1)+
  # labs(x='Beta slope of SBS17a', y='Beta slope of SBS17b')+
  guides(col='none')

all_diagREDMDL_betas_softmax_allsigs_v2 <- all_diagREDMDL_betas_softmax_allsigs
rownames(all_diagREDMDL_betas_softmax_allsigs_v2) <- all_diagREDMDL_betas_softmax_allsigs_v2$ct
all_diagREDMDL_betas_softmax_allsigs_v2 <- all_diagREDMDL_betas_softmax_allsigs_v2[,(colnames(all_diagREDMDL_betas_softmax_allsigs_v2) != 'ct')]


###------ correlation of signatures ------###
all_diagREDMDL_betas_softmax_allsigs_v2_cors <- outer(1:ncol(all_diagREDMDL_betas_softmax_allsigs_v2), 1:ncol(all_diagREDMDL_betas_softmax_allsigs_v2),
                                                      Vectorize(function(i,j){
  if( sum(!is.na(all_diagREDMDL_betas_softmax_allsigs_v2[,i]) & !is.na(all_diagREDMDL_betas_softmax_allsigs_v2[,j])) <=2){
    NA ## if there are 2 or fewer points in common. if there are 2 the correlation is possible but always 1
  }else{
    try(cor(x = unlist(all_diagREDMDL_betas_softmax_allsigs_v2[,i]), y = unlist(all_diagREDMDL_betas_softmax_allsigs_v2[,j]), use = "pairwise.complete.obs"))
  }
}))
colnames(all_diagREDMDL_betas_softmax_allsigs_v2_cors) <- rownames(all_diagREDMDL_betas_softmax_allsigs_v2_cors) <- colnames(all_diagREDMDL_betas_softmax_allsigs_v2)
rm_cols_cor_sigs <- !rowSums(is.na(all_diagREDMDL_betas_softmax_allsigs_v2_cors)) == ncol(all_diagREDMDL_betas_softmax_allsigs_v2_cors)
rm_rows_cor_sigs <- !colSums(is.na(all_diagREDMDL_betas_softmax_allsigs_v2_cors)) == nrow(all_diagREDMDL_betas_softmax_allsigs_v2_cors)
all_diagREDMDL_betas_softmax_allsigs_v2_cors <- all_diagREDMDL_betas_softmax_allsigs_v2_cors[rm_cols_cor_sigs,]
all_diagREDMDL_betas_softmax_allsigs_v2_cors <- all_diagREDMDL_betas_softmax_allsigs_v2_cors[,rm_rows_cor_sigs]

all_diagREDMDL_betas_softmax_allsigs_v2_cors_hclust <- hclust(dist(all_diagREDMDL_betas_softmax_allsigs_v2_cors))
pheatmap::pheatmap(all_diagREDMDL_betas_softmax_allsigs_v2_cors)

num_common_sigs <- outer(1:ncol(all_diagREDMDL_betas_softmax_allsigs_v2), 1:ncol(all_diagREDMDL_betas_softmax_allsigs_v2), Vectorize(function(i,j){
  try(sum(!is.na(unlist(all_diagREDMDL_betas_softmax_allsigs_v2[,i])+unlist(all_diagREDMDL_betas_softmax_allsigs_v2[,j]))))}))
colnames(num_common_sigs) <- rownames(num_common_sigs) <- colnames(all_diagREDMDL_betas_softmax_allsigs_v2)
num_common_sigs <- num_common_sigs[rm_cols_cor_sigs,]
num_common_sigs <- num_common_sigs[,rm_rows_cor_sigs]

cors_softmax_melt_all_sigs <- cbind.data.frame(cors_softmax=melt(all_diagREDMDL_betas_softmax_allsigs_v2_cors),
                                      num_common_sigs=melt(num_common_sigs))

cors_softmax_meltv2_all_sigs <- cors_softmax_melt_all_sigs[!is.na(cors_softmax_melt_all_sigs$cors_softmax.value),]
ggplot(cors_softmax_meltv2_all_sigs, aes(x=factor(num_common_sigs.Var1, levels = all_diagREDMDL_betas_softmax_allsigs_v2_cors_hclust$labels[all_diagREDMDL_betas_softmax_allsigs_v2_cors_hclust$order]),
                                         y=factor(num_common_sigs.Var2, levels=all_diagREDMDL_betas_softmax_allsigs_v2_cors_hclust$labels[all_diagREDMDL_betas_softmax_allsigs_v2_cors_hclust$order]),
                                         col=cors_softmax.value, size=num_common_sigs.value))+
  geom_point()+scale_color_viridis() + theme_bw()+theme(legend.position = "bottom", legend.box="vertical")+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='',y='', size='Number of signatures in common', col='Pearson correlation')
ggsave("../../code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/comparison_order_coefficients_heatmap_all_sigs.pdf", height=5.5, width = 5.2)
saveRDS(cors_softmax_meltv2_all_sigs, "../../data/pcawg_robjects_cache/cors_softmax_meltv2_all_sigs.RDS")

tikzDevice::tikz("../../code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/comparison_order_coefficients_heatmap_all_sigs.tex", height=6.5, width = 6)
ggplot(cors_softmax_meltv2_all_sigs, aes(x=factor(num_common_sigs.Var1, levels = all_diagREDMDL_betas_softmax_allsigs_v2_cors_hclust$labels[all_diagREDMDL_betas_softmax_allsigs_v2_cors_hclust$order]),
                                         y=factor(num_common_sigs.Var2, levels=all_diagREDMDL_betas_softmax_allsigs_v2_cors_hclust$labels[all_diagREDMDL_betas_softmax_allsigs_v2_cors_hclust$order]),
                                         col=cors_softmax.value, size=num_common_sigs.value))+
  geom_point()+scale_color_viridis() + theme_bw()+theme(legend.position = "bottom", legend.box="vertical")+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='',y='', size='Number of signatures in common', col='Pearson correlation')
dev.off()

###------ correlation of samples ------###
all_diagREDMDL_betas_softmax_allsigs_v2_corssamps <- outer(1:nrow(all_diagREDMDL_betas_softmax_allsigs_v2), 1:nrow(all_diagREDMDL_betas_softmax_allsigs_v2),
                                                      Vectorize(function(i,j){
                                                        if( sum(!is.na(all_diagREDMDL_betas_softmax_allsigs_v2[i,]) & !is.na(all_diagREDMDL_betas_softmax_allsigs_v2[j,])) <=2){
                                                          NA ## if there are 2 or fewer points in common. if there are 2 the correlation is possible but always 1
                                                        }else{
                                                          try(cor(x = unlist(all_diagREDMDL_betas_softmax_allsigs_v2[i,]), y = unlist(all_diagREDMDL_betas_softmax_allsigs_v2[j,]), use = "pairwise.complete.obs"))
                                                        }
                                                      }))
colnames(all_diagREDMDL_betas_softmax_allsigs_v2_corssamps) <- rownames(all_diagREDMDL_betas_softmax_allsigs_v2_corssamps) <- rownames(all_diagREDMDL_betas_softmax_allsigs_v2)
all_diagREDMDL_betas_softmax_allsigs_v2_corssamps_hclust <- hclust(dist(all_diagREDMDL_betas_softmax_allsigs_v2_corssamps))

num_common_sigssamps <- outer(1:nrow(all_diagREDMDL_betas_softmax_allsigs_v2), 1:nrow(all_diagREDMDL_betas_softmax_allsigs_v2), Vectorize(function(i,j){
  try(sum(!is.na(unlist(all_diagREDMDL_betas_softmax_allsigs_v2[i,])+unlist(all_diagREDMDL_betas_softmax_allsigs_v2[j,]))))}))
colnames(num_common_sigssamps) <- rownames(num_common_sigssamps) <- colnames(all_diagREDMDL_betas_softmax_allsigs_v2_corssamps)
cors_softmax_melt_all_sigs_samps <- cbind.data.frame(cors_softmax=melt(all_diagREDMDL_betas_softmax_allsigs_v2_corssamps),
                                               num_common_sigs=melt(num_common_sigssamps))
cors_softmax_melt_all_sigs_samps <- cors_softmax_melt_all_sigs_samps[!is.na(cors_softmax_melt_all_sigs_samps$cors_softmax.value),]

ggplot(cors_softmax_melt_all_sigs_samps, aes(x=factor(num_common_sigs.Var1, levels = all_diagREDMDL_betas_softmax_allsigs_v2_corssamps_hclust$labels[all_diagREDMDL_betas_softmax_allsigs_v2_corssamps_hclust$order]),
                                         y=factor(num_common_sigs.Var2, levels=all_diagREDMDL_betas_softmax_allsigs_v2_corssamps_hclust$labels[all_diagREDMDL_betas_softmax_allsigs_v2_corssamps_hclust$order]),
                                         col=cors_softmax.value, size=num_common_sigs.value))+
  geom_point()+scale_color_viridis() + theme_bw()+theme(legend.position = "bottom", legend.box="vertical")+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='',y='', size='Number of signatures in common', col='Pearson correlation')
ggsave("../../code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/correlations_from_beta_slopes_softmax2_all_sigs.pdf", height=6, width = 5.7)


ggplot(all_diagREDMDL_betas_softmax_1_5_40[!is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS3) | is.na(all_diagREDMDL_betas_softmax_1_5_40$SBS2),],
       aes(x=SBS3, y=SBS2, col=ct, group=1))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
  theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
  geom_smooth(method = "lm", size=0.2)+
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.1, vjust=0.2, hjust=-.02, size=3.5, label.sep='\n')+
  labs(x='Beta slope of SBS3', y='Beta slope of SBS2')+guides(col='none')

rownames(all_diagREDMDL_betas_softmax_1_5_40) <- all_diagREDMDL_betas_softmax_1_5_40$ct

pdf("../../code/2_inference_TMB/summary_TMB_PCAWG_SP_files/figure-latex/scatterplot_from_beta_slopes_softmax2_all_sigs.pdf", height=3, width = 3)
for(i in 1:ncol(all_diagREDMDL_betas_softmax_1_5_40)){
  for(j in i:ncol(all_diagREDMDL_betas_softmax_1_5_40)){
    isig <- colnames(all_diagREDMDL_betas_softmax_1_5_40)[i]
    jsig <- colnames(all_diagREDMDL_betas_softmax_1_5_40)[j]
    print(ggplot(all_diagREDMDL_betas_softmax_1_5_40[!is.na(all_diagREDMDL_betas_softmax_1_5_40[,isig]) | is.na(all_diagREDMDL_betas_softmax_1_5_40),],
           aes_string(x=isig, y=jsig, col='ct', group="1"))+
      geom_abline(slope = 1, intercept = 0, lty='dashed')+  geom_point(color='black', size=2)+ geom_point()+ 
      theme_bw()+ scale_color_manual(values = pcawg_palette)+theme(legend.position = "bottom")+
      geom_smooth(method = "lm", size=0.2)+
      ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 0.1, vjust=0.2, hjust=-.02, size=3.5, label.sep='\n')+
      labs(x=paste0('Beta slope of ', isig), y=paste0('Beta slope of ', jsig))+guides(col='none'))
  }
}
dev.off()

##-----------------------------------------------------------------------------------------------------##
library(YAPSA)
## Extract signatures with YAPSA (which contain confidence intervals) and check if the most extreme
## betas correspond to signatures with high confidence intervals

ct = enough_samples[1]

YAPSA_res <- YAPSA::LCD(in_mutation_catalogue_df = t(read_info_list[[ct]]$dataset_nucleotidesubstitution3$Y),
           in_signatures_df = sigs_cosmic0)
idx_sig=1
idx_sample = 1
CIs <- outer(1:ncol(YAPSA_res), 1:nrow(YAPSA_res), function(idx_sample, idx_sig){
  confIntExp(in_ind = idx_sig, in_exposure_vector = YAPSA_res[,idx_sample],
           in_catalogue_vector = t(read_info_list[[ct]]$dataset_nucleotidesubstitution3$Y)[,idx_sample],
           in_sig_df=sigs_cosmic0,  in_sigLevel = 0.05, in_delta = 0.4)
})
YAPSA_res$`f221c897-6ad0-0df9-e040-11ac0c4813ef`
