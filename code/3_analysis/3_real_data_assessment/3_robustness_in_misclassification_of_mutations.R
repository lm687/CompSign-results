##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
source("helper.R")
source("../../1_create_ROO/helper_1_create_ROO.R") ## for QPsig

library(TMB)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(dplyr)
library(jcolors)
library(viridis)
library(mutSigExtractor)
library(reshape2)
library(CompSign)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
##' Including all signatures, even with artefacts (SBS_SIGNATURE_PROFILES_V2 and SBS_SIGNATURE_PROFILES_V3 from
##' mutSigExtractor do not include artefactual signatures)
sigs_cosmic0 <- read.table(paste0( "../../../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
                           stringsAsFactors = FALSE, sep = ',', header = TRUE)
rownames(sigs_cosmic0) <- paste0(substr(sigs_cosmic0$SubType, 1, 1),'[',
                                 sigs_cosmic0$Type, ']', substr(sigs_cosmic0$SubType, 3, 3))
sigs_cosmic0 <- sigs_cosmic0[-c(1,2)];
sigs_cosmic <- colnames(sigs_cosmic0)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../../data/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
source("../../3_analysis/helper/pcawg.colour.palette.R") ## pcawg_palette
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
read_info <- function(ct){
  .x <- list(diagRE_DMDL_SP = try(readRDS(paste0("../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS"))),
             diagRE_DMDL_nucleotides = try(readRDS(paste0("../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_nucleotidesubstitution1.RDS"))),
             # fullRE_DMDL_SP = try(readRDS(paste0("../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS"))),
             dataset_active_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../data/", override_warning_X_Z = T),
             dataset_nucleotidesubstitution3 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution3", path_to_data = "../../../data/", override_warning_X_Z = T),
             dataset_nucleotidesubstitution1 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution1", path_to_data = "../../../data/", override_warning_X_Z = T))
  .x
}
read_info_list <- lapply(enough_samples, function(ct){
  read_info(ct)
}); names(read_info_list) <- enough_samples
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## From the PCAWG dataset, as shown in the paper
resDM_sigs <- sapply(read_info_list, function(i) i$diagRE_DMDL_SP, simplify = F)
resDM_sigs_tests <-  sapply(resDM_sigs, CompSign::wald_TMB_wrapper, simplify = F)
resDM_nucleotides <- sapply(read_info_list, function(i) i$diagRE_DMDL_SP, simplify = F)
resDM_nucleotides_tests <-  sapply(resDM_nucleotides, CompSign::wald_TMB_wrapper, simplify = F)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
pcawg_palette <- pcawg.colour.palette(x = gsub("\\..*", "", names(read_info_list)),  scheme = "tumour.subtype")
pcawg_palette_2 <- c(  "Bone-Osteosarc"="#FFD700", "Breast-AdenoCA"="#CD6090" , "CNS-GBM"="#3D3D3D",
                       "CNS-Medullo"="#D8BFD8", "CNS-PiloAstro"="#B0B0B0", "ColoRect-AdenoCA"="#191970",
                       "Eso-AdenoCA"="#1E90FF",  "Head-SCC" = "#8B2323", "Kidney-ChRCC"= "#B32F0B", 
                       "Kidney-RCC.clearcell"= "#FF4500", "Kidney-RCC.papillary"= "#FF4500",    "Liver-HCC"= "#006400",
                       "Lung-SCC" =  "#FDF5E6", "Lymph-BNHL"= "#698B22", "Lymph-CLL"= "#F4A35D", "Ovary-AdenoCA"= "#008B8B",
                       "Panc-AdenoCA" = "#7A378B", "Panc-Endocrine"= "#E066FF", "Prost-AdenoCA"= "#87CEFA",
                       "Skin-Melanoma.cutaneous"= "#000000",  "Stomach-AdenoCA"= "#BFEFFF" , "Thy-AdenoCA"= "#9370DB",
                       "Uterus-AdenoCA"=   "#FF8C69")
cbind(pcawg_palette, pcawg_palette_2)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
nreplicates <- 5
percentages_misclassification <- c(0.05, 0.10, 0.20, 0.4)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
give_missclassification_results <- function(percentages_misclassification_ct, ct=ct, signature_fitting_method){
  original_matrix = read_info_list[[ct]]$dataset_nucleotidesubstitution3
  ## this simulation only works if the data are organised in a specific way
  stopifnot(all(original_matrix$x[,2] == rep(c(0,1), each=nrow(original_matrix$x)/2)))
  stopifnot(all(original_matrix$z == rbind(diag(1, nrow = nrow(original_matrix$x)/2), diag(1, nrow = nrow(original_matrix$x)/2))))
  ## for each patient, select a percentage of counts and misclassify them, i.e. put them in the subclonal row if they are clonal, and in the clonal group if they are subclonal
  counts_to_misclassify = sample(1:sum(original_matrix$Y), size = round(sum(original_matrix$Y)*percentages_misclassification_ct),
                                 replace = F)
  new_dataset_obj_trinucleotide <- give_missclassified_matrix(original_matrix, counts_to_misclassify=counts_to_misclassify, missclassify = T)
  new_nucleotide_matrix = new_dataset_obj_trinucleotide
  ## get nucleotide changes from trinucleotides (i.e. aggregating trinucleotides)
  new_nucleotide_matrix$Y = sapply(sort(unique(sapply(colnames(new_nucleotide_matrix$Y), function(i) substr(i, 3, 5)))), function(SBS1_changes){
    rowSums(new_nucleotide_matrix$Y[,grepl(SBS1_changes, colnames(new_nucleotide_matrix$Y))])
  })
  
  ## total number of mutations per patient
  stopifnot(rowSums(matrix(rowSums(original_matrix$Y), ncol=2)) == rowSums(matrix(rowSums(new_nucleotide_matrix$Y), ncol=2)))
  
  list(signatures_matrix=extract_sigs_TMB_obj(dataset_obj_trinucleotide=new_dataset_obj_trinucleotide,
                                              subset_signatures = colnames(read_info_list[[ct]]$dataset_active_sigs$Y),
                                              signature_version=NULL, signature_definition = sigs_cosmic0,
                                              signature_fitting_method = signature_fitting_method),
       nucleotide_matrix=new_nucleotide_matrix)
}


give_missclassified_exposures_and_nucleotides <- function(signature_fitting_method){
  sapply(enough_samples, function(ct){
    ## Creating exposures with a percentage of misclassified counts
    x <- lapply(percentages_misclassification, give_missclassification_results, ct=ct, signature_fitting_method=signature_fitting_method)
    names(x) <- paste0('misclassified', percentages_misclassification)
    x
  }, simplify = F)
}
##-----------------------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------------------##
## Run the model

give_digRE_results <- function(misclassification_signature_exposures_list, model_arg="diagRE_DM"){
  sapply(enough_samples, function(ct){
    sapply(misclassification_signature_exposures_list[[ct]], function(x){
      try({CompSign::wrapper_run_TMB(object = x,
                                model = model_arg, use_nlminb=T, smart_init_vals=F)})
    }, simplify = F)
  }, simplify = F)
}

# repl=1
# add_to_name = 'QP'
# ct = 'Bone-Osteosarc'
# misclassification_signature_exposures <- list()
# misclassification_signature_exposures_signatures <- list()
# misclassification_signature_exposures_nucleotide <- list()
# diagDM_misclassification_signature_exposures_nucleotide <- list()
# diagDM_misclassification_signature_exposures_signatures <- list()
# for(repl in 1:nreplicates){
#   misclassification_signature_exposures[[repl]] <- readRDS(file = paste0("../../../data/assessing_models_real_data/simulated_datasets/3_robustness_in_misclassification_of_mutations/3_misclassification_signature_exposures_repl", repl, add_to_name, ".RDS"))
#   misclassification_signature_exposures_signatures[[repl]] = try(sapply(misclassification_signature_exposures[[repl]], function(i) sapply(i, `[[`, 'signatures_matrix', simplify = F), simplify = F))
#   misclassification_signature_exposures_nucleotide[[repl]] = try(sapply(misclassification_signature_exposures[[repl]], function(i) sapply(i, `[[`, 'nucleotide_matrix', simplify = F), simplify = F))
# 
#   diagDM_misclassification_signature_exposures_nucleotide[[repl]] <- readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/3_robustness_in_misclassification_of_mutations/3_diagDM_misclassification_signature_exposures_nucleotide_repl", repl, add_to_name, ".RDS"))
#   diagDM_misclassification_signature_exposures_signatures[[repl]] <- readRDS(file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/3_robustness_in_misclassification_of_mutations/3_diagDM_misclassification_signature_exposures_signatures_repl", repl, add_to_name, ".RDS"))
# }
# 
# diagDM_misclassification_signature_exposures_nucleotide[[repl]][[ct]]$misclassified0.4

generate_datasets <- function(signature_fitting_method){
  if(signature_fitting_method == 'mutSigExtractor'){
    add_to_name = ''
  }else{
    add_to_name = signature_fitting_method
  }
  misclassification_signature_exposures <- list()
  misclassification_signature_exposures_signatures <- list()
  misclassification_signature_exposures_nucleotide <- list()

  for(repl in 1:nreplicates){
    cat('Replicate: ', repl, '\n')
    
    ## extract signatures or new trinucleotide counts. if there is an error, continue until we get a good dataset (for ntrials trials)
    ntrials <- 5
    ntrials_ct <- 0
    misclassification_signature_exposures[[repl]] <- ''
    while(ntrials_ct==0 | class(misclassification_signature_exposures[[repl]]) == "try-error"){
      misclassification_signature_exposures[[repl]] <- try(give_missclassified_exposures_and_nucleotides(signature_fitting_method = signature_fitting_method))
      ntrials_ct <- ntrials_ct+1
    }
    ## save dataset
    saveRDS(misclassification_signature_exposures[[repl]],
            file = paste0("../../../data/assessing_models_real_data/simulated_datasets/3_robustness_in_misclassification_of_mutations/3_misclassification_signature_exposures_repl", repl, add_to_name, ".RDS"))
    
  }
}

wrapper_missclassified_exposures_and_nucleotides <- function(signature_fitting_method, arg_model='diagRE_DM', run_only_nucleotides=F){
  if(!run_only_nucleotides){
    if(signature_fitting_method == 'mutSigExtractor'){
      add_to_name = ''
    }else{
      add_to_name = signature_fitting_method
    }
  }else{
    add_to_name = ''
  }
  
  
  if(arg_model == 'diagRE_DM'){
    add_to_name_model = 'diagDM'
  }else if(arg_model == 'fullRE_DM'){
    add_to_name_model = 'fullDM'
  }else{
    stop('Incorrect <add_to_name_model>\n')
  }

  
  misclassification_signature_exposures <- list()
  diagDM_misclassification_signature_exposures_signatures <- list()
  diagDM_misclassification_signature_exposures_nucleotide <- list()
  misclassification_signature_exposures_signatures <- list()
  misclassification_signature_exposures_nucleotide <- list()
  
  for(repl in 1:nreplicates){
    cat('Replicate: ', repl, '\n')
    
    outnucleotide=paste0("../../../data/assessing_models_real_data/inference_results/TMB/3_robustness_in_misclassification_of_mutations/3_", add_to_name_model, "_misclassification_signature_exposures_nucleotide_repl", repl, add_to_name, ".RDS")
    outsignatures=paste0("../../../data/assessing_models_real_data/inference_results/TMB/3_robustness_in_misclassification_of_mutations/3_", add_to_name_model, "_misclassification_signature_exposures_signatures_repl", repl, add_to_name, ".RDS")
    
    if(file.exists(outnucleotide) & file.exists(outsignatures)){
      cat('File exists. Skipping.\n')
      next ## nothing to do
    }
    
    ## read dataset generated in generate_datasets()
    misclassification_signature_exposures[[repl]] <- readRDS(file = paste0("../../../data/assessing_models_real_data/simulated_datasets/3_robustness_in_misclassification_of_mutations/3_misclassification_signature_exposures_repl", repl, add_to_name, ".RDS"))
    
    ## get signatures and nucleotides separately
    misclassification_signature_exposures_signatures[[repl]] = try(sapply(misclassification_signature_exposures[[repl]], function(i) sapply(i, `[[`, 'signatures_matrix', simplify = F), simplify = F))
    misclassification_signature_exposures_nucleotide[[repl]] = try(sapply(misclassification_signature_exposures[[repl]], function(i) sapply(i, `[[`, 'nucleotide_matrix', simplify = F), simplify = F))
    
    ## fit the model
    
    if(!file.exists(outnucleotide)){
      diagDM_misclassification_signature_exposures_nucleotide[[repl]] <- (give_digRE_results(misclassification_signature_exposures_nucleotide[[repl]]))
      ## save inference results
      saveRDS(diagDM_misclassification_signature_exposures_nucleotide[[repl]],
              file = outnucleotide)
    }
    if(!file.exists(outsignatures) | run_only_nucleotides){
      diagDM_misclassification_signature_exposures_signatures[[repl]] <- (give_digRE_results(misclassification_signature_exposures_signatures[[repl]]))
      ## save inference results
      saveRDS(diagDM_misclassification_signature_exposures_signatures[[repl]],
              file = outsignatures)
    }
      
    ## Check that indeed, counts have only been re-arranged within patients
    stopifnot(rowSums(matrix(rowSums(misclassification_signature_exposures_nucleotide[[repl]]$`Bone-Osteosarc`$misclassified0.05$Y), ncol=2)) ==
                rowSums(matrix(rowSums(misclassification_signature_exposures_nucleotide[[repl]]$`Bone-Osteosarc`$misclassified0.1$Y), ncol=2)))
    
  }
}
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Run the functions to generate datasets
re_run <- F
if(re_run){
  ## generate the results
  generate_datasets(signature_fitting_method = 'QP')
  generate_datasets(signature_fitting_method = 'mutSigExtractor')

  ## extract signatures and run model
  wrapper_missclassified_exposures_and_nucleotides(signature_fitting_method = 'QP')
  wrapper_missclassified_exposures_and_nucleotides(signature_fitting_method = 'mutSigExtractor')
  
  wrapper_missclassified_exposures_and_nucleotides(signature_fitting_method = NA,
                                                   run_only_nucleotides = T, arg_model = 'fullRE_DM')
  
}
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Read the results in
diagDM_misclassification_signature_exposures_signatures_QP <- lapply(1:nreplicates, function(repl)
     readRDS(paste0("../../../data/assessing_models_real_data/inference_results/TMB/3_robustness_in_misclassification_of_mutations/3_diagDM_misclassification_signature_exposures_signatures_repl",
     repl, "QP.RDS")))
diagDM_misclassification_signature_exposures_signatures_mutSigExtractor <- lapply(1:nreplicates, function(repl)
  readRDS(paste0("../../../data/assessing_models_real_data/inference_results/TMB/3_robustness_in_misclassification_of_mutations/3_diagDM_misclassification_signature_exposures_signatures_repl",
                 repl, ".RDS")))

diagDM_misclassification_signature_exposures_nucleotides_QP <- lapply(1:nreplicates, function(repl)
  readRDS(paste0("../../../data/assessing_models_real_data/inference_results/TMB/3_robustness_in_misclassification_of_mutations/3_diagDM_misclassification_signature_exposures_nucleotide_repl",
                 repl, "QP.RDS")))
diagDM_misclassification_signature_exposures_nucleotides_mutSigExtractor <- lapply(1:nreplicates, function(repl)
  readRDS(paste0("../../../data/assessing_models_real_data/inference_results/TMB/3_robustness_in_misclassification_of_mutations/3_diagDM_misclassification_signature_exposures_nucleotide_repl",
                 repl, ".RDS")))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Several replicates for the below
give_test_results <-function(k) lapply(k, function(j) sapply(j, function(ct){
  sapply(ct, function(j) try(CompSign::wald_TMB_wrapper(j)))
}, simplify = FALSE))

diagDM_misclassification_signature_exposures_signatures_tests_QP <- give_test_results(diagDM_misclassification_signature_exposures_signatures_QP)
diagDM_misclassification_signature_exposures_signatures_tests_mutSigExtractor <- give_test_results(diagDM_misclassification_signature_exposures_signatures_mutSigExtractor)
diagDM_misclassification_signature_nucleotide_tests_QP <- give_test_results(diagDM_misclassification_signature_exposures_nucleotides_QP)
diagDM_misclassification_signature_nucleotide_tests_mutSigExtractor <- give_test_results(diagDM_misclassification_signature_exposures_nucleotides_mutSigExtractor)

##-----------------------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------##
## Plots
# plt3 <- give_barplot_agreement_in_DA(resDM_sigs_tests, diagDM_misclassification_signature_nucleotide_tests_QP[[4]], ylabel='Number of datasets with misassigned mutations')
# plt4 <- give_barplot_agreement_in_DA(resDM_sigs_tests, diagDM_misclassification_signature_exposures_signatures_tests_QP[[4]], ylabel='Number of datasets with misassigned mutations')
# 
# cowplot::plot_grid(plt3, plt4)
# 
# sapply(diagDM_misclassification_signature_exposures_signatures_tests_QP, function(j){
#   give_barplot_agreement_in_DA(resDM_sigs_tests, j, ylabel='Number of datasets with misassigned mutations')
# }, simplify = F)

# do.call('grid.arrange', c(nrow=1, lapply(diagDM_misclassification_signature_exposures_signatures_tests_QP, function(j) give_barplot_agreement_in_DA(j, resDM_sigs_tests, ylabel='Number of datasets with misassigned mutations'))))
# do.call('grid.arrange', c(nrow=1, lapply(diagDM_misclassification_signature_nucleotide_tests_QP, function(j) give_barplot_agreement_in_DA(j, resDM_nucleotides_tests, ylabel='Number of datasets with misassigned mutations'))))

give_agreement_in_DA_percentages_melt <- function(diagRE_DM_tests, df, return_df=F){
  diagDM_leave_one_out_exposures_tests_accordance_DA = melt(sapply(enough_samples, function(ct) t(t(rep(diagRE_DM_tests[[ct]] <= 0.05, length(df[[ct]])) + (df[[ct]] <= 0.05))), simplify = F, USE.NAMES = T))
  diagDM_leave_one_out_exposures_tests_accordance_DA$accordance = ifelse(diagDM_leave_one_out_exposures_tests_accordance_DA$value == 1, 'Discordance', 'Accordance')
  return(diagDM_leave_one_out_exposures_tests_accordance_DA)
}

give_agreement_in_DA_percentages <- function(diagRE_DM_tests, df, return_df=F){
  diagDM_leave_one_out_exposures_tests_accordance_DA <- give_agreement_in_DA_percentages_melt(diagRE_DM_tests, df)
  diagDM_leave_one_out_exposures_tests_accordance_DA$percentage_misclassified <- as.numeric(gsub("misclassified", "", diagDM_leave_one_out_exposures_tests_accordance_DA$Var1))
  diagDM_leave_one_out_exposures_tests_accordance_DA <- diagDM_leave_one_out_exposures_tests_accordance_DA %>%
    group_by(percentage_misclassified) %>% summarise(percent_accordance=mean(accordance == 'Accordance', na.rm=T))
  head(diagDM_leave_one_out_exposures_tests_accordance_DA)
  if(return_df){
    diagDM_leave_one_out_exposures_tests_accordance_DA
  }else{
    ggplot(diagDM_leave_one_out_exposures_tests_accordance_DA,
           aes(x=percentage_misclassified, y = percent_accordance))+
      geom_line()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
      # scale_fill_manual(values = c(pcawg_palette, "white"="white"), na.value="cyan")+
      # scale_color_manual(values = c('1'='black'))+guides(fill='none', color='none')+
      labs(y='Percentage of accordance across cancer types', x='Percentage of misclassified mutations')
  }
}

give_agreement_in_DA_percentages_perct <- function(diagRE_DM_tests, df, return_df=F){
  diagDM_leave_one_out_exposures_tests_accordance_DA <- give_agreement_in_DA_percentages_melt(diagRE_DM_tests, df)
  diagDM_leave_one_out_exposures_tests_accordance_DA$percentage_misclassified <- diagDM_leave_one_out_exposures_tests_accordance_DA$Var1
  diagDM_leave_one_out_exposures_tests_accordance_DA <- diagDM_leave_one_out_exposures_tests_accordance_DA %>%
    group_by(percentage_misclassified, L1) %>% summarise(percent_accordance=mean(accordance == 'Accordance', na.rm=T)) %>%
    ungroup()
  diagDM_leave_one_out_exposures_tests_accordance_DA$ct = diagDM_leave_one_out_exposures_tests_accordance_DA$L1
  diagDM_leave_one_out_exposures_tests_accordance_DA$L1 = NULL
  if(return_df){
    diagDM_leave_one_out_exposures_tests_accordance_DA
  }else{
    stop()
  }
}


# cowplot::plot_grid(give_agreement_in_DA_percentages(resDM_nucleotides_tests, diagDM_misclassification_signature_exposures_nucleotide_tests),
#                    give_agreement_in_DA_percentages(resDM_sigs_tests, diagDM_misclassification_signature_exposures_signatures_tests))


df_misclassification_percentage_QP <- rbind.data.frame(cbind.data.frame(melt(sapply(diagDM_misclassification_signature_nucleotide_tests_QP, function(replicate) give_agreement_in_DA_percentages(replicate, resDM_nucleotides_tests, return_df=T), simplify = F), id.vars=c('percentage_misclassified', 'percent_accordance')),
                                                                     category='nucleotides'),
                                                    cbind.data.frame(melt(sapply(diagDM_misclassification_signature_exposures_signatures_tests_QP, function(replicate) give_agreement_in_DA_percentages(replicate, resDM_sigs_tests, return_df=T), simplify = F), id.vars=c('percentage_misclassified', 'percent_accordance')),
                                                                     category='signatures'))
df_misclassification_percentage_QP$percentage_misclassified = percentages_misclassification[df_misclassification_percentage_QP$percentage_misclassified]
colnames(df_misclassification_percentage_QP)[colnames(df_misclassification_percentage_QP) == 'L1'] = 'Replicate'

df_misclassification_percentage_QP$percentage_misclassified <- as.character(df_misclassification_percentage_QP$percentage_misclassified)

df_misclassification_percentage_QP_summary <- df_misclassification_percentage_QP %>%
  group_by(percentage_misclassified, category) %>% 
  dplyr::summarise(median=median(percent_accordance),
                   min=min(percent_accordance), max=max(percent_accordance), .groups = 'drop') %>% ungroup()

df_misclassification_percentage_QP_summary$category <- stringr::str_to_title(df_misclassification_percentage_QP_summary$category)

ggplot(df_misclassification_percentage_QP_summary)+
  geom_ribbon(aes(x = percentage_misclassified, ymin = min, ymax=max, fill=category, group=category), alpha=0.2)+
  geom_line(aes(y=median, percentage_misclassified, col=category, group = category), lty='dashed')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'bottom')+
  labs(y='Percentage of accordance in DA across cancer types', x='Percentage of misclassified mutations', fill=NULL, color=NULL)
##-----------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------##
df_misclassification_percentage_QP_perct <- rbind.data.frame(cbind.data.frame(melt(sapply(diagDM_misclassification_signature_nucleotide_tests_QP, function(replicate) give_agreement_in_DA_percentages_perct(replicate, resDM_nucleotides_tests, return_df=T), simplify = F), id.vars=c('percentage_misclassified', 'percent_accordance', 'ct')),
                                                                              category='nucleotides'),
                                                             cbind.data.frame(melt(sapply(diagDM_misclassification_signature_exposures_signatures_tests_QP, function(replicate) give_agreement_in_DA_percentages_perct(replicate, resDM_sigs_tests, return_df=T), simplify = F), id.vars=c('percentage_misclassified', 'percent_accordance', 'ct')),
                                                                              category='signatures'))
head(df_misclassification_percentage_QP_perct)
df_misclassification_percentage_QP_perct <- df_misclassification_percentage_QP_perct %>% group_by(ct, category, percentage_misclassified) %>%
  summarise(percent_accordance_across_replicates=mean(percent_accordance, na.rm = T))
df_misclassification_percentage_QP_perct$percentage_misclassified = percentages_misclassification[df_misclassification_percentage_QP_perct$percentage_misclassified]

ggplot(df_misclassification_percentage_QP_perct,# %>% filter(ct == 'Bone-Osteosarc'),
       aes(x=percentage_misclassified, y=percent_accordance_across_replicates, col=ct))+
  geom_line(aes(group=ct))+facet_wrap(.~category)
##-----------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------##
## select one single run across replicates, making sure they have converged
test_run <- function(run_it){
  if(is.null(run_it)){
    # find another run which did converge
    'Fail'
  }else{
    if(!is.na(wald_TMB_wrapper(run_it))){
      ## it converged. All good
      'Good'
    }else{
      'Fail'
    }
  }
}

repl = 1
diagDM_misclassification_signature_exposures_nucleotides_QP_onerepl_all_converged = diagDM_misclassification_signature_exposures_nucleotides_QP[[repl]]
for(ct in enough_samples){
  for(missclassification_it in paste0('misclassified', percentages_misclassification)){
    run_it = diagDM_misclassification_signature_exposures_nucleotides_QP_onerepl_all_converged[[ct]][[missclassification_it]]
    if(test_run(run_it) == 'Fail'){
      while( (test_run(run_it) == 'Fail') & (repl < nreplicates)){
        repl = repl+1
        run_it = diagDM_misclassification_signature_exposures_nucleotides_QP[[repl]][[ct]][[missclassification_it]]
        test_run(run_it)
      }
    }
    if((test_run(run_it) == 'Fail')){
      cat('Failed for ct ', ct, ' ', missclassification_it, ' \n')
    }
    diagDM_misclassification_signature_exposures_nucleotides_QP_onerepl_all_converged[[ct]][[missclassification_it]] = run_it
  }
}
diagDM_misclassification_signature_exposures_nucleotides_QP_onerepl_all_converged$`Bone-Osteosarc`$misclassified0.4
diagDM_misclassification_signature_exposures_nucleotides_QP_onerepl_all_converged$`CNS-Medullo`$misclassified0.4

source("../4_figures_and_tables/helper_figures_manuscript.R")

for(it_misclassification in 1:length(percentages_misclassification)){
  give_plot_betas_all_ct(sapply(diagDM_misclassification_signature_exposures_nucleotides_QP_onerepl_all_converged, `[[`, it_misclassification, simplify = F),
                         names_trinucleotide=vector_cats_to_logR(colnames(read_info_list[[1]]$dataset_nucleotidesubstitution1$Y)), prepare_for_tikz = F)+
    ggtitle(paste0('Missclassification: ', percentages_misclassification[it_misclassification]*100, '%'))
  ggsave(paste0("../../../results/results_TMB/3_real_data_assessment/3_missclassification/cors_trinucleotide3sorted_v3_missclassification_", it_misclassification, ".pdf"), height = 5, width = 5.5)
  ggsave(paste0("../../../results/results_TMB/3_real_data_assessment/3_missclassification/cors_trinucleotide3sorted_v3_missclassification_", it_misclassification, ".png"), height = 5, width = 5.5)
}

##-----------------------------------------------------------------------------------------##

par(mfrow=c(1, length(percentages_misclassification)))
for(missclassification_it in paste0('misclassified', percentages_misclassification)){
  plot(as.vector(read_info_list$`Bone-Osteosarc`$dataset_nucleotidesubstitution1$Y),
       as.vector(misclassification_signature_exposures_onerepl$`Bone-Osteosarc`[[missclassification_it]]$nucleotide_matrix$Y))
}
