##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")

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
# nonexogenous = read.table("../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t",
#                           comment.char = "#", fill = F)
# nonexogenouswSBS1SBS5 = read.table("../../data/cosmic/exogenous_signatures_SBS_withSBS1SBS5.txt", sep = "\t",
#                                    comment.char = "#", fill = F)
# sigs_cosmic0 <- read.table(paste0( "../../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
#                            stringsAsFactors = FALSE, sep = ',', header = TRUE)
# rownames(sigs_cosmic0) <- paste0(substr(sigs_cosmic0$SubType, 1, 1),'[',
#                                  sigs_cosmic0$Type, ']', substr(sigs_cosmic0$SubType, 3, 3))
# sigs_cosmic0 <- sigs_cosmic0[-c(1,2)];
# sigs_cosmic <- colnames(sigs_cosmic0)
# 
# nucleotide_colours_logR <- c('C$>$A/T$>$G'= '#3cb371', 'C$>$G/T$>$G'= '#90ee90', 'C$>$T/T$>$G'= '#66cdaa',
#                              'T$>$A/T$>$G'= '#cd5c5c', 'T$>$C/T$>$G'= '#f4a460')
# nucleotide_colours <- c('C>A' = '#3cb371', 'C>G'= '#90ee90', 'C>T'= '#66cdaa',
#                         'T>A'= '#cd5c5c', 'T>C'= '#f4a460', 'T>G'='red')
# nucleotide_colours_dollar <- c('C$>$A' = '#3cb371', 'C$>$G'= '#90ee90', 'C$>$T'= '#66cdaa',
#                                'T$>$A'= '#cd5c5c', 'T$>$C'= '#f4a460', 'T$>$G'='red')
# 
# pcawg_meta <- read.table("../../data/restricted/pcawg/repository_1567600367.tsv", sep = "\t", head=T)
# metastatic_samples <- pcawg_meta[grepl('Metastatic', pcawg_meta$Specimen.Type),]
# pcawg_meta2 <- read.table("../../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt", sep = "\t", head=T)
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
             fullRE_DMDL_SP = try(readRDS(paste0("../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/fullREDM_", ct, "_signaturesPCAWG.RDS"))),
             dataset_active_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../data/", override_warning_X_Z = T),
             dataset_nucleotidesubstitution3 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution3", path_to_data = "../../../data/", override_warning_X_Z = T))
  .x
}
read_info_list <- lapply(enough_samples, function(ct){
  read_info(ct)
}); names(read_info_list) <- enough_samples
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
diagRE_DM <- sapply(read_info_list, function(i) i$diagRE_DMDL_SP, simplify = F)
diagRE_DM_tests <-  sapply(diagRE_DM, CompSign::wald_TMB_wrapper, simplify = F)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
pcawg_palette <- pcawg.colour.palette(x = gsub("\\..*", "", names(read_info_list)),  scheme = "tumour.subtype")
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Get signatures when active signatures are modified, and run the model to compare results
## leave_one_out_exposures: list of exposures, for each cancer type, where one of each active signature has been left out
## diagDM_leave_one_out_exposures: The model run on these exposures

colnames(read_info_list[[ct]]$dataset_active_sigs$Y)
leave_one_out_exposures <- sapply(enough_samples, function(ct){
  active_sigs_in_ct <- colnames(read_info_list[[ct]]$dataset_active_sigs$Y)
  x <- lapply(active_sigs_in_ct, function(active_sig_it){
    extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                                       subset_signatures = active_sigs_in_ct[!(active_sigs_in_ct %in% active_sig_it)],
                         signature_version=NULL, signature_definition = sigs_cosmic0)
  })
  names(x) <- paste0('-', active_sigs_in_ct)
  x
}, simplify = F)

add_one_signature_exposures <- sapply(enough_samples, function(ct){
  active_sigs_in_ct <- colnames(read_info_list[[ct]]$dataset_active_sigs$Y)
  ## sample 4 signatures from sigs_cosmic0
  additional_signatures <- sample(colnames(sigs_cosmic0)[! (colnames(sigs_cosmic0) %in% active_sigs_in_ct)], 4)
  x <- lapply(additional_signatures, function(additional_signature_ct){
    extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                         subset_signatures = c(additional_signature_ct, active_sigs_in_ct),
                         signature_version=NULL, signature_definition = sigs_cosmic0)
  })
  names(x) <- paste0('+', additional_signatures)
  x
}, simplify = F)


percentages_misclassification <- c(0.05, 0.10, 0.20, 0.4)
give_missclassified_exposures_and_nucleotides <- function(){
  sapply(enough_samples, function(ct){
    ## Creating exposures with a percentage of misclassified counts
    x <- lapply(percentages_misclassification, function(percentages_misclassification_ct){
      original_matrix = read_info_list[[ct]]$dataset_nucleotidesubstitution3
      ## this simulation only works if the data are organised in a specific way
      stopifnot(all(original_matrix$x[,2] == rep(c(0,1), each=nrow(original_matrix$x)/2)))
      stopifnot(all(original_matrix$z == rbind(diag(1, nrow = nrow(original_matrix$x)/2), diag(1, nrow = nrow(original_matrix$x)/2))))
      ## for each patient, select a percentage of counts and misclassify them, i.e. put them in the subclonal row if they are clonal, and in the clonal group if they are subclonal
      counts_to_misclassify = sample(1:sum(original_matrix$Y), size = round(sum(original_matrix$Y)*percentages_misclassification_ct),
                                     replace = F)
      cumsum_matrix = cumsum(original_matrix$Y) ## this is done by column
      new_dataset_obj_trinucleotide <- original_matrix
      for(i in counts_to_misclassify){
        idx_in_Y = (which(i < cumsum_matrix)[1])
        if(is.na(idx_in_Y)){
          idx_in_Y <- length(cumsum_matrix) ## i belongs to the last entry of cumsum_matrix
        }
        ## row
        x_in_Y = idx_in_Y %% nrow(new_dataset_obj_trinucleotide$Y)
        ## column
        y_in_Y = 1 + (idx_in_Y %/% nrow(new_dataset_obj_trinucleotide$Y))
        ## instead of 0, put last row
        # cat('nrow_dataset=', nrow(original_matrix$Y), ' idx=', idx_in_Y, ' x_in_Y=', x_in_Y, ' y_in_Y=', y_in_Y, '\n')
        if(x_in_Y == 0){
          x_in_Y = nrow(new_dataset_obj_trinucleotide$Y)
          y_in_Y = y_in_Y-1
        }
        stopifnot( ((y_in_Y-1)*nrow(new_dataset_obj_trinucleotide$Y)+x_in_Y) == idx_in_Y)
        ## add it to the row of the same patient but in the other group (clonal/subclonal)
        new_x_in_Y = ifelse(x_in_Y <= nrow(new_dataset_obj_trinucleotide$Y)/2, x_in_Y+nrow(new_dataset_obj_trinucleotide$Y)/2, x_in_Y-nrow(new_dataset_obj_trinucleotide$Y)/2)
        ## add it to the row of the same patient but in the other group (clonal/subclonal). Column remains the same.
        new_dataset_obj_trinucleotide$Y[new_x_in_Y,y_in_Y] = new_dataset_obj_trinucleotide$Y[new_x_in_Y,y_in_Y] + 1
        ## subtract the count from the original row
        new_dataset_obj_trinucleotide$Y[x_in_Y,y_in_Y] = new_dataset_obj_trinucleotide$Y[x_in_Y,y_in_Y] - 1
      }
      new_nucleotide_matrix = new_dataset_obj_trinucleotide
      new_nucleotide_matrix$Y = sapply(sort(unique(sapply(colnames(new_nucleotide_matrix$Y), function(i) substr(i, 3, 5)))), function(SBS1_changes){
        rowSums(new_nucleotide_matrix$Y[,grepl(SBS1_changes, colnames(new_nucleotide_matrix$Y))])
      })
      
      ## total number of mutations per patient
      stopifnot(rowSums(matrix(rowSums(original_matrix$Y), ncol=2)) == rowSums(matrix(rowSums(new_nucleotide_matrix$Y), ncol=2)))
      
      list(signatures_matrix=extract_sigs_TMB_obj(dataset_obj_trinucleotide=new_dataset_obj_trinucleotide,
                                                  subset_signatures = colnames(read_info_list[[ct]]$dataset_active_sigs$Y),
                                                  signature_version=NULL, signature_definition = sigs_cosmic0),
           nucleotide_matrix=new_nucleotide_matrix)
    })
    names(x) <- paste0('misclassified', percentages_misclassification)
    x
  }, simplify = F)
}
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Run model

diagDM_leave_one_out_exposures <- sapply(enough_samples, function(ct){
  sapply(leave_one_out_exposures[[ct]], function(leave_one_out_it){
    CompSign::wrapper_run_TMB(object = leave_one_out_it,
                              model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
  }, simplify = F)
}, simplify = F)

saveRDS(diagDM_leave_one_out_exposures,
        file = "../../../data/assessing_models_real_data/inference_results/TMB/leave_one_out_diagREDM.RDS")
  
# ct='Panc-AdenoCA'
# leave_one_out_it='SBS1'
# diagDM_leave_one_out_exposures[[ct]][[paste0('-', leave_one_out_it)]]
# extra_run = CompSign::wrapper_run_TMB(object = leave_one_out_exposures[[ct]][[paste0('-', leave_one_out_it)]],
#                                       model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
# extra_run
# diagDM_leave_one_out_exposures[[ct]][[paste0('-', leave_one_out_it)]] <- extra_run
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
misclassification_signature_exposures_signatures
misclassification_signature_exposures_nucleotide

give_digRE_results <- function(misclassification_signature_exposures_list){
  sapply(enough_samples, function(ct){
    sapply(misclassification_signature_exposures_list[[ct]], function(leave_one_out_it){
      CompSign::wrapper_run_TMB(object = leave_one_out_it,
                                model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
    }, simplify = F)
  }, simplify = F)
}

misclassification_signature_exposures <- list()
misclassification_signature_exposures_signatures <- list()
misclassification_signature_exposures_nucleotide <- list()
diagDM_misclassification_signature_exposures_signatures <- list()
diagDM_misclassification_signature_exposures_nucleotide <- list()
# for(repl in 1:4){
for(repl in 2:4){
  cat('Replicate: ', repl, '\n')
  misclassification_signature_exposures[[repl]] <- give_missclassified_exposures_and_nucleotides()
  misclassification_signature_exposures_signatures[[repl]] = sapply(misclassification_signature_exposures[[repl]], function(i) sapply(i, `[[`, 'signatures_matrix', simplify = F), simplify = F)
  misclassification_signature_exposures_nucleotide[[repl]] = sapply(misclassification_signature_exposures[[repl]], function(i) sapply(i, `[[`, 'nucleotide_matrix', simplify = F), simplify = F)
  diagDM_misclassification_signature_exposures_signatures[[repl]] <- give_digRE_results(misclassification_signature_exposures_signatures[[repl]])
  diagDM_misclassification_signature_exposures_nucleotide[[repl]] <- give_digRE_results(misclassification_signature_exposures_nucleotide[[repl]])
  
  ## save dataset
  saveRDS(misclassification_signature_exposures[[repl]],
          file = paste0("../../../data/assessing_models_real_data/simulated_datasets/misclassification_signature_exposures_repl", repl, ".RDS"))
  
  ## save inference results
  saveRDS(diagDM_misclassification_signature_exposures_nucleotide[[repl]],
          file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/diagDM_misclassification_signature_exposures_nucleotide_repl", repl, ".RDS"))
  saveRDS(diagDM_misclassification_signature_exposures_signatures[[repl]],
          file = paste0("../../../data/assessing_models_real_data/inference_results/TMB/diagDM_misclassification_signature_exposures_signatures_repl", repl, ".RDS"))
  
  ## Check that indeed, counts have only been re-arranged within patients
  stopifnot(rowSums(matrix(rowSums(misclassification_signature_exposures_nucleotide[[repl]]$`Bone-Osteosarc`$misclassified0.05$Y), ncol=2)) ==
              rowSums(matrix(rowSums(misclassification_signature_exposures_nucleotide[[repl]]$`Bone-Osteosarc`$misclassified0.1$Y), ncol=2)))
  
}

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
diagDM_add_one_signature_exposures <- sapply(enough_samples, function(ct){
  sapply(add_one_signature_exposures[[ct]], function(add_one_it){
    CompSign::wrapper_run_TMB(object = add_one_it,
                              model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
  }, simplify = F)
}, simplify = F)
saveRDS(diagDM_add_one_signature_exposures,
        file = "../../../data/assessing_models_real_data/inference_results/TMB/add_one_signature_diagREDM.RDS")

ct='Skin-Melanoma.cutaneous'
add_one_it='SBS45'
diagDM_add_one_signature_exposures[[ct]][[paste0('+', add_one_it)]]
extra_run = CompSign::wrapper_run_TMB(object = add_one_signature_exposures[[ct]][[paste0('+', add_one_it)]],
                                      model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
extra_run
# diagDM_add_one_signature_exposures[[ct]][[paste0('+', add_one_it)]] <- extra_run

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Tests for DA
diagDM_leave_one_out_exposures_tests <- sapply(diagDM_leave_one_out_exposures, function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE)


diagDM_add_one_signature_exposures_tests <- sapply(diagDM_add_one_signature_exposures, function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE)

diagDM_leave_one_out_exposures_tests <- sapply(diagDM_leave_one_out_exposures, function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE)


## Several replicates for the below
diagDM_misclassification_signature_exposures_signatures_tests <- lapply(diagDM_misclassification_signature_exposures_signatures, function(j) sapply(j, function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE))

diagDM_misclassification_signature_exposures_nucleotide_tests <- lapply(diagDM_misclassification_signature_exposures_nucleotide, function(j) sapply(j, function(ct){
  sapply(ct, function(j) CompSign::wald_TMB_wrapper(j))
}, simplify = FALSE))


##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## See if the DA results are in agreement with the original (non-modified) dataset
get_pcawg_name <- function(L1){
  return(tolower(make.names(gsub("[.].*", "", L1))))
}

give_barplot_agreement_in_DA <- function(diagRE_DM_tests, diagDM_leave_one_out_exposures_tests, ylabel='Number of leave-one-out datasets'){
  diagDM_leave_one_out_exposures_tests_accordance_DA = melt(sapply(enough_samples, function(ct) t(t(rep(diagRE_DM_tests[[ct]] <= 0.05, length(diagDM_leave_one_out_exposures_tests[[ct]])) + (diagDM_leave_one_out_exposures_tests[[ct]] <= 0.05))), simplify = F))
  diagDM_leave_one_out_exposures_tests_accordance_DA$accordance = ifelse(diagDM_leave_one_out_exposures_tests_accordance_DA$value == 1, 'Discordance', 'Accordance')
  
  ggplot(diagDM_leave_one_out_exposures_tests_accordance_DA, aes(x=L1, #alpha=(accordance== 'Accordance'),
                                                                 fill=ifelse(accordance== 'Accordance',
                                                                             get_pcawg_name(L1),  'white'), colour='1'))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    scale_fill_manual(values = c(pcawg_palette, "white"="white"), na.value="cyan")+
    scale_color_manual(values = c('1'='black'))+guides(fill='none', color='none')+
    labs(alpha='Accordance', x='Cancer types', y=ylabel)
}

give_indices_beta <- function(j, items_per_category=2){
  c(1+(j-1)*2, 2+(j-1)*2)
}

give_beta_cor <- function(df, diagRE_DM, mode='leave_one_out_all_betas'){
  
  if( grepl('leave_one_out', mode) | grepl('add_one', mode)){
    ## correlation can only be computed for the first n-1 signatures, as if we leave out the last category we don't have a shared baseline
    
    ## this is with both betas, intercept and slope
    sapply(enough_samples, function(ct) sapply(1:length(df[[ct]]),
         function(j){
           if(j == length(df[[ct]])){
             NA ## this is the last signature; it doesn't make any sense for beta coefs to be compared because we are using a different baseline
           }else{
             if(mode == 'leave_one_out_all_betas'){
              cor(plot_betas(df[[ct]][[j]], return_df = T, plot = F)[,'Estimate'],
                                                                 plot_betas(diagRE_DM[[ct]], return_df = T, plot = F)[,'Estimate'][-give_indices_beta(j, 2)])
             }else if(mode == 'leave_one_out_beta_slope'){
               ## select only even coefficients
               cor(plot_betas(df[[ct]][[j]], return_df = T, plot = F)[,'Estimate'][c(F,T)],
                   plot_betas(diagRE_DM[[ct]], return_df = T, plot = F)[,'Estimate'][-give_indices_beta(j, 2)][c(F,T)])
             }else if(mode == 'add_one'){
               cor(plot_betas(df[[ct]][[j]], return_df = T, plot = F)[,'Estimate'],
                   plot_betas(diagRE_DM[[ct]], return_df = T, plot = F)[,'Estimate'][-c(1:2)]) ## removing the first (new) signature
             }else if(mode == 'add_one_beta_slope'){
               cor(plot_betas(df[[ct]][[j]], return_df = T, plot = F)[,'Estimate'][c(F,T)],
                   plot_betas(diagRE_DM[[ct]], return_df = T, plot = F)[,'Estimate'][-c(1:2)][c(F,T)])
             }
             
           }
           }, simplify = F),
           simplify = F)
  }else{
    stop('Not implemented')
  }
}

plt1 <- give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_leave_one_out_exposures_tests, ylabel='Number of leave-one-out datasets')
plt2 <- give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_add_one_signature_exposures_tests, ylabel='Number of add-one datasets')
plt3 <- give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_misclassification_signature_exposures_nucleotide_tests, ylabel='Number of datasets with misassigned mutations')
plt4 <- give_barplot_agreement_in_DA(diagRE_DM_tests, diagDM_misclassification_signature_exposures_signatures_tests, ylabel='Number of datasets with misassigned mutations')

cowplot::plot_grid(plt1, plt2, nrow=1)

sapply(diagDM_misclassification_signature_exposures_signatures_tests, function(j){
  give_barplot_agreement_in_DA(diagRE_DM_tests, j, ylabel='Number of datasets with misassigned mutations')
}, simplify = F)


do.call('grid.arrange', c(nrow=1, lapply(diagDM_misclassification_signature_exposures_signatures_tests, function(j) give_barplot_agreement_in_DA(j, diagRE_DM_tests, ylabel='Number of datasets with misassigned mutations'))))
do.call('grid.arrange', c(nrow=1, lapply(diagDM_misclassification_signature_exposures_nucleotide_tests, function(j) give_barplot_agreement_in_DA(j, diagRE_DM_tests, ylabel='Number of datasets with misassigned mutations'))))

give_agreement_in_DA_percentages <- function(diagRE_DM_tests, df, return_df=F){
  diagDM_leave_one_out_exposures_tests_accordance_DA = melt(sapply(enough_samples, function(ct) t(t(rep(diagRE_DM_tests[[ct]] <= 0.05, length(df[[ct]])) + (df[[ct]] <= 0.05))), simplify = F))
  diagDM_leave_one_out_exposures_tests_accordance_DA$accordance = ifelse(diagDM_leave_one_out_exposures_tests_accordance_DA$value == 1, 'Discordance', 'Accordance')
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



cowplot::plot_grid(give_agreement_in_DA_percentages(diagRE_DM_tests, diagDM_misclassification_signature_exposures_nucleotide_tests),
                   give_agreement_in_DA_percentages(diagRE_DM_tests, diagDM_misclassification_signature_exposures_signatures_tests))



df_misclassification_percentage <- rbind.data.frame(cbind.data.frame(melt(sapply(diagDM_misclassification_signature_exposures_nucleotide_tests, function(j) give_agreement_in_DA_percentages(j, diagRE_DM_tests, return_df=T), simplify = F), id.vars=c('percentage_misclassified', 'percent_accordance')),
                                  category='nucleotides'),
                 cbind.data.frame(melt(sapply(diagDM_misclassification_signature_exposures_signatures_tests, function(j) give_agreement_in_DA_percentages(j, diagRE_DM_tests, return_df=T), simplify = F), id.vars=c('percentage_misclassified', 'percent_accordance')),
                                  category='signatures'))
df_misclassification_percentage$percentage_misclassified = percentages_misclassification[df_misclassification_percentage$percentage_misclassified]
colnames(df_misclassification_percentage)[colnames(df_misclassification_percentage) == 'L1'] = 'Replicate'

ggplot(df_misclassification_percentage,
       aes(x=percentage_misclassified, y = percent_accordance, col=category))+
  geom_line(data = df_misclassification_percentage %>% group_by(percentage_misclassified, category) %>% 
              summarise(median=median(percent_accordance)) %>% ungroup(), aes(y=median))+
  geom_boxplot(aes(group=interaction(percentage_misclassified, category)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(y='Percentage of accordance across cancer types', x='Percentage of misclassified mutations')

##-----------------------------------------------------------------------------------------------------##
## Beta correlation
theme_set(theme_bw())
table(is.na(beta_cor_leave_one_out_exposures$value))

plot_beta_cor <- function(title_arg, ...){
  beta_cor_leave_one_out_exposures <- give_beta_cor(...)
  beta_cor_leave_one_out_exposures <- melt(beta_cor_leave_one_out_exposures)
  ## remove NAs
  beta_cor_leave_one_out_exposures <- beta_cor_leave_one_out_exposures[!is.na(beta_cor_leave_one_out_exposures$value),]
  beta_cor_leave_one_out_exposures <- beta_cor_leave_one_out_exposures %>% group_by(L1) %>% arrange(-value, .by_group = TRUE)
  beta_cor_leave_one_out_exposures$x = paste0(beta_cor_leave_one_out_exposures$L2, beta_cor_leave_one_out_exposures$L1)
  beta_cor_leave_one_out_exposures$x <- factor(beta_cor_leave_one_out_exposures$x, levels=beta_cor_leave_one_out_exposures$x)
  ggplot(beta_cor_leave_one_out_exposures, aes(x=x, y=value, col=get_pcawg_name(L1)))+geom_point()+geom_line(aes(group=L1))+
    # facet_wrap(.~L1, scales = 'free_x', drop = T, nrow=1)+
    scale_color_manual(values = c(pcawg_palette, "white"="white"), na.value="cyan")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    labs(y='Correlation', x=title_arg)
}

beta_cor_leave_one_out_all_betas <-  plot_beta_cor(diagDM_leave_one_out_exposures, diagRE_DM, mode='leave_one_out_all_betas',
                                                   title_arg=latex2exp::TeX(r"(Datasets with one signature removed prior to parameter estimation ($\beta_0, \beta_1$))"))
beta_cor_leave_one_out_beta_slope <-  plot_beta_cor(diagDM_leave_one_out_exposures, diagRE_DM, mode='leave_one_out_beta_slope',
                                                    title_arg=latex2exp::TeX(r"(Datasets with one signature removed prior to parameter estimation ($\beta_1$))"))
beta_cor_add_one <-  plot_beta_cor(diagDM_leave_one_out_exposures, diagRE_DM, mode='add_one',
                                              title_arg=latex2exp::TeX(r"(Datasets with one signature added prior to parameter estimation ($\beta_0, \beta_1$))"))
beta_cor_add_one_beta_slope <-  plot_beta_cor(diagDM_leave_one_out_exposures, diagRE_DM, mode='add_one_beta_slope',
                                                    title_arg=latex2exp::TeX(r"(Datasets with one signature added prior to parameter estimation ($\beta_1$))"))

cowplot::plot_grid(beta_cor_leave_one_out_all_betas,
                   beta_cor_leave_one_out_beta_slope,
                   beta_cor_add_one,
                   beta_cor_add_one_beta_slope)
