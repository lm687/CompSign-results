#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
library(ggplot2)
require(R.utils)
require(dplyr)
library(parallel)
library(RColorBrewer)
library(jcolors)
library(reshape2)
source("mm_multinomial/helper_functions.R")
source("helper_TMB.R")
# set.seed(1234)
re_run_inference = FALSE ## use cache or not
give_summary_runs = FALSE ## whether to run the section to see what has converged, what hasn't, etc.
folder_robjs = "../../data/pcawg_robjects_cache/tmb_results/"
folder_robjs_nlminb = "../../data/pcawg_robjects_cache/tmb_results/nlminb/"
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
# TMB::compile("mm_multinomial/ME_multinomial.cpp")
# dyn.load(dynlib("mm_multinomial/ME_multinomial"))
# TMB::compile("mm_multinomial/ME_dirichletmultinomial.cpp", "-std=gnu++17")
# dyn.load(dynlib("mm_multinomial/ME_dirichletmultinomial"))
# TMB::compile("mm_multinomial/ME_LNM.cpp", "-std=gnu++17")
# dyn.load(dynlib("mm_multinomial/ME_LNM"))
TMB::compile("mm_multinomial/fullRE_ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_multinomial"))
TMB::compile("mm_multinomial/fullRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_dirichletmultinomial"))
# TMB::compile("mm_multinomial/fullRE_ME_dirichletmultinomial_altpar.cpp", "-std=gnu++17")
# dyn.load(dynlib("mm_multinomial/fullRE_ME_dirichletmultinomial_altpar"))
TMB::compile("mm_multinomial/fullRE_ME_singlelambda_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_singlelambda_dirichletmultinomial"))
TMB::compile("mm_multinomial/diagRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/diagRE_ME_dirichletmultinomial"))
TMB::compile("mm_multinomial/diagRE_ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/diagRE_ME_multinomial"))
# TMB::compile("mm_multinomial/fullRE_ME_multinomial_categorical.cpp", "-std=gnu++17")
# dyn.load(dynlib("mm_multinomial/fullRE_ME_multinomial_categorical"))
# TMB::compile("mm_multinomial/fullRE_ME_dirichletmultinomial_categorical.cpp", "-std=gnu++17")
# dyn.load(dynlib("mm_multinomial/fullRE_ME_dirichletmultinomial_categorical"))
TMB::compile("mm_multinomial/diagRE_dirichletmultinomial_single_lambda.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/diagRE_dirichletmultinomial_single_lambda"))
TMB::compile("mm_multinomial/fullRE_dirichletmultinomial_single_lambda.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_dirichletmultinomial_single_lambda"))
TMB::compile("mm_multinomial/diagRE_dirichletmultinomial_single_lambda.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/diagRE_dirichletmultinomial_single_lambda"))
TMB::compile("mm_multinomial/fullRE_dirichletmultinomial_single_lambda2.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_dirichletmultinomial_single_lambda2"))

#-------------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------------#
samples_files = data.frame(do.call('rbind', sapply(gsub("_ROO.RDS", "", list.files("../../data/roo/")),
                                        strsplit, split = "_")))
colnames(samples_files) = c('CT', 'type')
# table(samples_files[,1], samples_files[,2])
# ct = "Bladder-TCC" #samples_files[1,1]
# typedata =nucleotidesubstitution3  #"signatures" #samples_files[1,2]

samples_files2 = samples_files %>% filter(type != "nucleotidesubstitution3")
rownames(samples_files2) = rownames(samples_files)[samples_files$type != "nucleotidesubstitution3"]
#-------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------#
## run at random
if(re_run_inference){
  
  enough_samples = readLines("~/Desktop/CT_sufficient_samples.txt")

  mclapply(sample(which(is.na(match(rownames(samples_files2),
                                    gsub(".RDS", "", gsub("M_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
   function(idx){
     i = samples_files2[idx,]
     x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "M"),
                     timeout = 300, onTimeout = "warning")
     saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "M_", rownames(i), ".RDS"))
  })
  

  mclapply(sample(which(is.na(match(as.vector(sapply(enough_samples, function(i) c(paste0(i, '_signatures'), paste0(i, "_nucleotidesubstitution1")))),
                                    gsub(".RDS", "", gsub("fullRE_DMSL_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
           function(idx){
             i = samples_files2[idx,]
             x = withTimeout(wrapper_run_TMB(load_PCAWG(ct = i[1,1], typedata = i[1,2]), sort_columns = T, smart_init_vals = T,
                                             model = "fullREDMsinglelambda"),
                             timeout = 300, onTimeout = "warning")
             saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "fullRE_DMSL_", rownames(i), ".RDS"))
           })
  
  mclapply(sample(which(is.na(match(as.vector(sapply(enough_samples, function(i) c(paste0(i, '_signatures')))),
                                    gsub(".RDS", "", gsub("fullRE_DMSL2_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
           function(idx){
             i = samples_files2[idx,]
             x = withTimeout(wrapper_run_TMB(sort_columns_TMB(load_PCAWG(ct = i[1,1], typedata = i[1,2])), smart_init_vals = T,
                                             model = "fullREDMsinglelambda2"),
                             timeout = 300, onTimeout = "warning")
             saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "fullRE_DMSL2_", rownames(i), ".RDS"))
           })
  
  mclapply(sample(which(is.na(match(as.vector(sapply(enough_samples, function(i) c(paste0(i, '_signatures'), paste0(i, "_nucleotidesubstitution1")))),
                                    gsub(".RDS", "", gsub("diagRE_DMSL_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
           function(idx){
             i = samples_files2[idx,]
             x = wrapper_run_TMB(load_PCAWG(ct = i[1,1], typedata = i[1,2]), sort_columns = T, smart_init_vals = T,
                                             model = "diagREDMsinglelambda")
             saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "diagRE_DMSL_", rownames(i), ".RDS"))
           })
  
  # mclapply(sample(which(is.na(match(rownames(samples_files2),
  #                                   gsub(".RDS", "", gsub("DM_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
  #          function(idx){
  #            i = samples_files2[idx,]
  #            x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "DM"),
  #                            timeout = 300, onTimeout = "warning")
  #            saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "DM_", rownames(i), ".RDS"))
  #          })
  
  
  # mclapply(sample(which(is.na(match(rownames(samples_files2),
  #                                   gsub(".RDS", "", gsub("LNM_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
  #          function(idx){
  #            # mclapply(1:nrow(samples_files), function(idx){
  #            i = samples_files2[idx,]
  #            x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "LNM"),
  #                            timeout = 300, onTimeout = "warning")
  #            saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "LNM_", rownames(i), ".RDS"))
  #          })

  mclapply(sample(which(is.na(match(rownames(samples_files2),
                                    gsub(".RDS", "", gsub("fullRE_M_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
           function(idx){
             i = samples_files2[idx,]
             x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "fullRE_M"),
                             timeout = 300, onTimeout = "warning")
             saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "fullRE_M_", rownames(i), ".RDS"))
           })
  
  mclapply(sample(which(is.na(match(rownames(samples_files2),
                                     gsub(".RDS", "", gsub("fullRE_DM_altpar_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
            function(idx){
              i = samples_files2[idx,]
              x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "fullRE_DM_altpar"),
                              timeout = 300, onTimeout = "warning")
              saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "fullRE_DM_altpar_", rownames(i), ".RDS"))
            })

  ## diagonal M
  mclapply(sample(which(is.na(match(rownames(samples_files2),
                                    gsub(".RDS", "", gsub("diagRE_M_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/optim/"))))))),
           function(idx){
             outcome_inference="Not good"
             counter_tries = 0
             while(outcome_inference != "Good" & counter_tries < 6){
               i = samples_files2[idx,]
               x = withTimeout(wrapper_run_TMB(load_PCAWG(i[1,1], i[1,2]), model = "diagRE_M"),
                               timeout = 300, onTimeout = "warning")
               outcome_inference = give_summary_per_sample(x)
               counter_tries = counter_tries + 1
             }
             saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/optim/", "diagRE_M_", rownames(i), ".RDS"))
  })
  
  ## diagonal DM
  mclapply(sample(which(is.na(match(rownames(samples_files2),
                                    gsub(".RDS", "", gsub("diagRE_DM_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
           function(idx){
             outcome_inference="Not good"
             counter_tries = 0
             while(outcome_inference != "Good" & counter_tries < 6){
               i = samples_files2[idx,]
               x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "diagRE_DM"),
                               timeout = 300, onTimeout = "warning")
               outcome_inference = give_summary_per_sample(x)
               counter_tries = counter_tries + 1
             }
             saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "diagRE_DM_", rownames(i), ".RDS"))
           })
  
  
  ## Categorical: not over-dispersed
  # mclapply(sample(which(is.na(match(rownames(samples_files2),
  #                                   gsub(".RDS", "", gsub("fullRE_Mcat_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
  #          function(idx){
  #            i = samples_files2[idx,]
  #            x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "fullRE_Mcat"),
  #                            timeout = 300, onTimeout = "warning")
  #            saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "fullRE_Mcat_", rownames(i), ".RDS"))
  #          })
  # 
  # ## Categorical: over-dispersed
  # mclapply(sample(which(is.na(match(rownames(samples_files2),
  #                                   gsub(".RDS", "", gsub("fullRE_DMcat_", "", list.files("../../data/pcawg_robjects_cache/tmb_results/"))))))),
  #          function(idx){
  #            i = samples_files2[idx,]
  #            x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "fullRE_DMcat"),
  #                            timeout = 300, onTimeout = "warning")
  #            saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/", "fullRE_DMcat_", rownames(i), ".RDS"))
  #          })
  
}

#----------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------------------------------#
# re-run those that had NaN due to bad initial values
# ct = "Lymph-BNHL"
# type=  "signatures" # "nucleotidesubstitution1
# model = "diagRE_DM"
# rm(x)
# x = wrapper_run_TMB(ct, type, model, allow_new_LNM = TRUE, sort_columns = T)
# x
# # load_PCAWG(ct, type)$Y
# saveRDS(object = x, file=paste0(folder_robjs, model, "_", ct, "_", type, ".RDS"))
#------------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#
# results_TMB_M = lapply( python_like_select(list.files(folder_robjs), "^M_"), function(i) readRDS(paste0(folder_robjs, i)))
# names(results_TMB_M) = sapply(python_like_select(list.files(folder_robjs), "^M_"), clean_name)
# 
# results_TMB_DM = lapply( python_like_select(list.files(folder_robjs), "^DM_"), function(i) readRDS(paste0(folder_robjs, i)))
# names(results_TMB_DM) = sapply(python_like_select(list.files(folder_robjs), "^DM_"), clean_name)

# results_TMB_DM_dep = lapply( python_like_select(list.files("../../data/robjects_cache/tmb_results_dep/"), "^DM_"),
#                              function(i) readRDS(paste0("../../data/robjects_cache/tmb_results_dep/", i)))
# names(results_TMB_DM_dep) = sapply(python_like_select(list.files("../../data/robjects_cache/tmb_results_dep/"), "^DM_"), clean_name)

# results_TMB_LNM = lapply( python_like_select(list.files(folder_robjs), "^LNM_"), function(i) readRDS(paste0(folder_robjs, i)))
# names(results_TMB_LNM) = sapply(python_like_select(list.files(folder_robjs), "^LNM_"), clean_name_fullRE)

results_TMB_fullRE_M = lapply( python_like_select(list.files(folder_robjs), "^fullRE_M_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_M) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_M_"), clean_name_fullRE)

full_RE_DM = python_like_select(list.files(folder_robjs), "^fullRE_DM_"); full_RE_DM = full_RE_DM[-grep("_altpar_", full_RE_DM)]
results_TMB_fullRE_DM = lapply( full_RE_DM, function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_DM) = sapply(full_RE_DM, clean_name_fullRE)

results_TMB_diagRE_DM = lapply( python_like_select(list.files(folder_robjs), "^diagRE_DM_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_diagRE_DM) = sapply(python_like_select(list.files(folder_robjs), "^diagRE_DM_"), clean_name_fullRE)

results_TMB_diagRE_M = lapply( python_like_select(list.files(folder_robjs), "^diagRE_M_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_diagRE_M) = sapply(python_like_select(list.files(folder_robjs), "^diagRE_M_"), clean_name_fullRE)

results_TMB_fullRE_DMSL = lapply( python_like_select(list.files(folder_robjs), "^fullRE_DMSL_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_DMSL) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_DMSL_"), clean_name_fullRE)

results_TMB_diagRE_DMSL = lapply( python_like_select(list.files(folder_robjs), "^diagRE_DMSL_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_diagRE_DMSL) = sapply(python_like_select(list.files(folder_robjs), "^diagRE_DMSL_"), clean_name_fullRE)

results_TMB_fullRE_DM_nlminb = lapply( python_like_select(list.files(folder_robjs_nlminb), "^fullRE_DM_"), function(i) readRDS(paste0(folder_robjs_nlminb, i)))
names(results_TMB_fullRE_DM_nlminb) = sapply(python_like_select(list.files(folder_robjs_nlminb), "^fullRE_DM_"), clean_name_fullRE)

results_TMB_fullRE_DMSL_nlminb = lapply( python_like_select(list.files(folder_robjs_nlminb), "^fullRE_DMSL_"), function(i) readRDS(paste0(folder_robjs_nlminb, i)))
names(results_TMB_fullRE_DMSL_nlminb) = sapply(python_like_select(list.files(folder_robjs_nlminb), "^fullRE_DMSL_"), clean_name_fullRE)

results_TMB_diagRE_DMSL_nlminb = lapply( python_like_select(list.files(folder_robjs_nlminb), "^diagRE_DMSL_"), function(i) readRDS(paste0(folder_robjs_nlminb, i)))
names(results_TMB_diagRE_DMSL_nlminb) = sapply(python_like_select(list.files(folder_robjs_nlminb), "^diagRE_DMSL_"), clean_name_fullRE)


# results_TMB_fullRE_Mcat = lapply( python_like_select(list.files(folder_robjs), "^fullRE_Mcat_"), function(i) readRDS(paste0(folder_robjs, i)))
# names(results_TMB_fullRE_Mcat) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_Mcat_"), clean_name_fullRE)
# 
# results_TMB_fullRE_DMcat = lapply( python_like_select(list.files(folder_robjs), "^fullRE_DMcat_"), function(i) readRDS(paste0(folder_robjs, i)))
# names(results_TMB_fullRE_DMcat) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_DMcat_"), clean_name_fullRE)

# results_TMB_fullRE_DM = lapply( python_like_select(list.files(folder_robjs), "^fullRE_DM_altpar_"), function(i) readRDS(paste0(folder_robjs, i)))
# names(results_TMB_fullRE_DM) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_DM_altpar_"), clean_name_fullRE_2)
#----------------------------------------------------------------------------------------------------#

if(give_summary_runs){
  #----------------------------------------------------------------------------------------------------#
  ## checking how many errors were there
  ## a) timeout
  ## b) non-positive-definite hessian
  ## c) good convergence
  
  give_summary_of_runs(results_TMB_LNM, long_return = FALSE)
  give_summary_of_runs(results_TMB_M, long_return = FALSE)
  give_summary_of_runs(results_TMB_DM, long_return = FALSE)
  
  # long and short version of summaries
  sapply(give_summary_of_runs(results_TMB_LNM, long_return = TRUE), length)
  give_summary_of_runs(results_TMB_M, long_return = TRUE)
  
  results_TMB_M[['Biliary-AdenoCAnucleotidesubstitution1']]
  results_TMB_DM[['Biliary-AdenoCAnucleotidesubstitution1']]
  results_TMB_LNM[['Biliary-AdenoCAnucleotidesubstitution1']]
  
  ## for creating spreadsheet
  for(i in sort(unique(samples_files$CT))){cat(i, '\n')}
  
  ## see if nuc1 is positive semi-definite or has some problem
  for(i in sort(unique(samples_files$CT))){cat(i, '\n')}
  
  length(unique(samples_files$CT))
  length(results_TMB_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))])
  #------------------------------------------------------------------------------------------------------------------------#
  
  #------------------------------------------------------------------------------------------------------------------------#
  # for(i in sapply(results_TMB_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
  #        give_summary_per_sample)){cat(i,'\n')}
  # 
  # for(i in sapply(results_TMB_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution3", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  # 
  # for(i in sapply(results_TMB_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_fullRE_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_fullRE_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution3", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_fullRE_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  #------------------------------------------------------------------------------------------------------------------------#
  
  #------------------------------------------------------------------------------------------------------------------------#
  # for(i in sapply(results_TMB_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  # 
  # for(i in sapply(results_TMB_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution3", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  # 
  # for(i in sapply(results_TMB_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  
  # for(i in sapply(results_TMB_fullRE_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  # 
  # for(i in sapply(results_TMB_fullRE_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution3", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  # 
  # for(i in sapply(results_TMB_fullRE_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}

  for(i in sapply(results_TMB_fullRE_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_fullRE_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}

  ## uncorrelated RE
  
  for(i in sapply(results_TMB_diagRE_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_diagRE_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  
  for(i in sapply(results_TMB_diagRE_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_diagRE_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  #------------------------------------------------------------------------------------------------------------------------#
  
  #------------------------------------------------------------------------------------------------------------------------#
  # for(i in sapply(results_TMB_LNM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  # 
  # for(i in sapply(results_TMB_LNM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution3", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  # 
  # for(i in sapply(results_TMB_LNM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
  #                 give_summary_per_sample)){cat(i,'\n')}
  #------------------------------------------------------------------------------------------------------------------------#

  
  #------------------------------------------------------------------------------------------------------------------------#
  stan_results = read.table("../inference_diagnostics20200828.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
  
  for(i in get_summary_stan("M", "nucleotidesubstitution1")){cat(i,'\n')}
  for(i in get_summary_stan("M", "nucleotidesubstitution3")){cat(i,'\n')}
  for(i in get_summary_stan("M", "signatures")){cat(i,'\n')}
  
  for(i in get_summary_stan("DM", "nucleotidesubstitution1")){cat(i,'\n')}
  for(i in get_summary_stan("DM", "nucleotidesubstitution3")){cat(i,'\n')}
  for(i in get_summary_stan("DM", "signatures")){cat(i,'\n')}
  
  for(i in get_summary_stan("LNM", "nucleotidesubstitution1")){cat(i,'\n')}
  for(i in get_summary_stan("LNM", "nucleotidesubstitution3")){cat(i,'\n')}
  for(i in get_summary_stan("LNM", "signatures")){cat(i,'\n')}

  ## showing that I was not doing the parallelising correctly
  plot(stan_results$total.time..s., stan_results$real.time..s.)
  plot(density(na.omit(stan_results$real.time..s./stan_results$total.time..s.)))

  ## see if there is one run for each combination, and no more
  sort(table(apply(cbind(stan_results$CT, stan_results$features, stan_results$model), 1, paste0, collapse="")))
  ## there should be only one of each
  
  
  all_runs_summary = lapply(list(results_TMB_M, results_TMB_diagRE_M, results_TMB_fullRE_M, results_TMB_fullRE_Mcat,
                                 results_TMB_DM, results_TMB_diagRE_DM, results_TMB_fullRE_DM, results_TMB_fullRE_DMcat), 
                            give_summary_of_runs2, long_return = T)
  names(all_runs_summary) = c('results_TMB_M', 'results_TMB_diagRE_M', 'results_TMB_fullRE_M', 'results_TMB_fullRE_Mcat',
                                 'results_TMB_DM', 'results_TMB_diagRE_DM', 'results_TMB_fullRE_DM', 'results_TMB_fullRE_DMcat')
  all_runs_summary_melt = melt(all_runs_summary)
  all_runs_summary_melt$ct = sapply(as.character(all_runs_summary_melt$value), function(i) gsub("nucleotidesubstitution1", "", gsub("signatures", "", i)))
  all_runs_summary_melt$type = sapply(1:length(all_runs_summary_melt$ct), function(j) gsub(all_runs_summary_melt$ct[j], "", all_runs_summary_melt$value[j]))
  all_runs_summary_melt$L1 = factor(all_runs_summary_melt$L1, levels = (c('results_TMB_M', 'results_TMB_diagRE_M', 'results_TMB_fullRE_M', 'results_TMB_fullRE_Mcat',
                                          'results_TMB_DM', 'results_TMB_diagRE_DM', 'results_TMB_fullRE_DM', 'results_TMB_fullRE_DMcat')))
  ggplot(all_runs_summary_melt, aes(group=value, y=L2))+geom_bar(position = 'stack')+facet_wrap(.~L1)
  ggplot(all_runs_summary_melt, aes(x=interaction(ct, type), y=L2, col=type))+geom_point()+facet_wrap(.~L1)
  ggplot(all_runs_summary_melt, aes(y=interaction(type, ct), x=L1, fill=L2))+geom_tile()+
    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "bottom", legend.title = element_blank())+
    scale_fill_manual(values = c("red", "#aef7b3", "#f7baae"))+guides(fill=guide_legend(nrow=2))
  ggsave("../../results/results_TMB/pcawg/summary_runs.png", height = 12, width = 7)
  
}

#------------------------------------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------------------------------------#

betasM = load_posteriors("../../data/pcawg_robjects_cache/betas91ecb3fe-4ff0-4e91-b6f0-a2eaf027f91e_M_signatures.Rdata")
betasDM = load_posteriors("../../data/pcawg_robjects_cache/betas316eb9a5-31f9-4d4b-be88-1b0e5c184286_DM_signatures.Rdata")

# some example
betasDM$posteriors_betas$`Biliary-AdenoCA_signatures_20000_DMROO`
results_TMB_DM$`Biliary-AdenoCAsignatures`
set.seed(182)
#------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------------#
## they should have the same dimensions
length(select_slope_2(results_TMB_DM$`Biliary-AdenoCAsignatures`$par.fixed[grep("beta", names(results_TMB_DM$`Biliary-AdenoCAsignatures`$par.fixed))], verbatim = FALSE))
dim(betasDM$posteriors_betas$`Biliary-AdenoCA_signatures_20000_DMROO`[,python_like_select(colnames(betasDM$posteriors_betas$`Biliary-AdenoCA_signatures_20000_DMROO`), "beta\\[2,")])

## plotting the histogram of the posteriors and the point MLE
replot = FALSE
if(replot){
  for(nme in names(results_TMB_DM)[grepl('signatures', names(results_TMB_DM))]){
    if(typeof(results_TMB_DM[[nme]]) == "character"){next}
    results_TMB_DM[[nme]]$par.fixed[grep("beta", names(results_TMB_DM[[nme]]$par.fixed))]
    
    idx_posteriors = which(grepl(nme, sapply(names(betasDM$posteriors_betas), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = ""))))
    if(! (length(idx_posteriors) == 1)){next}
    
    d_min_1 = sum(grepl('beta\\[2,', colnames(betasDM$posteriors_betas[[idx_posteriors]])))
    
    pdf(paste0("../../results/TMB/MLE_stan_comparison_DM_double_lambda/", nme, ".pdf"))
    par(mfrow=c(4,4))
    for(idx in 1:d_min_1){
      beta_val = select_slope_2(results_TMB_DM[[nme]]$par.fixed[grep("beta", names(results_TMB_DM[[nme]]$par.fixed))])[idx]
      plot(density(betasDM$posteriors_betas[[idx_posteriors]][,python_like_select(colnames(betasDM$posteriors_betas[[idx_posteriors]]), "beta\\[2,")][,idx]),
           main=paste0(idx, ' ', round(beta_val, digits = 2)))
      abline(v=beta_val,
             col='blue', lty="dashed")
    }
    dev.off()
  }
}
#------------------------------------------------------------------------------------------------------------------------#

## Basic checks.
## 1. Check that there are negative entries in covariance matrix
sapply(results_TMB_fullRE_M, function(i) try(python_like_select_name(i$par.fixed, 'cov_par_RE')))
## (there are)
## 2. Check matrix u for random effects
results_TMB_fullRE_M[[1]]$par.random
.dmin1 = length(python_like_select_name(results_TMB_fullRE_M[[1]]$par.fixed, "beta"))/2 ## this is d-1
matrix(results_TMB_fullRE_M[[1]]$par.random, ncol=.dmin1) ## this is the matrix of random effects
matrix(results_TMB_fullRE_M[[1]]$par.random, ncol=.dmin1) ## this is the matrix of random effects

#------------------------------------------------------------------------------------------------------------------------#
## match stan and TMB
subset_results_TMB_M = results_TMB_M[match(sapply(names(betasM$posteriors_betas_slope), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = "")),
                                           names(results_TMB_M))]
subset_results_TMB_DM = results_TMB_DM[match(sapply(names(betasDM$posteriors_betas_slope), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = "")),
      names(results_TMB_DM))]

subset_results_TMB_DM_dep = results_TMB_DM_dep[match(sapply(names(betasDM$posteriors_betas_slope), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = "")),
                                             names(results_TMB_DM_dep))]

## compare betas (note that some gave NaN for their standard error)
unlisted_TMB_list = function(TMB_list){
  sapply(1:length(TMB_list), function(dataset){
  if(typeof(TMB_list[[dataset]]) != "character"){
    TMB_list[[dataset]]$par.fixed[grepl('beta', names(TMB_list[[dataset]]$par.fixed))]
  }else{
    NA
  }
})
}

betas_TMB_M_subset = unlisted_TMB_list(subset_results_TMB_M)
unlisted_TMB_M = sapply(betas_TMB_M_subset, select_slope_2, verbatim=FALSE)
betas_TMB_DM_subset = unlisted_TMB_list(subset_results_TMB_DM)
unlisted_TMB_DM = sapply(betas_TMB_DM_subset, select_slope_2, verbatim=FALSE)
names(unlisted_TMB_M) = names(subset_results_TMB_M)
names(unlisted_TMB_DM) = names(subset_results_TMB_DM)

unlisted_stan_mean_M = sapply(betasM$posteriors_betas_slope, function(i) colMeans(i))
unlisted_stan_mean_DM = sapply(betasDM$posteriors_betas_slope, function(i) colMeans(i))

unlisted_RE_coef_stan_mean_M = 0 ## I would need to get those from the cluster

for(i in which(is.na(names(unlisted_TMB_DM)))){
  ## need to find the number of features d for these missing tmb runs
  unlisted_TMB_DM[[i]] = rep(NA, ncol(betasDM$posteriors_betas_slope[[names(betasDM$posteriors_betas_slope)[i]]]))
}

## make sure they are the same length
length(unlisted_TMB_DM)
length(unlisted_stan_mean_DM)

## make sure they are the same length
length(unlisted_TMB_DM %>% unlist)
length(unlisted_stan_mean_DM %>% unlist)
## bad

## TRUE, FALSE, colour scheme
myColors <- c('#5ad3b7', '#d35a6c')
names(myColors) <- c(T, F)
colScale <- scale_colour_manual(name = "grp",values = myColors)

par(mfrow=c(1,1))
## scatterplot for multinomial draws of the mean beta from the posterior and the MLE for beta
plot(unlisted_TMB_M %>% unlist, unlisted_stan_mean_M %>% unlist,
     col=c('blue', 'red')[factor(sapply(subset_results_TMB_M, give_summary_per_sample) == "Good", levels=c(T,F))], pch=8, cex=0.5)
df_betas_slope_M = cbind.data.frame( TMB_estimate=(unlisted_TMB_M %>% unlist),
                                     stan_mean=(unlisted_stan_mean_M %>% unlist),
                                     Bool_good_convergence=factor(rep(sapply(subset_results_TMB_M, give_summary_per_sample) == "Good",
                                                                      as.vector(sapply(subset_results_TMB_M,
                                                                                       function(i) sum(names(i$par.fixed) == "beta")/2)))), ## /2 because of intercept+slope
                                     ct = rep(sapply(names(unlisted_stan_mean_M), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = " ")),
                                              sapply(unlisted_TMB_M, length)))
ggplot(df_betas_slope_M,
       aes(x=TMB_estimate, y=stan_mean, col=Bool_good_convergence))+
  geom_point()+ggtitle('Scatterplot of mean of\nposterior and TMB ML estimate')+
  geom_abline(slope = 1, intercept = 0)+theme(legend.position = "bottom")+colScale
# ggsave("../../results/TMB/betas_stan_TMB_M.pdf")

## same plot, but trying to see where it went well and where it did not
ggplot(df_betas_slope_M,
       aes(x=TMB_estimate, y=stan_mean, col=Bool_good_convergence))+
  geom_point()+ggtitle('Scatterplot of mean of posterior and TMB ML estimate')+
  geom_abline(slope = 1, intercept = 0)+theme(legend.position = "bottom")+facet_wrap(.~factor(ct))+colScale
# ggsave("../../results/TMB/betas_stan_TMB_M_per_ct.pdf", width = 12)


## scatterplot for dirichlet-multinomial draws of the mean beta from the posterior and the MLE for beta
plot(unlisted_TMB_DM %>% unlist, unlisted_stan_mean_DM %>% unlist, col=c('blue', 'red')[factor(sapply(subset_results_TMB_DM, give_summary_per_sample) == "Good", levels=c(T,F))], pch=8)
df_betas_slope_DM = cbind.data.frame( TMB_estimate=(unlisted_TMB_DM %>% unlist),
                                      stan_mean=(unlisted_stan_mean_DM %>% unlist),
                                      #TMB_estimate_dep=(unlisted_TMB_DM_dep %>% unlist),
                                      Bool_good_convergence=factor(rep(sapply(subset_results_TMB_DM, give_summary_per_sample) == "Good",
                                                                       as.vector(sapply(unlisted_TMB_DM, length)))), ## /2 because of the intercept and the slope
                                      ct = rep(sapply(names(unlisted_stan_mean_DM), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = " ")),
                                               sapply(unlisted_TMB_DM, length)))
ggplot(df_betas_slope_DM,
       aes(x=TMB_estimate, y=stan_mean, col=factor(Bool_good_convergence, levels=c(F, T))))+
  geom_point()+ggtitle('Scatterplot of mean of\nposterior and TMB ML estimate')+theme(legend.position = "bottom")+
  geom_abline(intercept = 0, slope = 1)+colScale
# ggsave("../../results/TMB/betas_stan_TMB_DM.pdf")

ggplot(df_betas_slope_DM %>% filter(Bool_good_convergence == TRUE),
       aes(x=TMB_estimate, y=stan_mean, col=factor(Bool_good_convergence, levels=c(F, T))))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of\nposterior and TMB ML estimate')+
  theme(legend.position = "bottom")+colScale
# ggsave("../../results/TMB/betas_stan_TMB_DM2.pdf")

## split by samples that had good convergence and samples that didn't
ggplot(df_betas_slope_DM,
       aes(x=TMB_estimate, y=stan_mean))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of\nposterior and TMB ML estimate')+
  theme(legend.position = "bottom")+facet_wrap(.~factor(Bool_good_convergence, levels=c(F, T)))

## checking the scatterplot of the previous version which didn't have an intercept in beta
ggplot(df_betas_slope_DM,
       aes(x=TMB_estimate, y=stan_mean))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of posterior and TMB ML\n (beta estimate for DM)')+
  theme(legend.position = "bottom")+facet_wrap(.~ct, scales = "free")


## same plot, but only with good convergence
ggplot(df_betas_slope_DM %>% filter(Bool_good_convergence == TRUE),
       aes(x=TMB_estimate, y=stan_mean))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of posterior and TMB ML estimate')+
  theme(legend.position = "bottom")+facet_wrap(.~ct, scales = "free")
# ggsave("../../results/TMB/betas_stan_TMB_DM_per_ct.pdf", width = 8)

ggplot(df_betas_slope_DM,
       aes(x=TMB_estimate_dep, y=stan_mean))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of posterior and TMB ML estimate')+
  theme(legend.position = "bottom")+facet_wrap(.~ct, scales = "free")

## For which cancer types is there a good correlation in DM?
df_betas_slope_DM %>% filter(Bool_good_convergence == TRUE) %>% select(ct) %>% unique() %>% unlist %>% as.vector

#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## when are random effects negligible?
## note: this is for single intercept for RE
RE_TMB_M = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){i$par.random}else{NA})
RE_TMB_DM = sapply(results_TMB_DM, function(i) if(typeof(i) == "list"){i$par.random}else{NA})

RE_TMB_M_subset = sapply(subset_results_TMB_M, function(i) if(typeof(i) == "list"){i$par.random}else{NA})
RE_TMB_DM_subset = sapply(subset_results_TMB_DM, function(i) if(typeof(i) == "list"){i$par.random}else{NA})
names(RE_TMB_M_subset) = names(subset_results_TMB_M); names(RE_TMB_DM_subset) = names(subset_results_TMB_DM)

## distribution of the random effects for each cancer type, and for all DM runs
ggplot(reshape2::melt(RE_TMB_DM), aes(x=value))+geom_density()+facet_wrap(.~L1, scales = "free")
ggplot(reshape2::melt(RE_TMB_DM), aes(x=value))+geom_density()+facet_wrap(.~L1, scales = "free_y")
## distribution of the random effects, pooling all cancer types and patients
ggplot(reshape2::melt(RE_TMB_DM), aes(x=log(abs(value))))+geom_density()
ggplot(reshape2::melt(RE_TMB_DM), aes(x=log(abs(value))))+geom_histogram()
## large enough RE are positive (difficult to see as most RE are quite low)

## mean of the absolute random effect coefficients, for each cancer type
RE_TMB_M_mean = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){mean(abs(i$par.random))}else{NA})
RE_TMB_DM_mean = sapply(results_TMB_DM, function(i) if(typeof(i) == "list"){mean(abs(i$par.random))}else{NA})
plot(density(na.omit(RE_TMB_DM_mean)))

## how to classify cancer types between negligible and non-negligible RE?
sort(RE_TMB_DM_mean)

## is there a correlation between RE in signatures and nucleotide substitutions?
plot(vector_to_ct_list(RE_TMB_M_mean))
plot(vector_to_ct_list(RE_TMB_M_mean), ylim = c(0,1))
plot(vector_to_ct_list(RE_TMB_DM_mean))
plot(vector_to_ct_list(RE_TMB_DM_mean), ylim = c(0,0.005))
## looks like an appalling correlation

## Maybe I should instead look at the ratio of random effects to fixed effects?
RE_TMB_M_median_ratio = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){
  median(abs(i$par.random))/median(abs(select_slope_2(python_like_select_name(i$par.fixed, 'beta'), verbatim = FALSE)))}else{NA})
plot(vector_to_ct_list(RE_TMB_M_median_ratio), ylim = c(0,10)) ## this is a marginally better correlation

## compare directly the random effects variance
RE_TMB_M_REvar = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){python_like_select_name(i$par.fixed, 'logSigma_RE')}else{NA})
## normalised by size of fixed effects (here I am normalising by ALL beta; slope and intercept)
RE_TMB_M_REvar_norm = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){python_like_select_name(i$par.fixed, 'logSigma_RE')/median(abs(python_like_select_name(i$par.fixed, 'beta')))}else{NA})
RE_TMB_M_REvar_norm_slope = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){python_like_select_name(i$par.fixed, 'logSigma_RE')/median(abs(select_slope_2(python_like_select_name(i$par.fixed, 'beta'), verbatim = FALSE)))}else{NA})
plot(vector_to_ct_list(RE_TMB_M_REvar)) ## this is the log
plot(vector_to_ct_list(RE_TMB_M_REvar), xlim = c(-5,4), ylim = c(-4, 3)) ## this is the log
## there doesn't seem to be a good agreement either
plot(vector_to_ct_list(RE_TMB_M_REvar_norm))
plot(vector_to_ct_list(RE_TMB_M_REvar_norm_slope))

#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## compute r*2 for betas, random effects

## we need to subset the random effects too, to match it with stan
RE_TMB_M_subset
betas_TMB_DM_subset
unlisted_stan_mean_M


#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Determine differential abundance
results_TMB_M

# pvals_M = sapply(results_TMB_M, wald_TMB_wrapper)
# pvals_DM = sapply(results_TMB_DM, wald_TMB_wrapper)
# pvals_LNM = sapply(results_TMB_LNM, wald_TMB_wrapper)
pvals_M_fullRE = sapply(results_TMB_fullRE_M, wald_TMB_wrapper)
pvals_M_diagRE = sapply(results_TMB_diagRE_M, wald_TMB_wrapper)
pvals_M_fullRE_good = pvals_M_fullRE[sapply(results_TMB_fullRE_M, give_summary_per_sample) == "Good"]
pvals_DM_fullRE = sapply(results_TMB_fullRE_DM, function(i) try(wald_TMB_wrapper(i)))
pvals_DM_fullRE_good = pvals_DM_fullRE[sapply(results_TMB_fullRE_DM, give_summary_per_sample) == "Good"]
pvals_DM_diagRE = sapply(results_TMB_diagRE_DM, function(i) try(wald_TMB_wrapper(i)))
pvals_fullRE_DMSL = sapply(results_TMB_fullRE_DMSL, function(i) try(wald_TMB_wrapper(i)))
pvals_diagRE_DMSL = sapply(results_TMB_diagRE_DMSL, function(i) try(wald_TMB_wrapper(i)))
pvals_DM_fullRE_nlminb = sapply(results_TMB_fullRE_DM_nlminb, function(i) try(wald_TMB_wrapper(i)))
pvals_DMSL_fullRE_nlminb = sapply(results_TMB_fullRE_DMSL_nlminb, function(i) try(wald_TMB_wrapper(i)))
pvals_DMSL_diagRE_nlminb = sapply(results_TMB_diagRE_DMSL_nlminb, function(i) try(wald_TMB_wrapper(i)))

sapply(list(M_single=pvals_M, DM_single=pvals_DM, LNM_single=pvals_LNM, M_full=pvals_M_fullRE, DM_full=pvals_DM_fullRE), max, na.rm = TRUE)

cat(paste0(apply(cbind(names(pvals_M_fullRE), pvals_M_fullRE < 0.05), 1, paste0, collapse='\t'), collapse = '\n'))
cat(paste0(sapply(c(pvals_M_fullRE < 0.05)[match(names(pvals_M_fullRE), names(pvals_DM_fullRE_good))],
                 paste0, collapse='\t'), collapse = '\n'))
cat(paste0(sapply(c(pvals_M_diagRE < 0.05)[match(names(pvals_M_fullRE), names(pvals_M_diagRE))],
                  paste0, collapse='\t'), collapse = '\n'))
cat(paste0(sapply(c(pvals_DM_diagRE < 0.05)[match(names(pvals_M_fullRE), names(pvals_DM_diagRE))],
                  paste0, collapse='\t'), collapse = '\n'))
cat(paste0(sapply(c(pvals_fullRE_DMSL < 0.05)[match(names(pvals_M_fullRE), names(pvals_fullRE_DMSL))],
                  paste0, collapse='\t'), collapse = '\n'))
cat(paste0(sapply(c(pvals_diagRE_DMSL < 0.05)[match(names(pvals_M_fullRE), names(pvals_diagRE_DMSL))],
                  paste0, collapse='\t'), collapse = '\n'))
cat(paste0(sapply(c(pvals_DM_fullRE_nlminb < 0.05)[match(names(pvals_M_fullRE), names(pvals_DM_fullRE_nlminb))],
                  paste0, collapse='\t'), collapse = '\n'))
cat(paste0(sapply(c(pvals_DMSL_fullRE_nlminb < 0.05)[match(names(pvals_M_fullRE), names(pvals_DMSL_fullRE_nlminb))],
                  paste0, collapse='\t'), collapse = '\n'))
cat(paste0(sapply(c(pvals_DMSL_diagRE_nlminb < 0.05)[match(names(pvals_M_fullRE), names(pvals_DMSL_diagRE_nlminb))],
                  paste0, collapse='\t'), collapse = '\n'))

par(mfrow=c(1,2))
plot(pvals_M, pvals_M_fullRE)
plot(pvals_DM, pvals_DM_fullRE[match(names(results_TMB_DM), names(results_TMB_fullRE_DM))])


names(pvals_DM)[is.na(match(names(pvals_DM), names(pvals_DM_fullRE)))]

all_names = unique(c(names(results_TMB_fullRE_M), names(results_TMB_fullRE_DM)))
betasM = sapply(all_names, function(i) try(python_like_select_name(results_TMB_fullRE_M[[i]]$par.fixed, 'beta')))
betasDM = sapply(all_names, function(i) try(python_like_select_name(results_TMB_fullRE_DM[[i]]$par.fixed, 'beta')))
betasDM[sapply(betasM, length) == 0] = NA
betasM[sapply(betasDM, length) == 0] = NA

length(betasM)
length(betasDM)

length(unlist(betasM))
length(unlist(betasDM))

## Differential abundance for full RE
table(pvals_M_fullRE_good <= 0.05)
table(pvals_DM_fullRE_good <= 0.05)

## Plot the betas
df_all_beta_slopes = melt(list(M=sapply(betasM[sapply(betasM, typeof) == "double"], select_slope_2, v=F),
                               DM=sapply(betasDM[sapply(betasM, typeof) == "double"], select_slope_2, v=F)))
df_all_beta_slopes$ct = gsub("signatures|nucleotidesubstitution1", "", df_all_beta_slopes$L2)
df_all_beta_slopes$type = sapply(1:nrow(df_all_beta_slopes), function(i) gsub(df_all_beta_slopes$ct[i], "", df_all_beta_slopes$L2[i]) )
df_all_beta_slopes$converged = 3

df_all_beta_slopes[df_all_beta_slopes$L1 == 'M','converged'] = unlist((sapply(results_TMB_fullRE_M, `[`, "pdHess"))[match(df_all_beta_slopes[df_all_beta_slopes$L1 == 'M',]$L2, names(results_TMB_fullRE_M))])
df_all_beta_slopes[df_all_beta_slopes$L1 == 'DM','converged'] = unlist((sapply(results_TMB_fullRE_DM, `[`, "pdHess"))[match(df_all_beta_slopes[df_all_beta_slopes$L1 == 'DM',]$L2, names(results_TMB_fullRE_DM))])

ggplot(df_all_beta_slopes,
       aes(x=ct, y=value, col=type, shape=converged))+geom_point()+facet_wrap(.~interaction(type, L1), nrow=2, scales = "free_y")+
  geom_abline(intercept = 0, slope = 0, lty='dashed')+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.5, hjust=1))+theme(legend.position = "bottom")+labs(x='')
ggsave("../../results/results_TMB/pcawg/betas_TMB_summary.pdf", width = 10, height = 7)

#----------------------------------------------------------------------------------------------------------

table(pvals_M <= 0.05, pvals_DM <= 0.05)

sapply(list(pvals_M, pvals_DM, pvals_LNM), function(i) sum(na.omit(i<= 0.05))/sum(!is.na(i)))

## why??!
names(pvals_M)[!(names(pvals_M) %in% names(pvals_LNM))]


#----------------------------------------------------------------------------------------------------------

list_results = list(results_TMB_M, results_TMB_DM, results_TMB_LNM)
cts = sapply(list_results, function(i){
  gsub("nucleotidesubstitution1|signatures", "", names(i))
})

df_idx = lapply(1:3, function(idx){
  t(sapply(unique(cts[[idx]]),
           function(j) c(which(cts[[idx]] == j & grepl("nucleotidesubstitution1", names(list_results[[idx]]))),
                         which(cts[[idx]] == j & grepl("signature", names(list_results[[idx]]))))))
})

clnm <- colnames(df_idx[[3]])
df_idx[[3]] = do.call('rbind', df_idx[[3]])
rownames(df_idx[[3]]) = clnm

#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## full RE
## we have the covariance matrix of the random effects (cov_par_RE) and the scaling factors (logs_sd_RE)

plot_random_effects = function(idx_ct){
  if(typeof(results_TMB_fullRE_M[[idx_ct]]) == "character"){
    plot(0, main='Run failed')
  }else{
    cov_RE = give_UNSTRUCTURED_CORR_t_matrix(vec = python_like_select_name(results_TMB_fullRE_M[[idx_ct]]$par.fixed, 'cov_par_RE'),
                                    dim_mat = length(python_like_select_name(results_TMB_fullRE_M[[idx_ct]]$par.fixed, 'logs_sd_RE')))
    ## scsaling (check if this is really how it's done!)
    cov_RE = cov_RE %*% diag(python_like_select_name(results_TMB_fullRE_M[[idx_ct]]$par.fixed, 'logs_sd_RE'))
    ## see how these variances compare to the 
    
    ## d-1
    # length(python_like_select_name(results_TMB_fullRE_M[[1]]$par.fixed, 'beta'))/2
    # cov_RE[,1]
    
    if(length(cov_RE[,1]) > 15){
      par(mfrow=c(1,1))
      plot(0, main='Too many features')
    }else{
      if(length(cov_RE[,1]) > 6){
        par(mfrow=c(1+ (length(cov_RE[,1]) %/% 6), 6))
      }else{
        par(mfrow=c(1,length(cov_RE[,1])))
      }
      sapply(1:length(cov_RE[,1]), function(j){
        plot(density(rnorm(6000, mean = 0, sd = exp(cov_RE[j,j]))), lty='dashed', col='blue')
        lines(density(matrix(results_TMB_fullRE_M[[idx_ct]]$par.random,
                             nrow=length(python_like_select_name(results_TMB_fullRE_M[[idx_ct]]$par.fixed,
                                                                 'logs_sd_RE')))[j,]), col='blue')
        
        ## and print the rest as background
        for(j2 in (1:length((cov_RE)[,1]))[-j]){
          lines(density(rnorm(6000, mean = 0, sd = exp(cov_RE[j2,j2]))), lty='dotted', col=alpha('black', 0.8))
        }
      })
    }
  }
}

pdf("../../results/assessing_models/RE_PCAWG_TMB_fullRE_M.pdf", width = 30, height = 5)
sapply(1:length(results_TMB_fullRE_M), plot_random_effects)
# sapply(67, plot_random_effects)
dev.off()
#----------------------------------------------------------------------------------------------------------

## Load objects
datasets = list.files("../../data/roo/", full.names = T)
datasets = datasets[!grepl("nucleotidesubstitution3", datasets)]
datasets_loaded = sapply(datasets, readRDS)

results_TMB_fullRE_DM_matched = results_TMB_fullRE_DM[match(gsub("_", "", gsub("_ROO.RDS", "", basename(datasets))),
                                                            names(results_TMB_fullRE_DM))]
names(results_TMB_fullRE_DM_matched) = datasets
results_TMB_fullRE_M_matched = results_TMB_fullRE_M[match(gsub("_", "", gsub("_ROO.RDS", "", basename(datasets))),
                                                            names(results_TMB_fullRE_M))]
names(results_TMB_fullRE_M_matched) = datasets

give_cov_matrix_annotated = function(.dataset_object, dataset_idx, verbatim=T){
  if(!is.null(.dataset[[dataset_idx]]) & !(typeof(.dataset[[dataset_idx]]) == "character")){
    if(verbatim) stop('Note: by default, the logratio is taken using the last category as baseline, but that is not necessarily, because it can be changed in the wrapper TMB')
    .obj = datasets_loaded[[dataset_idx]]
    if(!(length(attr(.obj,"count_matrices_active")[[1]]) == 0)){
      .colnames = colnames(attr(.obj,"count_matrices_active")[[1]])
    }else{
      .colnames = colnames(attr(.obj,"count_matrices_all")[[1]])
    }
    .cov_mat= fill_covariance_matrix(arg_d = length(python_like_select_name(.dataset[[dataset_idx]]$par.fixed, 'beta'))/2,
                                     arg_entries_var = (python_like_select_name(.dataset[[dataset_idx]]$par.fixed, 'logs_sd_RE'))**2,
                                     arg_entries_cov = python_like_select_name(.dataset[[dataset_idx]]$par.fixed, 'cov_par_RE'))
    colnames(.cov_mat) = rownames(.cov_mat) = .colnames[-length(.colnames)]
    .cov_mat
  }else{
    NA
  }
}

## Scatterplot of covariances
covariance_matrices = lapply(list(results_TMB_fullRE_M_matched,
                                  results_TMB_fullRE_DM_matched), function(.dataset){
   .x = sapply(1:length(.dataset), give_cov_matrix_annotated, .dataset_object = .dataset, verbatim=FALSE)
   names(.x) = gsub("_ROO.RDS", "", basename(names(.dataset)))
   return(.x)
})
names(covariance_matrices) = c('fullRE_M', 'fullRE_DM')

covariance_matrices_melt = as.data.frame(melt(covariance_matrices))
covariance_matrices_melt = covariance_matrices_melt[!is.na(covariance_matrices_melt$value),]
covariance_matrices_melt = rbind(covariance_matrices_melt, cbind(Var1=covariance_matrices_melt[,2], Var2=covariance_matrices_melt[,1],
                                                                 covariance_matrices_melt[,3:5]))
covariance_matrices_melt_substit = covariance_matrices_melt[grepl('nucleotidesubstitution1', covariance_matrices_melt$L2),]
covariance_matrices_melt_sigs = covariance_matrices_melt[grepl('signatures', covariance_matrices_melt$L2),]

## removing extreme values
ggplot(covariance_matrices_melt_sigs, aes(x=value))+geom_density()
extreme_vals_cov_sigs = quantile(covariance_matrices_melt_sigs$value, probs = c(0.01, 0.99))
covariance_matrices_melt_sigs[ !((covariance_matrices_melt_sigs$value > extreme_vals_cov_sigs[1]) & (covariance_matrices_melt_sigs$value < extreme_vals_cov_sigs[2])),'value'] = NA
covariance_matrices_melt_sigs = covariance_matrices_melt_sigs[!is.na(covariance_matrices_melt_sigs$ value),]
head(covariance_matrices_melt_sigs)
ggplot(covariance_matrices_melt_sigs, aes(x=Var1, y=Var2, col=value))+geom_jitter(size=.2)+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 1)+labs(col = "Covariance")+
  theme_bw()+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../../results/results_TMB/pcawg/covariances_signatures.pdf", width = 10, height = 10)

subset_sigs = names(tail(sort(table(c(as.character(covariance_matrices_melt_sigs$Var1),
                                      as.character(covariance_matrices_melt_sigs$Var2)))), n=20))
ggplot(covariance_matrices_melt_sigs[!apply(covariance_matrices_melt_sigs[,1:2], 1, function(i){
  any(!(i %in% subset_sigs))}),],
  aes(x=Var1, y=Var2, col=value))+geom_jitter(size=.2)+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 1)+labs(col = "Covariance")+
  theme_bw()+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../../results/results_TMB/pcawg/covariances_signatures_subset.pdf", width = 7, height = 6.7)

ggplot(covariance_matrices_melt_substit, aes(x=Var1, y=Var2, col=value))+geom_jitter()+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 4)+labs(col = "Covariance")
  

