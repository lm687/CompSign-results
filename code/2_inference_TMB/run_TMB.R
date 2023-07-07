## Creating dataset to assess the performance of models for parameter recovery

rm(list = ls())

library(TMB)
library(optparse)

source("1_create_ROO/roo_functions.R")
source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("2_inference_TMB/mm_multinomial/helper_functions.R")
source("2_inference_TMB/helper_TMB.R")

debug <- F
if(debug){
  opt <- list()
  opt$input = '../data/roo/Liver-HCC_signaturesMSE_ROO.RDS'
  opt$output = '../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_Liver-HCC_signaturesMSE.RDS'
  opt$model = 'diagREDM' 
  opt$feature_type = 'signaturesMSE' 
  opt$optimiser = 'nlminb' 
  opt$simulation_bool = F 
  opt$read_directly = T 
  opt$use_previous_run_startingvals  = T
  
}


option_list = list(
  make_option(c("--model"), type="character", default=NA,
              help="Which model to use for inference", metavar="character"),
  make_option(c("--input"), type="character", default=NA,
              help="Input file with dataset (RDS)", metavar="character"),
  make_option(c("--simulation_bool"), type="logical", default=T,
              help="Is ct the name of the file to read?", metavar="logical"),
  make_option(c("--read_directly"), type="logical", default=T,
              help="Should the opt$input file be read directly with load_PCAWG?", metavar="logical"),
  make_option(c("--feature_type"), type="character", default="signatures",
              help="Type of feature: signatures, signaturesPCAWG, etc", metavar="character"),
  make_option(c("--output"), type="character", default=NA,
              help="Output file in which to write the results of the inference (RDS file)", metavar="character"),
  make_option(c("--optimiser"), type="character", default="optim",
              help="Which optimiser to use", metavar="character"),
  make_option(c("--nonexo_bool"), type="logical", default=F,
              help="Should only nonexogenous signatures be selected?", metavar="logical"),
  make_option(c("--use_previous_run_startingvals"), type="logical", default=F,
              help="Should we use any previous run, if available, to use the estimated values as starting values? The previous run should be called the same as the output, but with _NC.RDS instead of .RDS", metavar="logical"),
  make_option(c("--return_opt_bool"), type="logical", default=F,
              help="Should the RDS file contain the opt (TRUE for LRT, FALSE OTHERWISE)", metavar="logical")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## by default
simulation_bool = opt$simulation_bool

if(opt$optimiser == "optim"){
  use_nlminb = F
}else{
  use_nlminb = T
}

cat('Model:', opt$model, '\n')
cat('Feature type:', opt$feature_type, '\n')
cat('Using nlminb:', use_nlminb, '\n')
cat('Simulation boolean:', simulation_bool, '\n')

if(opt$nonexo_bool | grepl('nonexo', opt$output)){
  opt$model <- gsub("wSBS1SBS5nonexo", "", opt$model)
  opt$model <- gsub("nonexo", "", opt$model)
}
cat('Model:', opt$model, '\n')

if(opt$model == "fullREM"){
  TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_multinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_multinomial"))
  mod_model_name = "fullRE_M"
}else if(opt$model == "diagREM"){
  TMB::compile("2_inference_TMB/mm_multinomial/diagRE_ME_multinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_ME_multinomial"))
  mod_model_name = "diagRE_M"
}else if(opt$model == "fullREDM"){
  TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial"))
  mod_model_name = "fullRE_DM"
}else if(opt$model == "diagREDM"){
  TMB::compile("2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial"))
  mod_model_name = "diagRE_DM"
}else if(opt$model =="fullREDMsinglelambda"){
  TMB::compile("2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda"))
  mod_model_name = "fullREDMsinglelambda"
  # use_nlminb=T
}else if(opt$model =="diagREDMsinglelambda"){
  TMB::compile("2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda"))
  mod_model_name = "diagREDMsinglelambda"
  # use_nlminb=T
}else if(opt$model =="FEDMsinglelambda"){
  TMB::compile("2_inference_TMB/mm_multinomial/FE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/FE_dirichletmultinomial_single_lambda"))
  mod_model_name = "FEDMsinglelambda"
}else if(opt$model =="fullREDMnoscaling"){
  TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomialnoscaling.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomialnoscaling"))
  mod_model_name = "fullREDMnoscaling"
}else if(opt$model =="fullREDMonefixedlambda"){
  TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial_onefixedlambda.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial_onefixedlambda"))
  mod_model_name = "fullRE_DMonefixedlambda"
}else if(opt$model =="fullREDMonefixedlambda2"){
  TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial_onefixedlambda2.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial_onefixedlambda2"))
  mod_model_name = "fullRE_DMonefixedlambda2"
}else if(opt$model =="fullREDMonefixedlambda3"){
  TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial_onefixedlambda3.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial_onefixedlambda3"))
  mod_model_name = "fullRE_DMonefixedlambda3"
}else if(opt$model == "singleREDM"){
  TMB::compile("2_inference_TMB/mm_multinomial/singleRE_dirichlet_multinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/singleRE_dirichlet_multinomial"))
  mod_model_name = "singleRE_DM"
}else if(opt$model == "diagREDMreducedmodel"){
  TMB::compile("2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial"))
  mod_model_name = "diagRE_DMreducedmodel"
}else if(opt$model == "diagREcorDM"){
  TMB::compile("2_inference_TMB/mm_multinomial/diagREcor_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagREcor_ME_dirichletmultinomial"))
  mod_model_name = "diagREcor_DM"
}else{
  stop('Specifiy a valid <model>')
}


cat(opt$input)
dataset = load_PCAWG(ct = opt$input, typedata = opt$feature_type, simulation = simulation_bool,
                     path_to_data = NA, read_directly=opt$read_directly, override_warning_X_Z=T)

if(opt$nonexo_bool | grepl('nonexo', opt$output)){
  if(grepl('wSBS1SBS5nonexo', opt$output)){
    ## including SBS1, SBS5
    nonexogenous = read.table("../data/cosmic/exogenous_signatures_SBS_withSBS1SBS5.txt", sep = "\t",
                              comment.char = "#", fill = F)
  }else{
    ## select only nonexogenous signatures
    ## not including SBS1, SBS5
    nonexogenous = read.table("../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t",
                              comment.char = "#", fill = F)
  }
  dataset <- give_subset_sigs_TMBobj(dataset, sigs_to_remove = nonexogenous$V1)
}

## use previous values as starting values
if(opt$use_previous_run_startingvals){
  cat('Using initial values from previous runs\n')
  outfile_not_converged <- gsub(".RDS", "_NC.RDS", opt$output)
  
  ## read a file that says how many times we have tried to make it converge, already. If it exceeds a threshold, save
  ## the output file with the original name, instead of the _NC.RDS name.
  ## get the number of previous tries, and add one to the <num_tries_for_convergence> document
  
  threshold_num_tries = 3 ###10
  header_num_tries='## Write down the number of tries of getting results that converge. If the number exceeds some limit, save the non-converged result\n'

  #file_num_tries <- "logs/num_tries_for_convergence.txt"
  file_num_tries <- paste0('logs/num_tries_for_convergence/', basename(opt$input))

  if(!file.exists(file_num_tries)){
    #system(paste0('touch ', file_num_tries))
    write(paste0(header_num_tries, 'dummy\tdummy\n' ), file_num_tries)
  }
  cat('Reading file', file_num_tries, '\n')
  num_tries_for_convergence <- read.table(file_num_tries, sep = '\t', comment.char = '#', stringsAsFactors=F, fill=T, row.names=NULL)
  
  ## check if there is an entry for the dataset under consideration
  which_num_tries <- which(num_tries_for_convergence$V1 == opt$output)
  which_num_tries = which_num_tries[length(which_num_tries)]
  cat('The current number of tries for this run is ', num_tries_for_convergence[which_num_tries, 2], '\n')
  print(num_tries_for_convergence[which_num_tries,2])
  if(length(which_num_tries)>0){
    num_tries_for_convergence[which_num_tries,2] = as.numeric(num_tries_for_convergence[which_num_tries,2]) + 1
    num_tries_for_convergence_append = t(matrix(c(num_tries_for_convergence[which_num_tries,1], as.character(as.numeric(num_tries_for_convergence[which_num_tries,2]) + 1))))
  }else{
    ## if not, append it
    cat('Appending\n')
    cat(opt$output, '\n')
    num_tries_for_convergence[which_num_tries,2] = 1 
    num_tries_for_convergence <- rbind(num_tries_for_convergence, c(opt$output, 1))
    num_tries_for_convergence_append = num_tries_for_convergence
  }

  
  ## save updated file
  cat(header_num_tries, file = file_num_tries)
  write.table(num_tries_for_convergence_append, quote = F, sep = "\t", file = file_num_tries, append = T, col.names = FALSE, row.names = FALSE)

  
  good_previous_nonconvergedf_file1 <- file.exists(outfile_not_converged)
  if(good_previous_nonconvergedf_file1){
    good_previous_nonconvergedf_file2 = typeof(readRDS(outfile_not_converged)) != "character"
  }else{
    good_previous_nonconvergedf_file2 = F
  }
  cat(' Previous file exists? (bool): ', good_previous_nonconvergedf_file1, '\n',
      'Was it a non-error run? (bool)', good_previous_nonconvergedf_file2, '\n')
  good_previous_nonconvergedf_file = good_previous_nonconvergedf_file1 & good_previous_nonconvergedf_file2
  ## use previous values as initial values
  if(good_previous_nonconvergedf_file){
    ## run with given initial parameters
    results_inference_previous <- readRDS(outfile_not_converged)
    
    ## depending on the model, the list of initial params is different
    ## need to compute dmin1
    ## need to compute dmin1
    if(opt$return_opt_bool){
      results_inference_previous <- results_inference_previous$return_report
    }
    
    if (!opt$model  %in% "diagREDMreducedmodel"){
      dmin1 <- length(python_like_select_name(results_inference_previous$par.fixed, 'beta'))/2
    }else{
      dmin1 <- length(python_like_select_name(results_inference_previous$par.fixed, 'beta'))
    }
    
    if(opt$model == "fullREM"){
      list_initial_params <- list(
        beta = matrix(python_like_select_name(results_inference_previous$par.fixed, 'beta'), nrow=2),
        u_large = matrix(results_inference_previous$par.random, ncol=dmin1),
        logs_sd_RE=python_like_select_name(results_inference_previous$par.fixed, 'logs_sd_RE'),
        cov_par_RE = python_like_select_name(results_inference_previous$par.fixed, 'cov_par_RE'))
    }else if(opt$model == "diagREM"){
      list_initial_params <- list(
        beta = matrix(python_like_select_name(results_inference_previous$par.fixed, 'beta'), nrow=2),
        u_large = matrix(results_inference_previous$par.random, ncol=dmin1),
        logs_sd_RE=python_like_select_name(results_inference_previous$par.fixed, 'logs_sd_RE'))
    }else if(opt$model == "fullREDM"){
      list_initial_params <- list(
        beta = matrix(python_like_select_name(results_inference_previous$par.fixed, 'beta'), nrow=2),
        u_large = matrix(results_inference_previous$par.random, ncol=dmin1),
        logs_sd_RE=python_like_select_name(results_inference_previous$par.fixed, 'logs_sd_RE'),
        cov_par_RE = python_like_select_name(results_inference_previous$par.fixed, 'cov_par_RE'),
        log_lambda = matrix(c(2,2)))
      }else if(opt$model == "fullREDMnoscaling"){
        list_initial_params <- list(
          beta = matrix(python_like_select_name(results_inference_previous$par.fixed, 'beta'), nrow=2),
          u_large = matrix(results_inference_previous$par.random, ncol=dmin1),
          cov_par_RE = python_like_select_name(results_inference_previous$par.fixed, 'cov_par_RE'),
          log_lambda = matrix(c(2,2)))
      }else if(opt$model %in% c("diagREDM", "diagREDMreducedmodel")){
          if(opt$model  == "diagREDM"){
            beta = matrix(python_like_select_name(results_inference_previous$par.fixed, 'beta'), nrow=2)
          }else if (opt$model  == "diagREDMreducedmodel"){
            beta = matrix(python_like_select_name(results_inference_previous$par.fixed, 'beta'), nrow=1)
          }
        list_initial_params <- list(
          beta=beta,
          u_large = matrix(results_inference_previous$par.random, ncol=dmin1),
          logs_sd_RE=python_like_select_name(results_inference_previous$par.fixed, 'logs_sd_RE'),
          log_lambda = matrix(c(2,2)))
      }else if(opt$model =="fullREDMsinglelambda"){
        list_initial_params <- list(
          beta = matrix(python_like_select_name(results_inference_previous$par.fixed, 'beta'), nrow=2),
          u_large = matrix(results_inference_previous$par.random, ncol=dmin1),
          logs_sd_RE=python_like_select_name(results_inference_previous$par.fixed, 'logs_sd_RE'),
          cov_par_RE = python_like_select_name(results_inference_previous$par.fixed, 'cov_par_RE'),
          log_lambda = 2,2)
        cat('Check log_lambda = 2,2. It should not affect results; it has simply added a "2" in the list')
    }else if(opt$model %in% c("fullREDMonefixedlambda", "fullREDMonefixedlambda2", "fullREDMonefixedlambda3")){
      list_initial_params <- list(
        beta = matrix(python_like_select_name(results_inference_previous$par.fixed, 'beta'), nrow=2),
        u_large = matrix(results_inference_previous$par.random, ncol=dmin1),
        logs_sd_RE=python_like_select_name(results_inference_previous$par.fixed, 'logs_sd_RE'),
        cov_par_RE = python_like_select_name(results_inference_previous$par.fixed, 'cov_par_RE'),
        log_lambda = 2)
    }else if(opt$model =="diagREDMsinglelambda"){
      list_initial_params <- list(
        beta = matrix(python_like_select_name(results_inference_previous$par.fixed, 'beta'), nrow=2),
        u_large = matrix(results_inference_previous$par.random, ncol=dmin1),
        logs_sd_RE=python_like_select_name(results_inference_previous$par.fixed, 'logs_sd_RE'),
        log_lambda = 2,2)
      cat('Check log_lambda = 2,2. It should not affect results; it has simply added a "2" in the list')
    }else if(opt$model =="FEDMsinglelambda"){
        stop('Custom initial values for FEDMsinglelambda: Not implemented')
    }else{
    }
    
    results_inference = try(wrapper_run_TMB(object = dataset, model = mod_model_name, use_nlminb=use_nlminb,
                                            initial_params = list_initial_params, return_opt=opt$return_opt_bool))
  }else{
    ## if there is no previous file, run without specified starting values
    results_inference = try(wrapper_run_TMB(object = dataset, model = mod_model_name, use_nlminb=use_nlminb,
                                            return_opt=opt$return_opt_bool))
  }
  
  if( typeof(results_inference) == "character"){
    not_converged = T
    error_run = T
  }else{
    error_run = F
    if(opt$return_opt_bool){
      if(!results_inference$return_report$pdHess){
        not_converged = T
      }else{
        not_converged = F
      }
    }else{
      if(!results_inference$pdHess){
        not_converged = T
      }else{
        not_converged = F
      }
    }
  }
  if(not_converged){
    ## if it still doesn't converge, save it as _NC
    cat('The run has not converged yet\n')
    print(which_num_tries)
    print(num_tries_for_convergence[which_num_tries,])
    cat('The current number of tries for this run is ', num_tries_for_convergence[which_num_tries, 2], '\n')
    
    if(num_tries_for_convergence[which_num_tries,2] > threshold_num_tries){
      if(typeof(results_inference) == "character"){
        ## if it was an <error> run: do not save unless it's the first run or the last one
        if((num_tries_for_convergence[which_num_tries,2] == 0) | (num_tries_for_convergence[which_num_tries,2] >= threshold_num_tries)){
          saveRDS(object = results_inference, file = opt$output)
        }
      }else{
        ## if it was a good run without convergence
        saveRDS(object = results_inference, file = opt$output)
      }
    }else{
      saveRDS(object = results_inference, file = outfile_not_converged)
    }
  }else{
    ## if it has converged, save it with the proper name
    saveRDS(object = results_inference, file = opt$output)
  }
}else{
  # dataset = sort_columns_TMB(dataset)
  # results_inference = try(wrapper_run_TMB(opt$input, model = mod_model_name, typedata = "simulation", simulation = TRUE))
  results_inference = try(wrapper_run_TMB(object = dataset, model = mod_model_name, use_nlminb=use_nlminb))
  
  saveRDS(object = results_inference, file = opt$output)
}


