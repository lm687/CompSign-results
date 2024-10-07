## replot the output of analyse_inference_simulations_integrate.R, if more plots are needed
## this way the p-values for competing methods do not need to be re-computed

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(ggrepel)
library(dplyr)
library(latex2exp)
library(reshape2)

source("../../3_analysis/2_simulation_model_assessment/2_analyse_inference_simulations/helper_model_assessment.R")
source("../../2_inference_TMB/helper_TMB.R")
source("../../1_create_ROO/roo_functions.R")

# flder_out <- "../../../results/results_TMB/simulated_datasets/mixed_effects_models_multiple_replot/"
# flder_in <- "../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries_multiple/"


find_string_files <- function(j) paste0("../../../data/assessing_models_simulation/inference_results/TMB/", j$optimiser, "/",
       'multiple_', j$dataset_generation, '_', j$n, '_',
       j$nlambda, '_', j$lambda, '_', j$d, '_', j$beta_gamma_shape, '_', j$model, '_',
       j$betaintercept, '_',
       j$betaslope, '_', j$cov, '/')

find_input_from_list <- function(j, keep_NC=F){
  require(TMB)
  
  string_files <- find_string_files(j)
  cat('Looking for files that match <', string_files, '>...\n')
  lst <- list.files(string_files, full.names = T)
  lst <- lst[grepl(paste0('multiple_', j$dataset_generation, '_'), lst)]
  if(!keep_NC){
    lst <- lst[!grepl('_NC.RDS', lst)]
  }
  lst <- lst[grepl(j$model, lst)]
  lst <- lst[grep(paste0(j$betaintercept, "_", j$betaslope, "_", j$cov, "*"), lst)]
  lst <- lst[grepl(paste0(j$n, '_', j$nlambda, '_', j$lambda, '_', j$d, '_', j$beta_gamma_shape), lst)]
  lst
}

opt_from_bias_fig <- function(i){
  isplit <- strsplit(i, '_')[[1]]
  opt <- list(model = isplit[10],
              dataset_generation = isplit[3],
              multiple_runs = T,
              run_nonconverged = F,
              optimiser = isplit[4],
              n = isplit[5],
              nlambda = isplit[6],
              lambda = isplit[7],
              d = isplit[8],
              beta_gamma_shape = isplit[9],
              betaintercept =  isplit[11],
              betaslope =  isplit[12],
              cov =  isplit[13]
  )
  opt$input = find_input_from_list(opt)
  opt$input_all = find_input_from_list(opt, keep_NC=T)
  opt$string_files = find_string_files(opt)
  
  return(opt)
}

opt_all_datasets <- list()
opt_all_datasets[['A1diag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_2_5_0_diagREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A1full']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_2_5_0_fullREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A1single']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_2_5_0_singleREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A2diag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_20_5_0_diagREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A2full']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_20_5_0_fullREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A2single']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_20_5_0_singleREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A3diag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_100_5_0_diagREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A3full']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_100_5_0_fullREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A3single']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_100_5_0_singleREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['baselinediag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_14072_80_6_0_diagREDM_betainterceptPCAWG4_betaslopeonechangingPCAWG4_covmatPCAWG4_onlyconverged_bias_betas_v2')
opt_all_datasets[['B1full']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_3401_18_6_0_fullREDM_betainterceptPCAWG2_betaslopePCAWG2_covmatPCAWG2_onlyconverged_coverage_beta_v2')
opt_all_datasets[['B1diag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_3401_18_6_0_diagREDM_betainterceptPCAWG2_betaslopePCAWG2_covmatPCAWG2_onlyconverged_coverage_beta_v2')
opt_all_datasets[['B2full']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_14072_80_6_0_fullREDM_betainterceptPCAWG4_betaslopePCAWG4_covmatPCAWG4_onlyconverged_coverage_beta_v2')
opt_all_datasets[['B2diag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_14072_80_6_0_diagREDM_betainterceptPCAWG4_betaslopePCAWG4_covmatPCAWG4_onlyconverged_coverage_beta_v2')
opt_all_datasets[['B3full']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_14072_87_6_0_fullREDM_betainterceptPCAWG4_betaslopePCAWG4_covmatFULLPCAWG4_onlyconverged_coverage_beta_v3')
opt_all_datasets[['B3diag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_14072_87_6_0_diagREDM_betainterceptPCAWG4_betaslopePCAWG4_covmatFULLPCAWG4_onlyconverged_coverage_beta_v3')
opt_all_datasets[['B4full']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_nPCAWG6_nlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_fullREDM_betainterceptPCAWG6_betaslopePCAWG6_covmatFULLPCAWG6_onlyconverged_bias_betas_v2')
opt_all_datasets[['B4diag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_nPCAWG6_nlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_diagREDM_betainterceptPCAWG6_betaslopePCAWG6_covmatFULLPCAWG6_onlyconverged_bias_betas_v2')
opt_all_datasets[['B4full_lownlambda']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_nPCAWG6_lownlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_fullREDM_betainterceptPCAWG6_betaslopePCAWG6_covmatFULLPCAWG6_onlyconverged_bias_betas_v2')
opt_all_datasets[['B4diag_lownlambda']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_nPCAWG6_lownlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_diagREDM_betainterceptPCAWG6_betaslopePCAWG6_covmatFULLPCAWG6_onlyconverged_bias_betas_v2')
opt_all_datasets[['B4full_low2nlambda']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_nPCAWG6_low2nlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_fullREDM_betainterceptPCAWG6_betaslopePCAWG6_covmatFULLPCAWG6_onlyconverged_bias_betas_v2')
opt_all_datasets[['B4diag_low2nlambda']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_nPCAWG6_low2nlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_diagREDM_betainterceptPCAWG6_betaslopePCAWG6_covmatFULLPCAWG6_onlyconverged_bias_betas_v2')
sapply(opt_all_datasets, function(i) length(i$input_all)) ## if zero, no files were found
## ------------------------------------------------------------------- ##

## ------------------------------------------------------------------- ##
for(i in opt_all_datasets){
  ## for each non-converged run, check if we have converged files from a run, remove the non-converged files from the same fun
  sapply(grep('_NC.RDS', i$input_all, value = T), function(j){
    if(gsub("_NC", "", j) %in% i$input_all){
      # print('Present')
      # print(j)
      # print(i$input_all[i$input_all == gsub("_NC", "", j)])
      newfile <- gsub('../../../data/assessing_models_simulation/inference_results/TMB/nlminb/', ' ~/Desktop/DM_revisions/obsolete_runs/', j)
      system(paste0('mkdir -p ', dirname(newfile)))
      system(paste0('mv ', j, newfile))
      cat('\n')
    }
  })
}
## ------------------------------------------------------------------- ##

## ------------------------------------------------------------------- ##

return_df_plots_bias_coverage <- function(opt, path_to_data="../../../data/"){
  
  for(opt_param in c('d', 'lambda', 'n', 'nlambda')){
    opt_param_as_numeric <- as.numeric(opt[[opt_param]])
    if(is.na(opt_param_as_numeric)){
      ## is character: read file
      opt_param_file = paste0(path_to_data, 'assessing_models_simulation/additional_files/multiple_fixed_',
                              opt[[opt_param]], '.RDS')
      if(opt_param == 'nlambda'){
        ## it should be a list
        Nm_lambda = lapply(readRDS(opt_param_file), as.numeric)
      }else{
        assign(opt_param, as.numeric(readRDS(opt_param_file))) ## assigning the variable to the read value
      }
      if(opt_param == 'lambda'){
        stopifnot(length(lambda) == 2)
      }
      if(opt_param != 'nlambda'){
        cat('Parameter ', opt_param, '=', get(opt_param), '\n')
      }else{
        cat('Parameter ', opt_param, ' read\n')
      }
    }else{
      if(opt_param == 'lambda'){
        ## in the case of lambda, it might contain different lambdas if read.
        ## if numeric, lambda is one single value
        lambda = rep(opt$lambda, 2) ## overdispersion scalars. Lower value -> higher overdispersion
      }else{
        assign(opt_param, as.numeric(opt[[opt_param]])) ## assigning the variable to the read value
        # d = opt$d ## number of signatures
        # n = opt$n ## number of samples
        if(opt_param == 'nlambda'){
          Nm_lambda = opt$nlambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)
        }
        
      }    
      
    }
  }
  
  input_run <- sub("_dataset[^_dataset]+$", "", basename(opt$input))
  if(length(table(input_run)) > 1){
    warning('There are runs of mutiple parameters included in this analysis\n')
  }
  
  plot_only_converged <- !opt$run_nonconverged
  
  if(plot_only_converged){
    add_convergence <- '_onlyconverged_'
  }else{
    add_convergence <- 'withnonconverged_'
  }
  
  model <- opt$model
  if(opt$multiple_runs){
    name_dataset <- paste0('multiple_', opt$dataset_generation, '_') 
  }else{
    name_dataset <- paste0(opt$dataset_generation, '_') 
  }
  
  appender <- function(string){
    sapply(string, function(stringb){
      if(stringb == 'Intercept'){
        TeX(paste("$\\beta_0$")) 
      }else if (stringb == 'Slope'){
        TeX(paste("$\\beta_1$")) 
      }
    })
  }
  
  # first_part_output <- paste0("../results/results_TMB/simulated_datasets/bias_and_coverage-replot/",
  #                             # "multiple_", 
  #                             name_dataset, opt$optimiser, '_', opt$n, '_', opt$nlambda,  '_', opt$lambda,  '_', opt$d,
  #                             '_', opt$beta_gamma_shape,  '_', model,  '_', idx_dataset_betaintercept, '_',
  #                             idx_dataset_betaslope, '_', idx_dataset_cov, "/setsim_", 
  #                             name_dataset, opt$optimiser, '_', opt$n, '_', opt$nlambda,  '_', opt$lambda,  '_', opt$d,
  #                             '_', opt$beta_gamma_shape,  '_', model,  '_', idx_dataset_betaintercept, '_',
  #                             idx_dataset_betaslope, '_', idx_dataset_cov, add_convergence)
  # system(paste0('mkdir -p ', first_part_output))
  
  # name_dataset0 <- paste0(strsplit(name_dataset, '_')[[1]][1:2], sep = '_', collapse = '')
  
  ## Load runs
  print(opt$input)
  cat("Loading runs\n")
  
  lst <- opt$input
  
  all_pd <- lapply(lst, function(i){x <- readRDS(i); try(x$pdHess)})
  all_pd[sapply(all_pd, typeof) == 'character'] = FALSE
  all_pd_list <- as.vector(unlist(all_pd))
  
  table(all_pd_list)
  
  lst0 <- lst
  if(plot_only_converged) lst = lst[all_pd_list]
  runs <- lapply(lst, readRDS)
  tryerror <- sapply(runs, typeof) == "list"
  runs <- runs[tryerror]
  lst <- lst[tryerror]
  
  # print(paste0(path_to_data, "assessing_models_simulation/datasets/",
  #              gsub(paste0(opt$model, "_"), "", basename(lst[1]))))
  
  if(opt$dataset_generation %in% c("GenerationJnormMax", "GenerationJnormInv")){
    ## as we do not have the true values, because we have changed the column order (for betas we could easily softmax and re-take the ALR, but for the covariances it would be more difficult)
    # ggplot()
    # ggsave(paste0(first_part_output, "coverage_beta.pdf"), width = 2.5, height = 2.8)
  }else{
    
    ### load datasets
    cat("Loading datasets...\n")
    x <- lapply(lst, function(i) (readRDS(paste0(path_to_data, "assessing_models_simulation/datasets/",
                                                 gsub("_dataset.*", "", gsub(paste0(opt$model, "_"), "", basename(i))), '/',
                                                 gsub(paste0(opt$model, "_"), "", basename(i))))))
    cat("... datasets loaded.\n")
    
    cat("Creating summaries of runs\n")
    summaries = lapply(runs, function(i){
      summary <- summary.sdreport(i)
      summary <- summary[!grepl("u_large", rownames(summary)),]
      summary <- summary[!grepl("u_random_effects", rownames(summary)),] ## for  "singleREDM"
      summary
    })
    
    summaries2 = sapply(summaries, function(j) j[,1])
    
    cat("Melting summaries of runs\n")
    
    summaries_melt = data.frame(melt(summaries2), stringsAsFactors = F)
    print(summaries_melt)
    summaries_melt$Var1 = as.character(summaries_melt$Var1)
    summaries_melt$idx_param = rep(1:sum(summaries_melt$Var2 == 1), length(lst))
    
    cat("Transforming estimates of runs\n")
    
    ## log_lambda simply need to be exponentiated (e^x) to compare it to the simulated value
    summaries_melt[summaries_melt$Var1 == "log_lambda", "value"] = exp(summaries_melt[summaries_melt$Var1 == "log_lambda", "value"])
    
    ## logs_sd_RE simply need to be exponentiated (e^x) to compare it to the simulated value
    ## if we want the variances, we need it to the power of two
    if("logs_sd_RE" %in% unique(summaries_melt$Var1)){
      summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"] = (exp(summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"]))
      summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"] = exp(summaries_melt[summaries_melt$Var1 == "logs_sd_RE",
                                                                                        "value"])**2
    }
    
    d <- length(python_like_select_name(runs[[1]]$par.fixed, 'logs_sd_RE'))
    if("cov_par_RE" %in% unique(summaries_melt$Var1)){
      cat('Getting estimated covariances\n')
      for(idx in unique(summaries_melt$Var2)){
        .x <- L_to_cov(summaries_melt[(summaries_melt$Var1 == "cov_par_RE") & (summaries_melt$Var2 == idx), "value"], d=d)
        summaries_melt[(summaries_melt$Var1 == "cov_par_RE") & (summaries_melt$Var2 == idx), "value"] <- .x[upper.tri(.x)]
      }
      ## checking positive semi-definiteness
      give_first_col <- function(i){
        if(is.null(dim(i))){i[1]}else{i[,1]}
      }
      mvtnorm::rmvnorm(n=10, mean = rep(0, d), sigma = L_to_cov(give_first_col(python_like_select_rownames(summaries[[1]], 'cov_par_RE')), d=d))
      
    }
    
    summaries_melt[summaries_melt$Var1 == "log_lambda", "Var1"] = "lambda"
    summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "Var1"] = "var_RE"
    
    
    if(opt$dataset_generation %in% c("GenerationJnorm")){
      ## full, unconstrained, random effects matrix
      cat('Reading covariance matrix\n')
      covmat_path = paste0(path_to_data, "assessing_models_simulation/additional_files/multiple_fixed_", opt$cov, ".RDS")
      cat('Reading ', covmat_path, '\n')
      cov_mat_true = readRDS(covmat_path)
      sds <- diag(cov_mat_true)
      covs <- cov_mat_true[upper.tri(cov_mat_true)]
      covs_true <- TRUE
    }else{
      stop('Specify correct <dataset_generation>\n')
    }
    
    if(model %in% c("fullREM")){
      true_vals = c(as.vector(x[[1]]$beta), ## betas
                    rep(0,((x[[1]]$d-1)**2-(x[[1]]$d-1))/2), ## covariances RE
                    sds) ##sd RE ## this is particular to GenerationCnorm
      # }else if(model %in% c("diagREM")){
      #   true_vals = c(as.vector(x[[1]]$beta), ## betas
      #                 sds) ##sd RE ## this is particular to GenerationCnorm
    }else if(model %in% c("fullREDM", "fullREDMsinglelambda")){
      if(model == "fullREDM"){
        overdisp <- x[[1]]$lambda
        if(opt$dataset_generation %in% c("GenerationMGnorm", "generationMGnorm")){
          overdisp <- c(NA,NA)
        }
      }else if(model=="fullREDMsinglelambda"){
        if(!is.na(x[[1]]$lambda)){
          stopifnot(x[[1]]$lambda[1] == x[[1]]$lambda[2])
        }
        overdisp <- x[[1]]$lambda[1]
      }
      
      if(covs_true){
        ## we have covariances
        covs <- covs
      }else{
        ## if we don't have covariances, we put zeros in the off-diagonal elements of the covariance matrix
        covs <- rep(0,((x[[1]]$d-1)**2-(x[[1]]$d-1))/2)
      }
      true_vals = c(as.vector(x[[1]]$beta), ## betas
                    sds, ##sd RE ## this is particular to GenerationCnorm
                    covs, ## covariances RE
                    overdisp)
    }else if(model %in% c("diagREDM", "diagREDMsinglelambda" )){
      if(model == "diagREDM"){
        if(opt$dataset_generation %in% c("GenerationMGnorm", "generationMGnorm")){
          overdisp = c(NA, NA)
        }else{
          overdisp <- x[[1]]$lambda
        }
      }else if(model=="diagREDMsinglelambda"){
        stopifnot(x[[1]]$lambda[1] == x[[1]]$lambda[2])
        overdisp <- x[[1]]$lambda[1]
      }
      true_vals = c(as.vector(x[[1]]$beta), ## betas
                    sds, ##sd RE ## this is particular to GenerationCnorm
                    overdisp)
    }else if(model=="singleREDM"){
      if(opt$dataset_generation %in% c("GenerationMGnorm", "generationMGnorm")){
        overdisp = c(NA, NA)
      }else{
        overdisp <- x[[1]]$lambda
      }
      true_vals = c(as.vector(x[[1]]$beta), ## betas
                    NA, ## single logsd for the single random effects. They are not simulated like this in most datasets, so I add NA. If I am using a generation where random effects are simulated with a single value, I should add it here
                    overdisp)
    }else if(model %in% c("FEDMsinglelambda" )){
      overdisp <- x[[1]]$lambda
      true_vals = c(as.vector(x[[1]]$beta), ## betas
                    overdisp)
    }else if(model == "fullREDMonefixedlambda"){
      true_vals <- c(as.vector((x[[1]]$beta)),
                     rep(NA, length(python_like_select_name(runs[[1]]$par.fixed, 'cov_par_RE'))),
                     sds, c((x[[1]]$lambda)))
    }else if(model == "fullREDMonefixedlambda3"){
      true_vals <- c(as.vector((x[[1]]$beta)),
                     rep(NA, length(python_like_select_name(runs[[1]]$par.fixed, 'cov_par_RE'))),
                     sds, c((x[[1]]$lambda)))
    }else{
      stop('<Specify correct model>')
    }
    
    summaries_melt$true = rep(true_vals, length(lst))
    summaries_melt$subtract = summaries_melt$value - summaries_melt$true
    
    # plot(summaries_melt$true[summaries_melt$Var1 == 'beta'],
    #      summaries_melt$value[summaries_melt$Var1 == 'beta']) ### correlation between true and estimated beta
    # abline(coef = c(0,1), lty='dashed', col='blue')
    
    # cat('Creating ', paste0(first_part_output, "betacorrelation.pdf"), '\n')
    # ggplot(cbind.data.frame(summaries_melt[summaries_melt$Var1 == 'beta',], intslope=c('Intercept', 'Slope') ),
    #        aes(x=true, y=value, shape=intslope))+geom_point()+
    #   geom_abline(intercept = 0, slope = 1, lty='dashed', col='black')+theme_bw()+labs(shape="", x='True beta coefficient', y='Estimated beta coefficient')+
    #   theme(legend.position = "bottom")#+scale_color_manual(values = c(''))
    # ggsave(paste0(first_part_output, "betacorrelation.pdf"), width = 2.5, height = 2.8)
  
    # cat('Creating ', paste0(first_part_output, "bias_betas_v2.pdf"), '\n')
    plot_bias <- ggplot(cbind.data.frame(summaries_melt[summaries_melt$Var1 == 'beta',], 
                            intslope=c('Intercept', 'Slope')),
           aes(x=(1+idx_param) %/% 2, y=subtract, group=idx_param))+
      geom_abline(slope = 0, intercept = 0, lty='dashed', col='blue')+
      geom_boxplot()+
      labs(x="Beta", y="Bias", lty="")+
      theme_bw()+
      facet_wrap(.~intslope, drop = T, labeller = as_labeller(appender, 
                                                              default = label_parsed))+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+theme(legend.position = "bottom")
    # ggsave(paste0(first_part_output, "bias_betas_v2.pdf"), width = 2.2, height = 2)
    # cat('... done\n')
    
    confints <- sapply(summaries, function(it_run){
      sapply(1:length(true_vals), function(i){
        confintint = c(it_run[i,1]-1.96*it_run[i,2],it_run[i,1]+1.96*it_run[i,2])
        (true_vals[i] > confintint[1]) & (true_vals[i] < confintint[2])
      } )
    })
    rownames(confints) = rownames(summaries[[1]])
    
    df_coverage <- cbind.data.frame(parameter=make.names(rownames(confints), unique = T), CI=apply(confints, 1, mean, na.rm=T))
    df_coverage$type_param = gsub("\\..*", "", df_coverage$parameter)
    df_coverage$type_param[df_coverage$type_param == "beta"] = c('beta_intercept', 'beta_slope')
    
    df_coverage$type_param[df_coverage$type_param == "beta_intercept"] <- "Intercept"
    df_coverage$type_param[df_coverage$type_param == "beta_slope"] <- "Slope"

    draws_binom <- rbinom(n = 1000, size = sum(all_pd_list), prob = 0.95)
    hist(draws_binom)
    var(draws_binom)
    pbinomprobs <- data.frame(t(quantile(draws_binom, probs = c(0.025, 0.975))/ sum(all_pd_list)))
    colnames(pbinomprobs) <- c('x1', 'x2')
    
    df_coverage[df_coverage$parameter == 'beta','parameter'] <- 'beta.0'
    
    pbinomprobs <- data.frame(t(apply(df_coverage[df_coverage$type_param %in% c("Intercept",  "Slope"),c('parameter', 'type_param')],
                                      1, function(i) unlist(c(i, pbinomprobs)))))
    colnames(pbinomprobs) <- c('parameter', 'type_param', 'x1', 'x2')
    pbinomprobs$x1 <- as.numeric(pbinomprobs$x1)
    pbinomprobs$x2 <- as.numeric(pbinomprobs$x2)
    
    transformparam <- function(parameter) 1+as.numeric(gsub('.*[.]', '', parameter)) %/% 2
    
    if(min(df_coverage[grepl('beta', df_coverage$parameter),'CI']) < 0.3){
      min_y = min(df_coverage[grepl('beta', df_coverage$parameter),'CI'])
      min_y_2 = min(df_coverage[grepl('beta', df_coverage$parameter),'CI'])
    }else{
      min_y = 0.3
      min_y_2 = 0.7
    }
    
    A <- ggplot(df_coverage[df_coverage$type_param %in% c("Intercept"),])+
      geom_ribbon(data = pbinomprobs[pbinomprobs$type_param %in% c("Intercept"),], aes(x=transformparam(parameter), ymin=x1, ymax=x2, group=type_param),
                  fill='yellow', alpha=0.9)+
      geom_abline(slope = 0, intercept = 0.95, col='blue', lty='dashed')+
      geom_line(aes(x=(1+as.numeric(gsub('.*[.]', '', parameter)) %/% 2), y=CI,
                    group=1))+geom_point(aes(x=transformparam(parameter), y=CI,
                                             group=1))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_bw()+#lims(y=c(0,1))+
      theme(legend.position = "bottom")+
      facet_wrap(.~type_param, drop = T, labeller = as_labeller(appender, 
                                                                default = label_parsed))+
      labs(shape="", x='Beta', y='Coverage')+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      ylim(min=c(min_y, 1))
    # if(min_y == 0.3){
    #   A <- A+annotate(geom="label",x=-Inf,y=0.315,label="   //", fill="white", col='black', label.size = NA)+
    #     annotate(geom="label",x=Inf,y=0.315,label="//   ", fill="white", col='black', label.size = NA)+
    #     geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=0.31, ymax=0.32), col='black', fill="white", alpha=0.2, lty=1)
    # }
    B <- ggplot(df_coverage[df_coverage$type_param %in% c("Slope"),])+
      geom_ribbon(data = pbinomprobs[pbinomprobs$type_param %in% c("Slope"),], aes(x=transformparam(parameter), ymin=x1, ymax=x2, group=type_param),
                  fill='yellow', alpha=0.9)+
      geom_abline(slope = 0, intercept = 0.95, col='blue', lty='dashed')+
      geom_line(aes(x=(1+as.numeric(gsub('.*[.]', '', parameter)) %/% 2), y=CI,
                    group=1))+geom_point(aes(x=transformparam(parameter), y=CI,
                                             group=1))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_bw()+#lims(y=c(0,1))+
      theme(legend.position = "bottom")+
      scale_y_continuous(position = "right", limits = c(min_y_2, 1))+
      facet_wrap(.~type_param, drop = T, labeller = as_labeller(appender, 
                                                                default = label_parsed))+
      labs(shape="", x='Beta', y='Coverage')+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), axis.title.y=element_blank())
    # if(min_y == 0.3){
    #   B <- B+  annotate(geom="label",x=-Inf,y=0.715,label="   //", fill="white", col='black', label.size = NA)+
    #     annotate(geom="label",x=Inf,y=0.715,label="//   ", fill="white", col='black', label.size = NA)+
    #     geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=0.71, ymax=0.715), col='black', fill="white", alpha=0.2, lty=1)
    # }
    
    # try(dev.off())
    # cat('Creating ', paste0(first_part_output, "coverage_beta_v3.pdf"), '\n')
    # pdf(paste0(first_part_output, "coverage_beta_v3.pdf"), width = 2.2, height = 2)
    print(cowplot::plot_grid(A, B, nrow=1, rel_widths = c(1.1, 0.95)))
    # dev.off()
  return(list(bias_df=plot_bias, coverage_df=list(A,B)))
    
  }
}

# return_df_plots_bias_coverage(opt_all_datasets$A1diag) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$A1full) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$A1single) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$A2diag) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$A2full) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$A2single) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$A3diag) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$A3full) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$A3single) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$baselinediag) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B1diag) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B1full) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B2diag) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B2full) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B3diag) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B3full) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B4diag) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B4full) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B4full_lownlambda) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B4diag_lownlambda) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B4full_low2nlambda) ## pass
# return_df_plots_bias_coverage(opt_all_datasets$B4diag_low2nlambda) ## pass

# opt_all_datasets_subset <- opt_all_datasets[!grepl("B4|B2", names(opt_all_datasets))] ## already done
# df_plots_bias_coverage_subset <- lapply(opt_all_datasets_subset, function(i) try(return_df_plots_bias_coverage(i)))
df_plots_bias_coverage_all <- lapply(opt_all_datasets, function(i) try(return_df_plots_bias_coverage(i)))


## maybe there should be two shared limits, one for intercept and one for slope - https://stackoverflow.com/questions/51735481/ggplot2-change-axis-limits-for-each-individual-facet-panel
common_lims <- list(
                    baselinediag = list(bias=c(-0.35, 0.35), coverage=list(c(0.71, 1), c(0.71, 1))), ## only one run
                    # A1 =list(bias=c(-3, 7)),
                    # A2 =list(bias=c(-2.1, 1.8)),
                    # A3 =list(bias=c(-2.5, 1.5)),
                    A = list(bias = c(-3, 7), coverage=list(c(0,1), c(0.25,1))),
                    B1 =list(bias=c(-4.2, 1), coverage=list(c(0.71,1), c(0.71,1))),
                    B2=list(bias=c(-0.4, 0.4), coverage=list(c(0.71,1), c(0.71,1))),
                    B3=list(bias=c(-.5, .3), coverage=list(c(0.51, 1), c(0.71,1))),
                    B4 = list(bias = c(-1, 1), coverage=list(c(0.71,1), c(0.61,1)))
)
names(df_plots_bias_coverage_all)

for(sim_name in names(common_lims)){
  require(gridExtra)
  
  cat('Sim name: ', sim_name, '\n')
  
  idx_sim = grep(sim_name, names(df_plots_bias_coverage_all))
  for(i in idx_sim){
    if(sim_name %in% c('A1', 'A2', 'A3')){
      sim_name2 = 'A'
    }else{
      sim_name2 = sim_name
    }
    ## bias
    df_plots_bias_coverage_all[[i]]$bias_df = df_plots_bias_coverage_all[[i]]$bias_df+lims(y=common_lims[[sim_name2]]$bias)
    
    
    add_line <- function(input, idx){
      input$coverage_df[[idx]] +
        annotate(geom="label",x=-Inf,y=common_lims[[sim_name2]]$coverage[[idx]][1]+0.01,label="   //", fill="white", col='black', label.size = NA)+
        annotate(geom="label",x=Inf,y=common_lims[[sim_name2]]$coverage[[idx]][1]+0.01, label="//   ", fill="white", col='black', label.size = NA)+
        geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, 
                              ymin=common_lims[[sim_name2]]$coverage[[idx]][1]+0.005,
                              ymax=common_lims[[sim_name2]]$coverage[[idx]][1]+0.015),
                  col='black', fill="white", alpha=0.2, lty=1)
    }
    ## coverage
    for(idx in c(1,2)){
      if(common_lims[[sim_name2]]$coverage[[idx]][1] != 0){
        df_plots_bias_coverage_all[[i]]$coverage_df[[idx]] = add_line(df_plots_bias_coverage_all[[i]], idx=idx)
      }
      if(idx == 1){
        df_plots_bias_coverage_all[[i]]$coverage_df[[idx]] = df_plots_bias_coverage_all[[i]]$coverage_df[[idx]]+lims(y=common_lims[[sim_name2]]$coverage[[idx]])
      }else if(idx == 2){
        df_plots_bias_coverage_all[[i]]$coverage_df[[idx]] = df_plots_bias_coverage_all[[i]]$coverage_df[[idx]]+ 
          scale_y_continuous(position = "right", limits = common_lims[[sim_name2]]$coverage[[idx]])+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(), axis.title.y=element_blank())
      }
    }
  }
  
  # do.call('grid.arrange', lapply(df_plots_bias_coverage_all[idx_sim], `[[`, 'bias_df'))
  # do.call('grid.arrange', lapply(lapply(df_plots_bias_coverage_all[idx_sim], `[[`, 'coverage_df'), `[[`, 1))
  for(i in idx_sim){
    print(df_plots_bias_coverage_all[[i]]$bias_df)
    ggsave(paste0("../../../results/figures_paper/bias_coverage/", names(df_plots_bias_coverage_all)[i], "_bias.pdf"),
           width = 2.2, height = 2)
  }
  
  for(i in idx_sim){
    print(cowplot::plot_grid(df_plots_bias_coverage_all[[i]]$coverage_df[[1]],
                             df_plots_bias_coverage_all[[i]]$coverage_df[[2]], nrow=1, rel_widths = c(1.1, 0.95)))
    ggsave(paste0("../../../results/figures_paper/bias_coverage/", names(df_plots_bias_coverage_all)[i], "_coverage.pdf"),
           width = 2.2, height = 2)
  }
}

