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

find_input_from_list <- function(j){
  string_files <- paste0("../../../data/assessing_models_simulation/inference_results/TMB/", j$optimiser, "/",
         'multiple_', j$dataset_generation, '_', j$n, '_',
         j$nlambda, '_', j$lambda, '_', j$d, '_', j$beta_gamma_shape, '_', j$model, '_',
         j$betaintercept, '_',
         j$betaslope, '_', j$cov, '/')
  cat('Looking for files that match <', string_files, '>...\n')
  lst <- list.files(string_files, full.names = T)
  lst <- lst[grepl(paste0('multiple_', j$dataset_generation, '_'), lst)]
  lst <- lst[!grepl('_NC.RDS', lst)]
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
  return(opt)
}

opt_all_datasets <- list()
opt_all_datasets[['A1diag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_2_5_0_diagREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A1full']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_2_5_0_fullREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
opt_all_datasets[['A2diag']] <- opt_from_bias_fig('setsim_multiple_GenerationJnorm_nlminb_200_180_20_5_0_diagREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_bias_betas_v2')
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

sapply(opt_all_datasets, function(i) length(i$input)) ## if zero, no files were found

## ------------------------------------------------------------------- ##

return_df_plots_bias_coverage <- function(opt_list){
  
  for(opt_param in c('d', 'lambda', 'n', 'nlambda')){
    opt_param_as_numeric <- as.numeric(opt[[opt_param]])
    if(is.na(opt_param_as_numeric)){
      ## is character: read file
      opt_param_file = paste0('../data/assessing_models_simulation/additional_files/multiple_fixed_',
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
  
  first_part_output <- paste0("../results/results_TMB/simulated_datasets/bias_and_coverage-replot/",
                              # "multiple_", 
                              name_dataset, opt$optimiser, '_', opt$n, '_', opt$nlambda,  '_', opt$lambda,  '_', opt$d,
                              '_', opt$beta_gamma_shape,  '_', model,  '_', idx_dataset_betaintercept, '_',
                              idx_dataset_betaslope, '_', idx_dataset_cov, "/setsim_", 
                              name_dataset, opt$optimiser, '_', opt$n, '_', opt$nlambda,  '_', opt$lambda,  '_', opt$d,
                              '_', opt$beta_gamma_shape,  '_', model,  '_', idx_dataset_betaintercept, '_',
                              idx_dataset_betaslope, '_', idx_dataset_cov, add_convergence)
  system(paste0('mkdir -p ', first_part_output))
  
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
  
  print(paste0("../data/assessing_models_simulation/datasets/",
               gsub(paste0(opt$model, "_"), "", basename(lst[1]))))
  
  if(opt$dataset_generation %in% c("GenerationJnormMax", "GenerationJnormInv")){
    ## as we do not have the true values, because we have changed the column order (for betas we could easily softmax and re-take the ALR, but for the covariances it would be more difficult)
    ggplot()
    ggsave(paste0(first_part_output, "coverage_beta.pdf"), width = 2.5, height = 2.8)
  }else{
    
    ### load datasets
    cat("Loading datasets...\n")
    x <- lapply(lst, function(i) (readRDS(paste0("../data/assessing_models_simulation/datasets/",
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
      cov_mat_true = readRDS(paste0("../data/assessing_models_simulation/additional_files/multiple_fixed_", idx_dataset_cov, ".RDS"))
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
    
    plot(summaries_melt$true[summaries_melt$Var1 == 'beta'],
         summaries_melt$value[summaries_melt$Var1 == 'beta']) ### correlation between true and estimated beta
    abline(coef = c(0,1), lty='dashed', col='blue')
    
    cat('Creating ', paste0(first_part_output, "betacorrelation.pdf"), '\n')
    ggplot(cbind.data.frame(summaries_melt[summaries_melt$Var1 == 'beta',], intslope=c('Intercept', 'Slope') ),
           aes(x=true, y=value, shape=intslope))+geom_point()+
      geom_abline(intercept = 0, slope = 1, lty='dashed', col='black')+theme_bw()+labs(shape="", x='True beta coefficient', y='Estimated beta coefficient')+
      theme(legend.position = "bottom")#+scale_color_manual(values = c(''))
    # ggsave(paste0(first_part_output, "betacorrelation.pdf"), width = 2.5, height = 2.8)
  
    cat('Creating ', paste0(first_part_output, "bias_betas_v2.pdf"), '\n')
    ggplot(cbind.data.frame(summaries_melt[summaries_melt$Var1 == 'beta',], 
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
    if(min_y == 0.3){
      A <- A+annotate(geom="label",x=-Inf,y=0.315,label="   //", fill="white", col='black', label.size = NA)+
        annotate(geom="label",x=Inf,y=0.315,label="//   ", fill="white", col='black', label.size = NA)+
        geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=0.31, ymax=0.32), col='black', fill="white", alpha=0.2, lty=1)
    }
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
    if(min_y == 0.3){
      B <- B+  annotate(geom="label",x=-Inf,y=0.715,label="   //", fill="white", col='black', label.size = NA)+
        annotate(geom="label",x=Inf,y=0.715,label="//   ", fill="white", col='black', label.size = NA)+
        geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=0.71, ymax=0.715), col='black', fill="white", alpha=0.2, lty=1)
    }
    
    # try(dev.off())
    cat('Creating ', paste0(first_part_output, "coverage_beta_v3.pdf"), '\n')
    # pdf(paste0(first_part_output, "coverage_beta_v3.pdf"), width = 2.2, height = 2)
    print(cowplot::plot_grid(A, B, nrow=1, rel_widths = c(1.1, 0.95)))
    # dev.off()
  return(list(bias_df=, coverage_df=))
    
  }
}

# for(generation in c(
#   # 'GenerationHnormtwolambdas', 'GenerationJnorm', 'GenerationK',
#   # 'GenerationMixturefewersignaturespairedKidneyRCCPCAWG', 'GenerationMixturefewersignaturespairedPCAWG',
#   #                   'GenerationMixturefewersignaturespairedstomachPCAWG', 'GenerationMixturefewersignaturesPCAWG',
#   #                   'GenerationMixturefewersmallsignaturespairedKidneyRCCPCAWG',
#   #                   'GenerationMixturefewersignaturespairedProstAdenoCAPCAWG',
#   #                   'GenerationMixturefewersignaturespairedCNSGBMPCAWG',
#   #                   'GenerationMixturefewersignaturespairedObsNmPancEndocrinePCAWG',
#   #                   'GenerationMixturefewersignaturespairedObsNmUterusAdenoCAPCAWG',
#   # 'GenerationJnormBTwoLambdasOneChangingBeta',
#   # 'GenerationPois',
#   'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGProstAdenoCAPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMOvaryAdenoCAPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMLungSCCPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMKidneyRCCpapillaryPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMPancEndocrinePCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMColoRectAdenoCAPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMLymphBNHLPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMHeadSCCPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMEsoAdenoCAPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMLymphCLLPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMPancAdenoCAPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMKidneyChRCCPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmObsDMCNSGBMPCAWG',
#   'GenerationMixtureallsignaturespairedObsNmCNSGBMPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmPoissonSigPCAWGColoRectAdenoCAPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGEsoAdenoCAPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGHeadSCCPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGLungSCCPCAWG',
#   'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGKidneyChRCCPCAWG'
#                     )){
#   
#   try(rm(a))
#   a <- readRDS(paste0("../../../data/assessing_models_simulation/summaries_synthetic_DA/", generation, ".RDS"))
#   
#   all_converged <- !apply(a$pvals_data_frame, 1, function(i) any(is.na(i)))
#   a$pvals_data_frame_adj <- a$pvals_data_frame
#   
#   a$pvals_data_frame_adj <- adjust_all(a$pvals_data_frame_adj, new_version = T, method = "bonferroni")
#   varying_betashape <-give_accuracies_with_varying_var(var = 'beta_gamma_shape',
#                                                        datasets_arg = a$datasets,
#                                                        pvals_data_frame_arg = a$pvals_data_frame)
#   varying_betashape_adj <-give_accuracies_with_varying_var(var = 'beta_gamma_shape',
#                                                        datasets_arg = a$datasets,
#                                                        pvals_data_frame_arg = a$pvals_data_frame_adj)
#   varying_betashapeAC <-give_accuracies_with_varying_var(var = 'beta_gamma_shape',
#                                                        datasets_arg = a$datasets[all_converged],
#                                                        pvals_data_frame_arg = a$pvals_data_frame[all_converged,])
#   
#   varying_betashape_n_adj <-give_accuracies_with_varying_var(var = c('beta_gamma_shape', 'n'), two_var = T,
#                                                        datasets_arg = a$datasets,
#                                                        pvals_data_frame_arg = a$pvals_data_frame_adj)
#   
#   varying_betashape_n <-give_accuracies_with_varying_var(var = c('beta_gamma_shape', 'n'), two_var = T,
#                                                              datasets_arg = a$datasets,
#                                                              pvals_data_frame_arg = a$pvals_data_frame)
#   varying_betashape_d <-give_accuracies_with_varying_var(var = c('beta_gamma_shape', 'd'), two_var = T,
#                                                          datasets_arg = a$datasets,
#                                                          pvals_data_frame_arg = a$pvals_data_frame)
#   
#   if((generation %in% c("GenerationMixturePCAWG", "GenerationMixturefewersignaturesPCAWG", "GenerationMixturefewersignaturespairedPCAWG")) | grepl('GenerationMixturefewersignaturespaired', generation) ){
#     varying_betashape$beta_gamma_shape <- signif(varying_betashape$beta_gamma_shape, 2)
#     varying_betashape_adj$beta_gamma_shape <- signif(varying_betashape_adj$beta_gamma_shape, 2)
#   }
#   
#   varying_betashape$model <- gsub("pvals_", "", varying_betashape$model)
#   varying_betashape_adj$model <- gsub("pvals_", "", varying_betashape_adj$model)
#   
#   .xx <- varying_betashape[which(varying_betashape$beta_gamma_shape == max(varying_betashape$beta_gamma_shape)),]
#   ## sometimes there are problems with numeric precision. if that is the case, select as the last value only one of the selected rows
#   .xx <- .xx[!duplicated(.xx$model),]
#   
#   remove_duplicated_rows <- function(.xx){
#     .xx <- .xx[!duplicated(.xx$model),]
#   }
#   
#   
#   if(generation == 'GenerationMixturefewersignaturespairedKidneyRCCPCAWG'){
#     title = 'Kidney-RCC (paired, larger sigs)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedPCAWG'){
#     title = 'Liver-HCC (paired, larger sigs)'
#   }else if ( generation  == 'GenerationMixturefewersignaturespairedstomachPCAWG'){
#     title = 'Stomach-AdenoCa (paired, larger sigs)'
#   }else if( generation == 'GenerationMixturefewersignaturesPCAWG'){
#     title = 'Liver-HCC (not paired, larger sigs)'
#   }else if( generation == 'GenerationMixturefewersmallsignaturespairedKidneyRCCPCAWG'){
#     title = 'Kidney-RCC (paired, small sigs)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedProstAdenoCAPCAWG'){
#     title = 'Prost-AdenoCA (paired, larger sigs)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedCNSGBMPCAWG'){
#     title = 'CNS-GBM (paired, larger sigs)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmPancEndocrinePCAWG'){
#     title = 'Panc-Endocrine (paired, larger sigs)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmUterusAdenoCAPCAWG'){
#     title = 'Uterus-AdenoCA (paired, larger sigs)'
#   }else if( generation ==  'GenerationMixtureallsignaturespairedObsNmCNSGBMPCAWG'){
#     title = 'CNS-GBM (paired, all sigs)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMOvaryAdenoCAPCAWG'){
#     title = 'Ovary-Adeno (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMLungSCCPCAWG'){
#     title = 'Lung-SCC (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMKidneyRCCpapillaryPCAWG'){
#     title = 'Kidney-RCCp (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMPancEndocrinePCAWG'){
#     title = 'Panc-Endocrine (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMColoRectAdenoCAPCAWG'){
#     title = 'ColoRect-AdenoCA (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMLymphBNHLPCAWG'){
#     title = 'Lymph-BNHL (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMHeadSCCPCAWG'){
#     title = 'Head-SCC (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMEsoAdenoCAPCAWG'){
#     title = 'Eso-Adeno (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMLymphCLLPCAWG'){
#     title = 'Lymph-CLL (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMPancAdenoCAPCAWG'){
#     title = 'Panc-AdenoCA (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMKidneyChRCCPCAWG'){
#     title = 'Kidney-ChRCC (paired, larger sigs, DM)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmObsDMCNSGBMPCAWG'){
#     title = 'CNS-GBM (paired, larger sigs, DM)'
#   }else if(generation == 'GenerationJnormBTwoLambdasOneChangingBeta'){
#     title = "Single category change"
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmPoissonSigPCAWGColoRectAdenoCAPCAWG'){
#     title = 'Colo-RectAdenoCA (paired, larger sigs, Poisson)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGEsoAdenoCAPCAWG'){
#     title = 'Eso-AdenoCA (paired, larger sigs, Norm)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGHeadSCCPCAWG'){
#     title = 'Head-SCC (paired, larger sigs, Gaussian)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGLungSCCPCAWG'){
#     title = 'Lung-SCC (paired, larger sigs, Gaussian)'
#   }else if( generation == 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWGKidneyChRCCPCAWG'){
#     title = 'Kidney-ChRCC (paired, larger sigs, Norm)'
#   }else{
#     title = generation
#   }
# 
#   varying_betashape$beta_gamma_shape <- signif(varying_betashape$beta_gamma_shape, 2)
#   varying_betashape_adj$beta_gamma_shape <- signif(varying_betashape_adj$beta_gamma_shape, 2)
#   
#   title_x <- 'Percentage of mixture'
#   if(generation %in% c("GenerationJnormBTwoLambdasOneChangingBeta", "GenerationPois", "GenerationJnorm", "GenerationHnormtwolambdas", "GenerationK")){
#     title_x <- TeX("$\\gamma$")
#   }
#   
#   ggplot(varying_betashape, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
#                                 lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_point()+geom_line()+theme_bw()+
#     scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
#     ggtitle(title)+
#     geom_label_repel(data = .xx,
#                      max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
#                      nudge_x=max(varying_betashape$beta_gamma_shape)+2,
#                      size=2.5)+
#     coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.3))+ ## used to be *1.5
#     theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
#   ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_palette2_factorv2.pdf"),
#          height = 3.0, width = 4.0)
#   
#   ggplot(varying_betashape_adj, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
#                                 lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_point()+geom_line()+theme_bw()+
#     scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
#     ggtitle(title)+
#     geom_label_repel(data = remove_duplicated_rows(varying_betashape_adj[which(varying_betashape_adj$beta_gamma_shape == max(varying_betashape_adj$beta_gamma_shape)),]),
#                      max.overlaps = Inf, aes(x=factor(max(varying_betashape_adj$beta_gamma_shape)),
#                                              lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')), direction = "y",
#                      nudge_x=max(varying_betashape_adj$beta_gamma_shape)+2,
#                      size=2.5)+
#     coord_cartesian(xlim = c(0, length(unique(varying_betashape_adj$beta_gamma_shape))*1.5))+
#     theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
#   ggsave(paste0(flder_out, generation, "/summaries/accuracy_adjbonferroni_with_betagammashape_palette2_factorv2.pdf"),
#          height = 3.0, width = 4.0)
#   
#   system(paste0("open ", flder_out, generation, "/summaries/"))
#   
#   ggplot(varying_betashapeAC, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
#                                 lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_point()+geom_line()+theme_bw()+
#     scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
#     ggtitle(title)+
#     geom_label_repel(data = .xx,
#                      max.overlaps = Inf, aes(x=factor(max(varying_betashapeAC$beta_gamma_shape))), direction = "y",
#                      nudge_x=max(varying_betashapeAC$beta_gamma_shape)+2,
#                      size=2.5)+
#     coord_cartesian(xlim = c(0, length(unique(varying_betashapeAC$beta_gamma_shape))*1.5))+
#     theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
#   
#   
#   varying_betashape_n_adj$beta_gamma_shape <- signif(varying_betashape_n_adj$beta_gamma_shape, 2)
#   varying_betashape_n$beta_gamma_shape <- signif(varying_betashape_n$beta_gamma_shape, 2)
#   varying_betashape_d$beta_gamma_shape <- signif(varying_betashape_d$beta_gamma_shape, 2)
#   ggplot(varying_betashape_n_adj, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
#                                 lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_point()+geom_line()+theme_bw()+
#     scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
#     ggtitle(title)+
#     geom_label_repel(data = .xx,
#                      max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
#                      nudge_x=max(varying_betashape$beta_gamma_shape)+2,
#                      size=2.5)+
#     coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.5))+
#     theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
#     facet_wrap(.~n)
#   
#   ggplot(varying_betashape_n_adj[varying_betashape_n_adj$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
#          aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
#                                       lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_point()+geom_line()+theme_bw()+
#     scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
#     ggtitle(title)+
#     geom_label_repel(data = .xx[varying_betashape_n_adj$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
#                      max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
#                      nudge_x=max(varying_betashape$beta_gamma_shape)+2,
#                      size=2.5)+
#     coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.5))+
#     theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
#     facet_wrap(.~n)
#   ggsave(paste0(flder_out, generation, "/summaries/accuracy_adjbonferroni_with_betagammashape_palette2_factorv2_subset.pdf"),
#          height = 3.0, width = 6.6)
#   
#   ggplot(varying_betashape_n[varying_betashape_n$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
#          aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
#              lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_point()+geom_line()+theme_bw()+
#     scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
#     ggtitle(title)+
#     geom_label_repel(data = .xx[varying_betashape_n$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
#                      max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
#                      nudge_x=max(varying_betashape$beta_gamma_shape)+2,
#                      size=2.5)+
#     coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.5))+
#     theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
#     facet_wrap(.~n)
#   ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_palette2_factorv2_subset.pdf"),
#          height = 3.0, width = 6.6)
#   
#   ggplot(varying_betashape_d[varying_betashape_d$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
#          aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
#              lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_point()+geom_line()+theme_bw()+
#     scale_color_manual(values=colours_models2)+labs(col='', x=title_x)+guides(col='none', lty='none')+
#     ggtitle(title)+
#     geom_label_repel(data = .xx[varying_betashape_d$model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM'),],
#                      max.overlaps = Inf, aes(x=factor(max(varying_betashape_d$beta_gamma_shape))), direction = "y",
#                      nudge_x=max(varying_betashape_d$beta_gamma_shape)+2,
#                      size=2.5)+
#     coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.5))+
#     theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
#     facet_wrap(.~d)
#   ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_facetD_palette2_factorv2_subset.pdf"),
#          height = 3.0, width = 6.6)
#   
#   ggplot(varying_betashape_n_adj[varying_betashape_n_adj$beta_gamma_shape == 0,],
#          aes(x=n, y = Accuracy, col=model, group=model, label=model,
#                                   lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_point()+geom_line()+theme_bw()+
#     scale_color_manual(values=colours_models2)+labs(col='', x='Number of samples')+
#     guides(lty='none')+
#     ggtitle(title)+
#     geom_label_repel(data = remove_duplicated_rows(varying_betashape_n_adj[which(varying_betashape_n_adj$beta_gamma_shape == 0),]),
#                      max.overlaps = Inf, aes(x=(max(varying_betashape_n_adj$n))), direction = "y",
#                      size=2.5)+
#     coord_cartesian(xlim = c(0, max(varying_betashape_n_adj$n)*1.5))+
#     theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
#   
#   pvals_noDA <- a$pvals_data_frame[as.vector(sapply(a$datasets, function(i) i$beta_gamma_shape)) == 0,]
#   pvals_DA <- a$pvals_data_frame[as.vector(sapply(a$datasets, function(i) i$beta_gamma_shape)) != 0,]
#   plot(pvals_noDA$pvals_diagREDM, pvals_noDA$pvals_fullREM); abline(coef = c(0,1))
#   
#   table(pvals_noDA$pvals_diagREDM < 0.05)
#   table(pvals_noDA$pvals_fullREM < 0.05)
#   
#   hist(pvals_noDA$pvals_diagREDM, breaks = 20)
#   hist(pvals_noDA$pvals_fullREM, breaks = 20)
#   
#   pvals_noDA_adj <- a$pvals_data_frame_adj[as.vector(sapply(a$datasets, function(i) i$beta_gamma_shape)) == 0,]
#   plot(pvals_noDA_adj$pvals_diagREDM, pvals_noDA_adj$pvals_fullREM); abline(coef = c(0,1), main=generation)
#   table(pvals_noDA_adj$pvals_diagREDM < 0.05)
#   table(pvals_noDA_adj$pvals_fullREDMSL < 0.05)
#   table(pvals_noDA_adj$pvals_fullREM < 0.05)
#   
# 
#   length_out_rocauc <- 50
#   rocauc <- lapply(seq(0, 1, length.out = length_out_rocauc), function(thresholdit) (as(t(sapply(colnames(a$pvals_data_frame), function(col_it){
#     summarise_DA_detection(true = a$pvals_data_frame$true, predicted = a$pvals_data_frame[,col_it] < thresholdit, verbose = F)[c('Specificity', 'Sensitivity')]
#   })), 'matrix')))
#   names(rocauc) <- seq(0, 1, length.out = length_out_rocauc)
#   # rocauc <- lapply(rocauc, as.matrix)
#   rocauc <- (reshape2::melt(rocauc, id.vars=c('Specificity', 'Sensitivity')))
#   rocauc <- dcast(rocauc, Var1+L1~Var2, value.var = "value")
#   rocauc$model <- gsub("pvals_", "", rocauc$Var1)
#   rocauc <- rocauc[rocauc$model != "true",]
#   unique(rocauc$model)
#   rocauc$model[rocauc$model == "ttest_ilr_adj"] = "ILR"
#   rocauc$model[rocauc$model == "ttest_props"] = "ttest"
#   
#   ggplot(rocauc, aes(x=1-Specificity, y= Sensitivity, col=model,
#                      lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_line()+theme_bw()+theme(legend.position = "bottom")+
#     scale_color_manual(values=colours_models2)+labs(col='', x='1- Specificity', y='Sensitivity')+guides(col='none', lty='none')
#   ggsave(paste0(flder_out, generation, "/summaries/ROCAUCcurve_palette.pdf"),
#          height = 3.0, width = 3)
#   
#   ggplot(rocauc, aes(x=1-Specificity, y= Sensitivity, col=model,
#                      lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_line()+theme_bw()+theme(legend.position = "bottom")+
#     scale_color_manual(values=colours_models2)+labs(col='', x='1- Specificity', y='Sensitivity')+guides(col='none', lty='none')+facet_wrap(.~model)
#   
#   ## for fixed cut-offs, changing effect sizes
#   ## information about effect sizes
#   beta_gamma_shapes <- as.vector(unlist(sapply(a$datasets, `[`, 'beta_gamma_shape')))
#   ns <- as.vector(unlist(sapply(a$datasets, `[`, 'n')))
#   thresholdit = 0.05
#   rocauc_2 <- lapply(sort(unique(ns)), function(betagammashapeit){
#     (as(t(sapply(colnames(a$pvals_data_frame), function(col_it){
#     summarise_DA_detection(true = a$pvals_data_frame$true[ns == betagammashapeit],
#                            predicted = a$pvals_data_frame[ns == betagammashapeit,col_it] < thresholdit, verbose = F)[c('Specificity', 'Sensitivity')]
#     })), 'matrix'))
#   })
#   names(rocauc_2) <- sort(unique(ns))
#   rocauc_2 <- (reshape2::melt(rocauc_2, id.vars=c('Specificity', 'Sensitivity')))
#   rocauc_2 <- dcast(rocauc_2, Var1+L1~Var2, value.var = "value")
#   rocauc_2$model <- gsub("pvals_", "", rocauc_2$Var1)
#   rocauc_2$model[rocauc_2$model == "ttest_ilr_adj"] = "ILR"
#   rocauc_2$model[rocauc_2$model == "ttest_props"] = "ttest"
#   ggplot(rocauc_2, aes(x=1-Specificity, y= Sensitivity, col=model,
#                      lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
#     geom_line()+theme_bw()+theme(legend.position = "bottom")+
#     scale_color_manual(values=colours_models2)+labs(col='', x='1- Specificity', y='Sensitivity')+guides(col='none', lty='none')
# 
#   
# }
# 
# a$pvals_data_frame$pvals_diagREDM[log(as.vector(sapply(a$datasets, function(i) i$beta_gamma_shape))) == -1.2]
