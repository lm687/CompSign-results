a <- readRDS("~/Desktop/comparison_betas_models_all.RDS")
b <- readRDS("~/Desktop/comparison_betas_models_nonexo.RDS")

table(a$fullRE_DMSL, a$diagRE_DMSL)
table(b$fullRE_DMSL)
table(b$diagRE_DMSL)

plot(a$fullRE_DMSL, a$diagRE_DMSL)

#######
library(reshape2)
library(ggplot2)
library(dplyr)

fullRE_DMSL_nonexo <- readRDS("~/Desktop/fullRE_DMSL_nonexoL.RDS")
diagRE_DMSL_nonexo <- readRDS("~/Desktop/diagRE_DMSL_nonexo.RDS")
fullRE_M_nonexo <- readRDS("~/Desktop/fullRE_M_nonexo.RDS")
model_fullRE_DMSL_list <- fullRE_DMSL_nonexo
model_diagRE_DMSL_list <- diagRE_DMSL_nonexo
model_fullRE_M_list <- fullRE_M_nonexo

enough_samples = read.table("~/Documents/PhD/GlobalDA/data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
renaming_pcawg <- read.table("~/Documents/PhD/GlobalDA/data/other/short_names_pcawg.txt", sep = "\t")
source("~/Documents/PhD/GlobalDA/code/2_inference_TMB/helper_TMB.R")


comparison_betas_models <- function(model_fullRE_DMSL_list, model_diagRE_DMSL_list, model_fullRE_M_list ){
  
  .x <- do.call('rbind.data.frame', lapply(enough_samples, function(ct){
    x_beta_fullRE_DMSL <- try(python_like_select_name(model_fullRE_DMSL_list[[ct]]$par.fixed, "beta"))
    x_beta_diagRE_DMSL <- try(python_like_select_name(model_diagRE_DMSL_list[[ct]]$par.fixed, "beta"))
    x_beta_fullRE_M <- try(python_like_select_name(model_fullRE_M_list[[ct]]$par.fixed, "beta"))
    if( (length(x_beta_fullRE_DMSL) != length(x_beta_diagRE_DMSL)) | (length(x_beta_fullRE_DMSL) != length(x_beta_fullRE_M)) ){
      ## if we don't have results for any, remove from the analysis
      list_betas <- list(x_beta_fullRE_DMSL, x_beta_diagRE_DMSL, x_beta_fullRE_M)
      typeofs_of_betas <- sapply(list_betas, typeof)
      if( all(typeofs_of_betas == "character")  ){
        return(NULL)
      }else{
        ## if we do have results for some, replace the error message by an NA string
        ## replace using the length of the first double entry
        
        if(typeofs_of_betas[1] == "character"){
          x_beta_fullRE_DMSL <- rep(NA, length(list_betas[[which(typeofs_of_betas == "double")[1]]]))
        }
        
        if(typeofs_of_betas[2] == "character"){
          x_beta_diagRE_DMSL <- rep(NA, length(list_betas[[which(typeofs_of_betas == "double")[1]]]))
        }
        
        if(typeofs_of_betas[3] == "character"){
          x_beta_fullRE_M <- rep(NA, length(list_betas[[which(typeofs_of_betas == "double")[1]]]))
        }
        
        if( (length(x_beta_fullRE_DMSL) != length(x_beta_diagRE_DMSL)) | (length(x_beta_fullRE_DMSL) != length(x_beta_fullRE_M)) ){
          warning(paste0(ct, ': the number of log-ratios is not consistent'))
          return(NULL)
        }
      }
      
    }else{
      cbind.data.frame(rbind.data.frame(
        cbind.data.frame(fullRE_DMSL=select_slope_2(x_beta_fullRE_DMSL),
                         diagRE_DMSL=select_slope_2(x_beta_diagRE_DMSL),
                         fullRE_M=select_slope_2(x_beta_fullRE_M),
                         beta_type='slope'),
        cbind.data.frame(fullRE_DMSL=select_intercept(x_beta_fullRE_DMSL),
                         diagRE_DMSL=select_intercept(x_beta_diagRE_DMSL),
                         fullRE_M=select_intercept(x_beta_fullRE_M),
                         beta_type='intercept')),
        ct=ct)
    }
  }))
  
  ## if something hasn't converged, remove the value
  bad_ct_fullRE_DMSL <- c(enough_samples[sapply(enough_samples,
                                                function(ct) (typeof(model_fullRE_DMSL_list[[ct]])) == 'character')],
                          enough_samples[sapply(enough_samples,
                                                function(ct)  try(model_fullRE_DMSL_list[[ct]]$pdHess)) != "TRUE"])
  bad_ct_diagRE_DMSL <- c(enough_samples[sapply(enough_samples,
                                                function(ct) (typeof(model_diagRE_DMSL_list[[ct]])) == 'character')],
                          enough_samples[sapply(enough_samples, function(ct)  try(model_diagRE_DMSL_list[[ct]]$pdHess)) != "TRUE"])
  bad_ct_fullRE_M <- c(enough_samples[sapply(enough_samples,
                                             function(ct) (typeof(model_fullRE_M_list[[ct]])) == 'character')],
                       enough_samples[sapply(enough_samples, function(ct)  try(model_fullRE_M_list[[ct]]$pdHess)) != "TRUE"])
  
  .x$fullRE_DMSL[(.x$ct %in% bad_ct_fullRE_DMSL)] = NA
  .x$diagRE_DMSL[(.x$ct %in% bad_ct_diagRE_DMSL)] = NA
  .x$fullRE_M[(.x$ct %in% bad_ct_fullRE_M)] = NA
  
  .x$ct2=renaming_pcawg[,2][match(.x$ct, renaming_pcawg[,1])]
  
  return(.x)
  
}


comparison_betas_models_all <- comparison_betas_models(model_fullRE_DMSL_list = fullRE_DMSL_nonexo,
                                                       model_diagRE_DMSL_list = diagRE_DMSL_nonexo,
                                                       model_fullRE_M_list = fullRE_M_nonexo)

plot(comparison_betas_models_all$diagRE_DMSL)


####
library(dplyr)
comparison_betas_models_rbind <- readRDS("~/Desktop/comparison_betas_models_rbind.RDS")

comparison_betas_models_rbind %>% filter(beta_type == 'Slope') %>%
  group_by(ct) %>%
  # filter(ct %in% "Liver-HCC") %>%
  summarise(
    # slope_diag_full_DMSL=try(coefficients(lm(y~x, data = cbind.data.frame(x=diagRE_DMSL, y=fullRE_M),
    #                                                  na.action = na.omit))[2]),
            slope_diag_full_DMSL=as.numeric(try(coefficients(lm(y~x, data = cbind.data.frame(x=fullRE_M, y=fullRE_DMSL),
                                                     na.action = na.omit))[2])))


comparison_betas_models_rbind_stats_per_ct <- rbind.data.frame(
  cbind.data.frame(comparison_betas_models_rbind %>% filter(beta_type == 'Slope') %>%
                     group_by(ct) %>%
                     summarise(rmse_diag_full_DMSL=sqrt(mean( (diagRE_DMSL-fullRE_DMSL)^2, na.rm = T )),
                               slope_diag_full_DMSL=as.numeric(try(coefficients(lm(y~x, data = cbind.data.frame(x=diagRE_DMSL, y=fullRE_DMSL), na.action = na.omit))[2])),
                               rmse_fullDMSL_fullM=sqrt(mean( (fullRE_M-fullRE_DMSL)^2, na.rm = T ))),
                   slope_fullDMSL_fullM=as.numeric(try(coefficients(lm(y~x, data = cbind.data.frame(x=fullRE_M, y=fullRE_DMSL), na.action = na.omit))[2])),
                   beta_type='Slope'),
  cbind.data.frame(comparison_betas_models_rbind %>% filter(beta_type == 'Intercept') %>%
                     group_by(ct) %>%
                     summarise(rmse_diag_full_DMSL=sqrt(mean( (diagRE_DMSL-fullRE_DMSL)^2, na.rm = T )),
                               slope_diag_full_DMSL=as.numeric(try(coefficients(lm(y~x, data = cbind.data.frame(x=diagRE_DMSL, y=fullRE_DMSL), na.action = na.omit))[2])),
                               rmse_fullDMSL_fullM=sqrt(mean( (fullRE_M-fullRE_DMSL)^2, na.rm = T ))),
                   slope_fullDMSL_fullM=as.numeric(try(coefficients(lm(y~x, data = cbind.data.frame(x=fullRE_M, y=fullRE_DMSL), na.action = na.omit))[2])),
                   beta_type='Intercept'))


######


comparison_randomintercepts_models <- function(model_fullRE_DMSL_list, model_diagRE_DMSL_list, model_fullRE_M_list ){
  
  .x <- do.call('rbind.data.frame', lapply(enough_samples, function(ct){
    x_RE_fullRE_DMSL <- try(python_like_select_name(model_fullRE_DMSL_list[[ct]]$par.random, "u_large"))
    x_RE_diagRE_DMSL <- try(python_like_select_name(model_diagRE_DMSL_list[[ct]]$par.random, "u_large"))
    x_RE_fullRE_M <- try(python_like_select_name(model_fullRE_M_list[[ct]]$par.random, "u_large"))
    if( (length(x_RE_fullRE_DMSL) != length(x_RE_diagRE_DMSL)) | (length(x_RE_fullRE_DMSL) != length(x_RE_fullRE_M)) ){
      ## if we don't have results for any, remove from the analysis
      list_RE <- list(x_RE_fullRE_DMSL, x_RE_diagRE_DMSL, x_RE_fullRE_M)
      typeofs_of_RE <- sapply(list_RE, typeof)
      if( all(typeofs_of_RE == "character")  ){
        return(NULL)
      }else{
        ## if we do have results for some, replace the error message by an NA string
        ## replace using the length of the first double entry
        
        if(typeofs_of_RE[1] == "character"){
          x_RE_fullRE_DMSL <- rep(NA, length(list_RE[[which(typeofs_of_RE == "double")[1]]]))
        }
        
        if(typeofs_of_RE[2] == "character"){
          x_RE_diagRE_DMSL <- rep(NA, length(list_RE[[which(typeofs_of_RE == "double")[1]]]))
        }
        
        if(typeofs_of_RE[3] == "character"){
          x_RE_fullRE_M <- rep(NA, length(list_RE[[which(typeofs_of_RE == "double")[1]]]))
        }
        
        if( (length(x_RE_fullRE_DMSL) != length(x_RE_diagRE_DMSL)) | (length(x_RE_fullRE_DMSL) != length(x_RE_fullRE_M)) ){
          warning(paste0(ct, ': the number of log-ratios is not consistent'))
          return(NULL)
        }
      }
      
    }
    
    ## put the coefficients in matrix form
    ## get the number of log-ratios, d-1
    dmin1 <- (names(table(sapply(list(model_fullRE_DMSL_list, model_diagRE_DMSL_list, model_fullRE_M_list), function(i)    as.numeric(try(length(python_like_select_name(i[[ct]]$par.fixed, 'beta'))/2))))))
    if(length(dmin1) == 1){
      ## there should only be one, shared, d-1
      dmin1 <- as.numeric(dmin1)
    }else{
      stop(paste0('Models do not agree on number of log-ratios. CT: ', ct))
    }
    
    x_RE_fullRE_DMSL <- matrix(x_RE_fullRE_DMSL, ncol=dmin1, byrow=F)
    x_RE_diagRE_DMSL <- matrix(x_RE_diagRE_DMSL, ncol=dmin1, byrow=F)
    x_RE_fullRE_M <- matrix(x_RE_fullRE_M, ncol=dmin1, byrow=F)
    
    bad_fullRE_DMSL=F
    bad_diagRE_DMSL=F
    bad_fullRE_M=F
    ## if something hasn't converged, set all the random coefficients to NA
    if((typeof(model_fullRE_DMSL_list[[ct]]) == "character") ){
      bad_fullRE_DMSL=T
    }else{
      if(try(!(model_fullRE_DMSL_list[[ct]]$pdHess))){
        bad_fullRE_DMSL=T
      }
    }
    if(bad_fullRE_DMSL){
      x_RE_fullRE_DMSL <- matrix(NA, nrow = nrow(x_RE_fullRE_DMSL), ncol=ncol(x_RE_fullRE_DMSL))
    }
    #----
    if((typeof(model_diagRE_DMSL_list[[ct]]) == "character") ){
      bad_diagRE_DMSL=T
    }else{
      if(try(!(model_diagRE_DMSL_list[[ct]]$pdHess))){
        bad_diagRE_DMSL=T
      }
    }
    if(bad_diagRE_DMSL){
      x_RE_diagRE_DMSL <- matrix(NA, nrow = nrow(x_RE_diagRE_DMSL), ncol=ncol(x_RE_diagRE_DMSL))
    }
    #-----
    if((typeof(model_fullRE_M_list[[ct]]) == "character") ){
      bad_fullRE_M=T
    }else{
      if(try(!(model_fullRE_M_list[[ct]]$pdHess))){
        bad_fullRE_M=T
      }
    }
    if(bad_fullRE_M){
      x_RE_fullRE_M <- matrix(NA, nrow = nrow(x_RE_fullRE_M), ncol=ncol(x_RE_fullRE_M))
    }

    ## for each patient using the x_RE_fullRE_DMSL intercepts, get the distance to the intercepts of the other two models
    dist_DMSLs <- sapply(1:nrow(x_RE_fullRE_DMSL), function(i){
      if(all(is.na(x_RE_fullRE_DMSL[i,])) | all(is.na(x_RE_diagRE_DMSL[i,]))){
        NA
      }else{
        dist(rbind(x_RE_fullRE_DMSL[i,], x_RE_diagRE_DMSL[i,]))
      }
    })
    dist_fullREs <- sapply(1:nrow(x_RE_fullRE_DMSL), function(i){
      if(all(is.na(x_RE_fullRE_DMSL[i,])) | all(is.na(x_RE_fullRE_M[i,]))){
        NA
      }else{
        dist(rbind(x_RE_fullRE_DMSL[i,], x_RE_fullRE_M[i,]))
      }
    })
    
    cbind.data.frame(melt(list(dist_DMSLs=dist_DMSLs, dist_fullREs=dist_fullREs)),
                     ct=ct)
  }
  ))
  
  # ## if something hasn't converged, remove the value
  # bad_ct_fullRE_DMSL <- c(enough_samples[sapply(enough_samples,
  #                                               function(ct) (typeof(model_fullRE_DMSL_list[[ct]])) == 'character')],
  #                         enough_samples[sapply(enough_samples,
  #                                               function(ct)  try(model_fullRE_DMSL_list[[ct]]$pdHess)) != "TRUE"])
  # bad_ct_diagRE_DMSL <- c(enough_samples[sapply(enough_samples,
  #                                               function(ct) (typeof(model_diagRE_DMSL_list[[ct]])) == 'character')],
  #                         enough_samples[sapply(enough_samples, function(ct)  try(model_diagRE_DMSL_list[[ct]]$pdHess)) != "TRUE"])
  # bad_ct_fullRE_M <- c(enough_samples[sapply(enough_samples,
  #                                            function(ct) (typeof(model_fullRE_M_list[[ct]])) == 'character')],
  #                      enough_samples[sapply(enough_samples, function(ct)  try(model_fullRE_M_list[[ct]]$pdHess)) != "TRUE"])
  # 
  # .x$fullRE_DMSL[(.x$ct %in% bad_ct_fullRE_DMSL)] = NA
  # .x$diagRE_DMSL[(.x$ct %in% bad_ct_diagRE_DMSL)] = NA
  # .x$fullRE_M[(.x$ct %in% bad_ct_fullRE_M)] = NA
  
  .x$ct2=renaming_pcawg[,2][match(.x$ct, renaming_pcawg[,1])]
  
  return(.x)
  
}


# comparison_randomintercepts_models_all <- comparison_randomintercepts_models(model_fullRE_DMSL_list = fullRE_DMSL,
#                                                                   model_diagRE_DMSL_list = diagRE_DMSL,
#                                                                   model_fullRE_M_list = fullRE_M)
comparison_randomintercepts_models_nonexo <- comparison_randomintercepts_models(model_fullRE_DMSL_list = fullRE_DMSL_nonexo,
                                                                     model_diagRE_DMSL_list = diagRE_DMSL_nonexo,
                                                                     model_fullRE_M_list = fullRE_M_nonexo)


head(comparison_randomintercepts_models_nonexo)

comparison_randomintercepts_models_nonexo[comparison_randomintercepts_models_nonexo$ct %in% (comparison_randomintercepts_models_nonexo %>% group_by(ct, L1) %>% summarise(good_dist=any(!is.na(value))) %>% group_by(ct) %>% summarise(good_dist2=sum(good_dist)) %>% filter(good_dist2 == 2) %>% select(ct) %>% unlist()),]

comparison_randomintercepts_models_nonexo <- as(comparison_randomintercepts_models_nonexo, 'data.frame')
comparison_randomintercepts_models_nonexo <- do.call(rbind.data.frame, comparison_randomintercepts_models_nonexo)
ldply (comparison_randomintercepts_models_nonexo, data.frame)
# comparison_randomintercepts_models_nonexo <- data.frame(apply(comparison_randomintercepts_models_nonexo, 2, function(i) i))
# comparison_randomintercepts_models_nonexo$value <- as.numeric(comparison_randomintercepts_models_nonexo$value)
ggplot((as.data.frame(comparison_randomintercepts_models_nonexo)),
       aes(x=(ct2),
           y=value))+geom_point()

