##------------------------------------------------------------------------##
add_metadata <- function(df, vector_from_DA_with_runname){
  df$params = sapply(vector_from_DA_with_runname, function(i) paste0(strsplit(i, '_')[[1]][1:6], collapse = '_'))
  df$model = sapply(vector_from_DA_with_runname, function(i) strsplit(i, '_')[[1]][1])
  df$n = sapply(vector_from_DA_with_runname, function(i) strsplit(i, '_')[[1]][2])
  df$nlambda = sapply(vector_from_DA_with_runname, function(i) strsplit(i, '_')[[1]][3])
  df$pi = as.numeric(sapply(vector_from_DA_with_runname, function(i) strsplit(i, '_')[[1]][6]))
  df$pi_softmax = sapply(df$pi, function(i) softmax(c(i, 0))[1])
  return(df)
}
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
modify_names_HiLDA <- function(j){
  sapply(j, function(i){
    i <- sapply(i, rename_datasets_fun)
    i <- paste0(gsub("HiLDA.*", "", i), gsub(".*dataset", "", i))
    return(i)
  })
}

get_inference_files <- function(folder_in, verbose=T, remove_HiLDAGlobal=F, runtime=F, return_filenames=F, HiLDA_return='all'){
  in_files <- unlist(sapply(grep('GenerationMixtureSimulation_|GenerationMixtureSimulationv4_|GenerationMixtureSimulationv7_|GenerationMixtureSimulationTwoCT', 
                     list.files(folder_in, full.names = T),
                     value = T),
                list.files, full.names = T))
  
  ## there are the results from inference and the time it took for these models to run - select only the results
  if(runtime){
    in_files <- grep('.time', in_files, value = T)
  }else{
    in_files <- grep('.RDS', in_files, value = T)
  }
  
  if(remove_HiLDAGlobal){ in_files <- in_files[grep("HiLDA_", in_files)]}
  
  
  if(verbose){
    cat(length(in_files), ' files found\n')
  }
  
  if(runtime){
    all_out_TMB <- lapply(in_files, function(j){
      if(grepl('.time', j) & grepl('TCSM', j)){
        read.table(j)
      }else{
        readRDS(j)
      }
    })
  }else{
    if(return_filenames){
      all_out_TMB <- in_files
    }else{
      if(HiLDA_return == 'all'){
        all_out_TMB <- lapply(in_files, readRDS)
      }else if (HiLDA_return == 'betas'){
        all_out_TMB <- sapply(in_files, function(i){
          x = readRDS(i)
          output_to_beta_matrix(x,  model='HiLDA')
        }, simplify = F)
      }else if (HiLDA_return == 'rhats'){
        all_out_TMB <- sapply(in_files, function(i){
          x = readRDS(i)
          try(python_like_select_rownames(x$BUGSoutput$summary, 'alpha|beta')[,'Rhat'])
        }, simplify = F)
        
      }else{
        stop('<HiLDA_return> not found\n')
      }
    }
  }
  
  names(all_out_TMB) <- paste0(gsub(".RDS", "", 
                                    gsub("multiple_GenerationMixtureSimulation", "", basename(in_files))))
  names(all_out_TMB) <-  gsub("diagREDM_NA_NA_NA_dataset", "", names(all_out_TMB))
  names(all_out_TMB) <-  gsub("TCSM_NA_NA_NA_dataset", "", names(all_out_TMB))
  names(all_out_TMB) <-  gsub("HilDA_NA_NA_NA_dataset", "", names(all_out_TMB))
  return(all_out_TMB)
}


# get_dataset_files <- function(folder_in){
#   in_files <- unlist(sapply(grep('GenerationMixtureSimulationv4|GenerationMixtureSimulationv7|GenerationMixtureSimulationTwoCT|GenerationMixtureSimulation_', 
#                                  list.files(folder_in, full.names = T),
#                                  value = T),
#                             list.files, full.names = T))
#   
#   ## there are the results from inference and the time it took for these models to run - select only the results
#   in_files <- grep('.RDS', in_files, value = T)
#   
#   all_out_TMB <- lapply(in_files, readRDS)
#   names(all_out_TMB) <- paste0(gsub(".RDS", "",  
#                                     gsub("multiple_GenerationMixtureSimulation", "", basename(in_files))))
#   names(all_out_TMB) <-  gsub("diagREDM_NA_NA_NA_dataset", "", names(all_out_TMB))
#   names(all_out_TMB) <-  gsub("TCSM_NA_NA_NA_dataset", "", names(all_out_TMB))
#   names(all_out_TMB) <-  gsub("HilDA_NA_NA_NA_dataset", "", names(all_out_TMB))
#   return(all_out_TMB)
# }
##------------------------------------------------------------------------##


##------------------------------------------------------------------------##
## Rename the datasets
rename_datasets_vec = c('TwoCT'= 'C1', 'v7'='C2', 'v4' = 'C3')
rename_datasets_fun <- function(i){
  for(k in 1:length(rename_datasets_vec)){
    i <- gsub(names(rename_datasets_vec)[k], rename_datasets_vec[k], x = i)
  }
  if(grepl('^_', i)){
    ## dataset from generation GenerationMixtureSimulation -- add the number of mutations as real
    ## they belong to Generation C3 (they are equivalent to v4, as there is patient-specific information)
    i = gsub('^_', 'C3_', i)
    isplit <- strsplit(i, '_')[[1]]
    isplit[3] <- 'Obs'
    i = paste0(isplit, collapse = '_')
  }
  i
}
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
hildaGlobal_from_beta <- function(jagsOutput){
  ## see if zero is not included in 95% credible interval of any signature (in which case there is DA)
  credible_intervals =   apply(jagsOutput$BUGSoutput$sims.list$beta, 2, function(j) quantile(j, probs=c(0.025, 0.975)))
  DA = !any(apply(credible_intervals, 2, function(j) (j[1] < 0) & (j[2] > 0)))
  if(DA){
    0
  }else{
    1
  }
}

hildaGlobal_from_alpha <- function(jagsOutput){
  ## see if zero is not included in 95% credible interval of any signature (in which case there is DA)
  credible_intervals =   lapply(1:dim(jagsOutput$BUGSoutput$sims.list$alpha)[3], function(k){
    apply(jagsOutput$BUGSoutput$sims.list$alpha[,,k], 2, function(j) quantile(j, probs=c(0.025, 0.975)))
  })
  DA = any(sapply(credible_intervals, function(j) all(j[,1] < j[2,1]) | all(j[,1] > j[2,2])))
  if(DA){
    0
  }else{
    1
  }
}

output_to_beta_matrix <- function(output, model){
  ## get the beta coefficients for each type of input data
  if(length(output) == 0){
    return(NA)
  }else{
    if(model =='TCSM'){
      return_obj <- output$effect
    }else if(model == 'TMB'){
      return_obj <- plot_betas(output, return_df = T, plot = F)
      return_obj <- dcast(return_obj, LogR~type_beta, value.var = 'Estimate')[,-1]
    }else if(model == 'HiLDA'){
      return_obj <- data.frame(output$BUGSoutput$summary) %>% slice(grep('alpha', row.names(.)))
      return_obj <- matrix((return_obj$mean), ncol=2, byrow=T)
      return_obj <- cbind(return_obj[,1], return_obj[,2]-return_obj[,1])
    }
    return(return_obj)
  }
}
##------------------------------------------------------------------------##
