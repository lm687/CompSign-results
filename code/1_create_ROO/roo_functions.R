##' to add:
##' - function to create the datasets from vcf, e.g. for creating the features
##'   objects
##' - add a slot in the object saying what are the name of the features (cosmic signatures, or features, etc.)
##'   and what are their names
##' - make function to plot the changes between conditions, including zeros as an extra facet, or something like that.

setClass("exposures_cancertype",
         representation = representation(
           cancer_type="character",
           type_classification = "character",
           number_categories = "numeric",
           id_categories = "character",
           active_signatures = "character", ## active signatures for this cancer type
           count_matrices_all = "list", ## for each of the categories
           count_matrices_active = "list", ## for each of the categories, only active signatures
           sample_names = "character",
           is_null_active = "logical",
           is_empty = "character",
           modification = "character"
         ))


#' Output a 2x2 matrix with TRUE, FALSE columns and rows indicating how many samples have, for trunk and branch,
#' a zero exposure
compute_activity_matrix_ROOSigs <- function(exposures_cancertype_obj, active_only, normalised){
  ## for each signature
  if(active_only){
    if(exposures_cancertype_obj@is_null_active){
      return_table <- NA
      # stop(paste0('No active signatures in cancer type ', exposures_cancertype_obj@cancer_type))
    }else{
      return_table <- list()
      for(sig in 1:ncol(exposures_cancertype_obj@count_matrices_active[[1]])){
        cbind_all_categories <- do.call('cbind', lapply(1:exposures_cancertype_obj@number_categories, function(j){
          exposures_cancertype_obj@count_matrices_active[[j]][,sig]}))
        cbind_all_categories <- cbind_all_categories > 0
        if(exposures_cancertype_obj@number_categories > 2 | exposures_cancertype_obj@number_categories == 1){
          stop('Cannot compute the matrix for more than two categories, or only one')
        }
        return_table[[sig]] <- table(assign(exposures_cancertype_obj@id_categories[1], factor(cbind_all_categories[,1], levels=c(T,F))),
                                     assign(exposures_cancertype_obj@id_categories[2], factor(cbind_all_categories[,2], levels=c(T,F))))
        if(normalised){
          return_table[[sig]] <- return_table[[sig]]/sum(return_table[[sig]])
        }
      }
      names(return_table) <- colnames(exposures_cancertype_obj@count_matrices_active[[1]])
    }
  }else{
    return_table <- list()
    for(sig in 1:ncol(exposures_cancertype_obj@count_matrices_all[[1]])){
      cbind_all_categories <- do.call('cbind', lapply(1:exposures_cancertype_obj@number_categories, function(j){
        exposures_cancertype_obj@count_matrices_all[[j]][,sig]}))
      cbind_all_categories <- cbind_all_categories > 0
      if(exposures_cancertype_obj@number_categories > 2 | exposures_cancertype_obj@number_categories == 1){
        stop('Cannot compute the matrix for more than two categories, or only one')
      }
      return_table[[sig]] <- table(assign(exposures_cancertype_obj@id_categories[1], factor(cbind_all_categories[,1], levels=c(T,F))),
                            assign(exposures_cancertype_obj@id_categories[2], factor(cbind_all_categories[,2], levels=c(T,F))))
      if(normalised){
        return_table[[sig]] <- return_table[[sig]]/sum(return_table[[sig]])
      }
    }
    names(return_table) <- colnames(exposures_cancertype_obj@count_matrices_all[[1]])
  }
  return_table
}

#' Push near-zero exposures to zeros
pushzeros_ROOSigs <- function(exposures_cancertype_obj, type_pushing, verbatim=FALSE, path_helper="",
                               slot_name=c('count_matrices_all', 'count_matrices_active')){
  source(paste0(path_helper, "functions_helper_DA.R"))

  ##' within each object, we have an 'active signature' and 'all signature' list,
  ##' and each list contains matrices for trunk and branch (or more categories)
  ## for each category, for each signature, independently, push to zero
  
  new_objects_sigs_per_CT <- exposures_cancertype_obj

  for(type_exposure in slot_name){
    
    ## looping over active, and all signatures
    ## check if there are active signatures
      
    ## loop over categories of exposures
    for(category_exposure in 1:length(slot(exposures_cancertype_obj, type_exposure))){
      if(exposures_cancertype_obj@is_null_active & type_exposure == 'count_matrices_active'){
        ## no active found
        new_objects_sigs_per_CT@count_matrices_active <- exposures_cancertype_obj@count_matrices_active
      }else{
        .x <- slot(exposures_cancertype_obj, type_exposure)[[category_exposure]]
        ## loop over signatures
          slot(new_objects_sigs_per_CT, type_exposure)[[category_exposure]] <-  apply(log(abs(.x)), 2, function(sig){
            if(type_pushing == 'FMM'){
              sig[is.nan(sig)] <- min(sig, na.rm = TRUE)
              fmm_res <- give_fmm_zeros(sig, maxiter=20)
              sig[fmm_res] <- 0
              if(verbatim){ cat('Exposures of', length(fmm_res), 'samples were set to 0\n') }
              sapply(sig, function(jj) if(jj == 0){0}else{exp(jj)})
            }else if(type_pushing == 'threshold'){
              sig[is.nan(sig)] <- min(sig, na.rm = TRUE)
              fmm_res <- which(sig < -10)
              sig[fmm_res] <- 0
              if(verbatim){ cat('Exposures of', length(fmm_res), 'samples were set to 0\n') }
              sapply(sig, function(jj) if(jj == 0){0}else{exp(jj)})
            }
          })
        }
      }
    }
    
    new_objects_sigs_per_CT@modification <- paste0('(', new_objects_sigs_per_CT@modification, '), zero removal with ', type_pushing, ' method')
    
    return(new_objects_sigs_per_CT)
}
  
  
#' Push near-zero exposures to zeros
normalise_anchor_ROOSigs <- function(exposures_cancertype_obj, combination_anchor_sigs){
    cat('To be implemented. Combination of anchor signatures can be any number of signatures')
}

#' Keep only samples whihc don't have any zero, for each signature, or for all combined signatures
keep_stricly_pos_ROOSigs <- function(exposures_cancertype_obj, types_exposure=c('count_matrices_all', 'count_matrices_active')){

  new_objects_sigs_per_CT <- exposures_cancertype_obj
  
  for(type_exposure in types_exposure){
    
    ## looping over active, and all signatures
    ## check if there are active signatures
    
    ## loop over categories of exposures
    ncats <- length(slot(exposures_cancertype_obj, type_exposure))
    if(length(exposures_cancertype_obj@sample_names) != nrow(slot(exposures_cancertype_obj, 'count_matrices_all')[[1]])){
      cat('Samples', exposures_cancertype_obj@sample_names[!exposures_cancertype_obj@sample_names %in% rownames(slot(exposures_cancertype_obj, 'count_matrices_all')[[1]])], ' are not found in the exposures matrices but are listed in the sample slot.\n')
      numzeros <- matrix(0, ncol=ncats, nrow=nrow(slot(exposures_cancertype_obj, 'count_matrices_all')[[1]]))
    }else{
      numzeros <- matrix(0, ncol=ncats, nrow=length(exposures_cancertype_obj@sample_names))
    }  
      
    
    for(category_exposure in 1:ncats){
      .x <- slot(exposures_cancertype_obj, type_exposure)[[category_exposure]]
      if(!is.null(.x)){
        numzeros[,category_exposure] <- rowSums(.x == 0)
      }
    }
    
    if(type_exposure == 'count_matrices_active' & exposures_cancertype_obj@is_null_active){
      ## no active found
      new_objects_sigs_per_CT@count_matrices_active <- exposures_cancertype_obj@count_matrices_active
    }else{
      ## now remove the rows which don't have only non-zero signatures
      for(category_exposure in 1:ncats){
        .x <- slot(exposures_cancertype_obj, type_exposure)[[category_exposure]]
        if(is.null(.x)){
          slot(new_objects_sigs_per_CT, type_exposure)[[category_exposure]] <- NULL
        }else{
          if(sum(rowSums(numzeros) == 0) > 0){
            slot(new_objects_sigs_per_CT, type_exposure)[[category_exposure]] <- .x[rowSums(numzeros) == 0,]
          }else{
            slot(new_objects_sigs_per_CT, type_exposure)[[category_exposure]] <- NA
          }
        }
      }
      ## change also the samples names
      if(sum(rowSums(numzeros) == 0) > 0){
        new_objects_sigs_per_CT@sample_names <- exposures_cancertype_obj@sample_names[rowSums(numzeros) == 0]
      }else{
        new_objects_sigs_per_CT@sample_names <- ""
      }
    }
  }
  new_objects_sigs_per_CT@modification <- paste0('(', new_objects_sigs_per_CT@modification, '), samples with any zero exposure removed')
  new_objects_sigs_per_CT
}

normalise_exposures_ROOSigs <- function(exposures_cancertype_obj){
  require(CompSign)
  
  new_objects_sigs_per_CT <- exposures_cancertype_obj
  
  for(type_exposure in c('count_matrices_all', 'count_matrices_active')){
    ncats <- length(slot(exposures_cancertype_obj, type_exposure))
    
    ## looping over active, and all signatures
    ## check if there are active signatures
    
    ## loop over categories of exposures

    if(exposures_cancertype_obj@is_null_active){
      ## no active found
      new_objects_sigs_per_CT@count_matrices_active <- exposures_cancertype_obj@count_matrices_active
    }else{
      if(ncats == 0){
        slot(new_objects_sigs_per_CT, type_exposure) <- list()
      }else{
        ## now remove the rows which don't have only non-zero signatures
        for(category_exposure in 1:ncats){
          .x <- slot(exposures_cancertype_obj, type_exposure)[[category_exposure]]
          if(is.null(.x) | all(is.na(.x))){
            slot(new_objects_sigs_per_CT, type_exposure)[[category_exposure]] <- NA
          }else{
            if(is.null(dim(.x))){
              ## one-dimensional
              slot(new_objects_sigs_per_CT, type_exposure)[[category_exposure]] <- .x/sum(.x)
            }else{
              slot(new_objects_sigs_per_CT, type_exposure)[[category_exposure]] <- normalise_rw(.x)
            }
          }
        }
      }
    }
  }
  new_objects_sigs_per_CT@modification <- paste0('(', new_objects_sigs_per_CT@modification, '), normalised exposures')
  new_objects_sigs_per_CT
}

#' infer the parameters of the dm from which signatures come
fit_Dirichlet_ROOSigs <- function(exposures_cancertype_obj, which_type_sigs=c('count_matrices_all', 'count_matrices_active')){
  #require(HMP)
  require(Compositional)
  require(DirichletReg)
  require(compositions)
  
  new_objects_sigs_per_CT <- exposures_cancertype_obj

  for(type_exposure in which_type_sigs){
    ncats <- length(slot(exposures_cancertype_obj, type_exposure))
    
    if(exposures_cancertype_obj@is_null_active){
      ## no active found
      return(0) ##new_objects_sigs_per_CT@count_matrices_active <- exposures_cancertype_obj@count_matrices_active
    }else{
      for(category_exposure in 1:ncats){
        ## compute the HMP (for all signatures) for each of the categories
        #fitting_hmp <- HMP::DM.MoM(round(objects_sigs_per_CT_only_nonzero[[1]]@count_matrices_active[[1]]))
        try({
          mat_data_estimate <- as(compositions::acomp(slot(exposures_cancertype_obj, type_exposure)[[category_exposure]]), 'matrix')
          ## potential initial params colSums(mat_data_estimate)/sum(mat_data_estimate)
          slot(new_objects_sigs_per_CT, paste0('fitted_distrib_', type_exposure))[[category_exposure]] <- Compositional::diri.est(mat_data_estimate, type = "mle")},
        )
      }
    
      ## which to plot? use function below
      #give_which_ternary()
      
    }
  }
  new_objects_sigs_per_CT
}
  
#slot='count_matrices_active'
fit_HMP_ROOSigs <- function(exposures_cancertype_obj, slot='count_matrices_active'){
  require(HMP)
  mat <- slot(exposures_cancertype_obj, slot)
  # lapply(mat, function(m) HMP::DM.MoM(t(m)))
  lapply(mat, HMP::DM.MoM)
}

#' Test for equality of two compositions
test_HMP_ROOSigs <- function(exposures_cancertype_obj, slot='count_matrices_active'){
  require(HMP)
  HMP::Xdc.sevsample(list(t(slot(exposures_cancertype_obj, slot)[[1]]),
                        t(slot(exposures_cancertype_obj, slot)[[2]])))
}

fit_DMM_ROOSigs <- function(exposures_cancertype_obj, slot_name='count_matrices_active', kmax=15, path_helper){
  
  res_dmm <- lapply(1:2, function(it){
    if(! (sum(sapply(slot(pushzeros_ROOSigs(exposures_cancertype_obj,
                                          path_helper = path_helper,
                                          type_pushing = 'threshold', slot=slot_name), slot_name), dim)) > 0) ){
      NA
    }else{
      lapply(1:kmax, function(k) DirichletMultinomial::dmn(slot(pushzeros_ROOSigs(exposures_cancertype_obj,
                                                                                   path_helper = path_helper,
                                                                                   type_pushing = 'threshold',
                                                                                   slot=slot_name), slot_name)[[it]],
                                                           k=k))
      }
  })
  # par(mfrow=c(2,5))
  # sapply(1:2, function(it) sapply(1:5, function(goodnessoffit_measure) plot(1:15, sapply(res_dmm[[it]], function(ft){ft@goodnessOfFit[goodnessoffit_measure]}), type='l')))
  lapply(1:2, function(it) t(sapply(res_dmm[[it]], function(ft){ft@goodnessOfFit})))
  
  # par(mfrow=c(1,2))
  # plot(res_dmm[[1]][[2]]@group[,1], type = 'l')
  # plot(res_dmm[[2]][[2]]@group[,1], type = 'l')
  # 
  # assigned_cluster <- sapply(1:2, function(it) apply(res_dmm[[it]][[2]]@group, 1, which.max))
  # 
  # table(assigned_cluster[,1], assigned_cluster[,2])
  # ##' most stays in the same cluster. Only two samples have have been assigned to to cluster 1 in the trunk and cluster 2
  # ##' and two other samples in cluster 2 in the trunk and cluster 1 in the branches (it needs not be reciprocal).
  # 
  # grid.arrange(createBarplot(normalise_rw(df_pushed[[1]][assigned_cluster[,1] == 1,]))+ggtitle('Cluster 1, trunk'),
  #              createBarplot(normalise_rw(df_pushed[[1]][assigned_cluster[,1] == 2,]))+ggtitle('Cluster 1, branch'),
  #              createBarplot(normalise_rw(df_pushed[[2]][assigned_cluster[,2] == 1,]))+ggtitle('Cluster 2, trunk'),
  #              createBarplot(normalise_rw(df_pushed[[2]][assigned_cluster[,2] == 2,]))+ggtitle('Cluster 2, branch'), ncol=2)
}

alr_ROOSigs <- function(exposures_cancertype_obj, baseline_sig='Signature.1', presence_zeros=FALSE, log_bool=TRUE, verbatim=TRUE){
  prev_object <- exposures_cancertype_obj
  
  if(verbatim){ cat('We are only looking at active signatures\n') }
  if(exposures_cancertype_obj@is_null_active| all(sapply(exposures_cancertype_obj@count_matrices_active, is.na))){
    ## no active found
    exposures_cancertype_obj@count_matrices_active <- exposures_cancertype_obj@count_matrices_active
  }else{
    if(!presence_zeros & any( sapply(exposures_cancertype_obj@count_matrices_active, function(x) sum(x == 0)) > 0)){
      stop('There are zeros when there should not be any. If you expect zeros to be present, change <presence_zeros> to TRUE.')
    }
    
    for(type_exposure in c('count_matrices_active')){
      ncats <- length(slot(exposures_cancertype_obj, type_exposure))
      
      ## looping over active, and all signatures
      ## check if there are active signatures
      
      ## loop over categories of exposures

      ## only active signatures
      for(cats in 1:prev_object@number_categories){
        if(all(is.na(slot(exposures_cancertype_obj, type_exposure)[[cats]])) | is.null(slot(exposures_cancertype_obj, type_exposure)[[cats]])){
          slot(exposures_cancertype_obj, type_exposure)[[cats]] <- NA
        }else{
          nonlogged_exposures <- (slot(exposures_cancertype_obj, type_exposure)[[cats]])
          if(!log_bool){
            logged_exposures <- nonlogged_exposures
          }else{
            logged_exposures <- log(nonlogged_exposures)
          }
          if(is.null(dim(logged_exposures))){
            ## only one row
            if(length(strsplit(baseline_sig, '[+]')[[1]]) > 1){
              ## it is the sum of multiples
              exposure_base <- nonlogged_exposures[strsplit(baseline_sig, '[+]')[[1]]]
              if(!log_bool){
                exposure_base <- sum(exposure_base)
              }else{
                exposure_base <- log(sum(exposure_base))
              }
            }else{
              exposure_base <- logged_exposures[baseline_sig]
            }
            zeros_baseline_bool <- (is.infinite(-exposure_base))
            if(zeros_baseline_bool){
              ## set all to NaN
              names_matrix <- names(slot(exposures_cancertype_obj, type_exposure)[[cats]])
              new_matrix <- rep(NaN, length(names_matrix))
              names(new_matrix) <- names_matrix
              slot(exposures_cancertype_obj, type_exposure)[[cats]] <- new_matrix
            }else{
              if(!log_bool){
                slot(exposures_cancertype_obj, type_exposure)[[cats]] <- logged_exposures/logged_exposures[baseline_sig]
              }else{
                if(verbatim) cat('Important! there used to be a bug here. ')
                ## code used to read
                ## slot(exposures_cancertype_obj, type_exposure)[[cats]] <- logged_exposures/logged_exposures[baseline_sig]
                ## when it should have been
                ## slot(exposures_cancertype_obj, type_exposure)[[cats]] <- logged_exposures-logged_exposures[baseline_sig]
                slot(exposures_cancertype_obj, type_exposure)[[cats]] <- logged_exposures-logged_exposures[baseline_sig]
              }
            }
          }else{
            ## more than one row
            ## baseline signature might be the sum of multiple signatures
            if(length(strsplit(baseline_sig, '[+]')[[1]]) > 1){
              ## it is the sum of multiple signatures
              exposure_base <- nonlogged_exposures[,strsplit(baseline_sig, '[+]')[[1]]]
              if(!log_bool){
                exposure_base <- (rowSums(exposure_base))
              }else{
                exposure_base <- log(rowSums(exposure_base))
              }
            }else{
              exposure_base <- logged_exposures[,baseline_sig]
            }
            if(!log_bool){
              slot(exposures_cancertype_obj, type_exposure)[[cats]] <- sweep(logged_exposures, 1, exposure_base, '/')
            }else{
              slot(exposures_cancertype_obj, type_exposure)[[cats]] <- sweep(logged_exposures, 1, exposure_base, '-')
            }
            ## zeros_baseline
            zeros_baseline <- which(is.infinite(-exposure_base))
            # slot(exposures_cancertype_obj, type_exposure)[[cats]][zeros_baseline,][1:4,1:4]
            # slot(prev_object, type_exposure)[[cats]][zeros_baseline,][1:4,1:4]
            
            ## set all ALRs to NaN if the baseline signature is zero
            if(length(zeros_baseline) > 0){
              slot(exposures_cancertype_obj, type_exposure)[[cats]][zeros_baseline,] <- NaN
            }
          }
        }
        ## set the zeros to NA
        if(is.null(dim(logged_exposures))){
          ## only one row
          slot(exposures_cancertype_obj, type_exposure)[[cats]] <- sapply(slot(exposures_cancertype_obj, type_exposure)[[cats]], function(x){x[is.infinite(-x)] <-NA; x })
        }else{
          slot(exposures_cancertype_obj, type_exposure)[[cats]] <- apply(slot(exposures_cancertype_obj, type_exposure)[[cats]], 2, function(x){x[is.infinite(-x)] <-NA; x })
        }
      }
    }
  }
  exposures_cancertype_obj@modification = paste0("ALR by ", baseline_sig)
  exposures_cancertype_obj
}


alr_subtract_ROOSigs <- function(exposures_cancertype_obj, baseline_sig='Signature.1', presence_zeros=FALSE, verbatim=TRUE){
  prev_object <- exposures_cancertype_obj
  if(exposures_cancertype_obj@number_categories != 2){
    stop('Number of categories different from 2 has not been implemented yet')
  }else{
    ## two categories
    ## alrs
    alr_res <- alr_ROOSigs(exposures_cancertype_obj, baseline_sig, presence_zeros, verbatim=verbatim)
    ## now, compute the difference
    prev_object@count_matrices_active <- list( alr_res@count_matrices_active[[2]] - alr_res@count_matrices_active[[1]] )
  }
  prev_object@id_categories <- paste0(exposures_cancertype_obj@id_categories[2], '-', exposures_cancertype_obj@id_categories[1])
  prev_object
}


# exposures_cancertype_obj = objects_sigs_per_CT_only_nonzero[[19]]
# slot_name="count_matrices_active"
give_perturbation_ROOSigs <- function(exposures_cancertype_obj, slot_name="count_matrices_active"){
  prev_object <- exposures_cancertype_obj
  if(exposures_cancertype_obj@number_categories != 2){
    stop('Number of categories different from 2 has not been implemented yet')
  }else{
    ## two categories
    ## we will return a matrix with as many columns as signatures, and as many rows as rows originally (e.g. one each sample)
    ## but a single matrix, instead of two
    ## this relies on zero having been removed from the same samples, i.e. if a sample has zeros in
    ## the trunk but not in the branch, the sample will be removed in both cases
    
    if(is.null(slot(exposures_cancertype_obj, slot_name)[[1]])){
      slot(prev_object, slot_name) = list()
    }else{
      if(length(slot(exposures_cancertype_obj, slot_name)[[1]]) == 1){
        if(is.na(slot(exposures_cancertype_obj, slot_name)[[1]])){
          slot(prev_object, slot_name) = list()
        }
      }else{
        if(is.null(nrow(slot(exposures_cancertype_obj, slot_name)[[1]]))){
          ## only one row
          perturbation = normalise_rw( normalise_rw(t(matrix(slot(exposures_cancertype_obj, slot_name)[[2]])))/    normalise_rw(t(matrix(slot(exposures_cancertype_obj, slot_name)[[1]]))) )
          colnames(perturbation) = names(slot(exposures_cancertype_obj, slot_name)[[1]])
          rownames(perturbation) = exposures_cancertype_obj@sample_names
        }else{
          perturbation = normalise_rw( normalise_rw(slot(exposures_cancertype_obj, slot_name)[[2]])/    normalise_rw(slot(exposures_cancertype_obj, slot_name)[[1]]) )
        }
        ## correct
        # normalise_rw(slot(exposures_cancertype_obj, slot_name)[[2]])[1,]
        # normalise_rw(t(as.matrix(normalise_rw(slot(exposures_cancertype_obj, slot_name)[[1]])[1,]*perturbation[1,])))
        
        slot(prev_object, slot_name) = list(perturbation)
      }
    }
  }
  prev_object@id_categories <- "perturbation 1->2"
  prev_object
}

give_perturbation_coefficient_ROOSigs <- function(exposures_cancertype_obj, slot_name="count_matrices_active",
                                                  norm="L2", verbose=TRUE){
  if(verbose) message('This function has been modified as there was a calculation error in the previous version')
  if(length(slot(exposures_cancertype_obj, slot_name)) == 0){
    slot(exposures_cancertype_obj, slot_name) = list()
  }else{
    perturbation = slot(exposures_cancertype_obj, slot_name)[[1]] ## only one element in the list
    if(norm == "L2"){
      coefficients_perturbation = apply(perturbation, 1, function(i) sum(i)/sum(i**2)/length(perturbation))
      slot(exposures_cancertype_obj, slot_name) = list(sweep(perturbation, 1, coefficients_perturbation, '*'))
    }else{
      stop('No norm but L2 has been implemented.')
    }
  }
  exposures_cancertype_obj
}

## debug
# alr_res@count_matrices_active[[1]][1:4,1:4]
# exposures_cancertype_obj@count_matrices_active[[1]][1:4,1:4]

#' Plot boxplots
plot_boxplots_ROOSigs <- function(exposures_cancertype_obj, slot="count_matrices_active", test_equality=TRUE, relax_criteria=FALSE){
  new_objects_sigs_per_CT <- exposures_cancertype_obj
  
  num_cats <- length(slot(exposures_cancertype_obj, slot))
  first_criterion <- exposures_cancertype_obj@is_null_active
  if(relax_criteria){
    first_criterion <- FALSE
  }
  if( !(first_criterion | all(sapply(slot(exposures_cancertype_obj, slot), is.na))) ){

    ## test for equality to zero
    if(test_equality){
      ## if the data are the subtraction
      if( !is.null(dim(slot(exposures_cancertype_obj, slot)[[1]])) ){
        colnames(slot(exposures_cancertype_obj, slot)[[1]]) <- gsub('Signature[.]', '', colnames(slot(exposures_cancertype_obj, slot)[[1]]))
        
        cols=apply(slot(exposures_cancertype_obj, slot)[[1]], 2, function(x) try(t.test(x)$p.value < 0.05))
        cols = sapply(cols, function(cl) if(class(cl[[1]]) == "try-error"){NA}else{cl[[1]]})
        
      }else{
        cols <- rep(FALSE, length(slot(exposures_cancertype_obj, slot)[[1]]))
      }
      # cols <- factor(cols, levels=c(NA, FALSE, TRUE))
    }else{
      cols <- rep(1, ncol(slot(exposures_cancertype_obj, slot)[[1]]))
    }

    par(mfrow=c(1,num_cats))
    for(num_cats_it in 1:num_cats){
      boxplot(slot(exposures_cancertype_obj, slot)[[num_cats_it]],
              main=paste0(exposures_cancertype_obj@cancer_type, ': ', exposures_cancertype_obj@id_categories[num_cats_it]),
              col=c('white', 'red')[cols])
      abline(h = 0, lty='dashed')
    }
  }

  # for(type_exposure in c('count_matrices_all', 'count_matrices_active')){
  # }
}

#' Similar to plot_boxplots_ROOSigs, but here giving the p-values or statistics from the t-test as a table
alr_subtract_test_ROOSigs <- function(exposures_cancertype_obj){
  
  ## this object should only contain one matrix in count_matrices_active, as it should be the subtraction
  stopifnot(length(exposures_cancertype_obj@count_matrices_active) == 1)

  num_cats <- length(exposures_cancertype_obj@count_matrices_active)
  if( !(exposures_cancertype_obj@is_null_active | all(sapply(exposures_cancertype_obj@count_matrices_active, is.na))) ){
    
    if( !is.null(dim(exposures_cancertype_obj@count_matrices_active[[1]])) ){
      cols = apply(exposures_cancertype_obj@count_matrices_active[[1]], 2, function(x) list(try(t.test(x)$p.value)))
      cols = sapply(cols, function(cl) if(class(cl[[1]]) == "try-error"){NA}else{cl[[1]]})
      list(df_pvals=data.frame(Signature=names(cols), pval=cbind(cols), CT=exposures_cancertype_obj@cancer_type),
           mean_difference=colMeans(exposures_cancertype_obj@count_matrices_active[[1]], na.rm = TRUE))
    }else{
      list(df_pvalsdata.frame(Signature=NA, pval=NA, CT=exposures_cancertype_obj@cancer_type),
           mean_difference=NA)
      }
  }
}


#' Plot, for a given set of samples, the exposures of each of the signatures in their categories
plot_ternary_categories_ROOSigs <- function(exposures_cancertype_obj){
  require(compositions)
  stop('Under constructions')
  N_MI <- 1e3 ## number of markov instances
  fits_categories_df <- do.call('rbind', lapply(lapply(fits_categories, function(k) k$b), function(fit_obk){
    as(compositions::acomp(MCMCpack::rdirichlet(N_MI, fit_obk)), 'matrix')
  }))
  
  # plot(DR_data(MCMCpack::rdirichlet(1e5, fitting_dir$param)), a2d = list(colored = TRUE, c.grid = FALSE, col.scheme = c("entropy")))
  
  ## comparing the two categories of exposures
  df <- cbind.data.frame(dat=fits_categories_df,
                         col=rep(1:length(fits_categories), each=N_MI))
  plot(compositions::acomp(df[,-ncol(df)]), col=df[,ncol(df)], pch=19, cex=0.2)
}


#' 
#' Similar to plot_boxplots_ROOSigs, but here giving the p-values or statistics from the t-test as a table
multinomial_logistic_test_ROOSigs <- function(exposures_cancertype_obj){
  exposures_cancertype_obj@count_matrices_active
  
  .x <- lapply(1:ncol(exposures_cancertype_obj@count_matrices_active[[1]]), function(sig){
    sig_name <- colnames(exposures_cancertype_obj@count_matrices_active[[1]])[sig]
    rbind_exposures <- cbind.data.frame(exposures=c(exposures_cancertype_obj@count_matrices_active[[1]][,sig],
                             exposures_cancertype_obj@count_matrices_active[[2]][,sig]),
                             group=c(rep(1, dim(exposures_cancertype_obj@count_matrices_active[[1]])[1]),
                                     rep(2, dim(exposures_cancertype_obj@count_matrices_active[[2]])[1])
                             ))
    rbind_exposures$exposures <- round(rbind_exposures$exposures)
    
    ## flatten and do the same for each signature
    .x2 <- do.call('rbind', lapply(1:nrow(rbind_exposures), function(rw){if(rbind_exposures[rw,'exposures'] > 0 ){
      cbind.data.frame(Sig=rep(sig_name, rbind_exposures[rw,'exposures']), group=rep(rbind_exposures[rw,'group'], each=rbind_exposures[rw,'exposures']),
                       id=rep(rw, rbind_exposures[rw,'exposures']), stringsAsFactors = FALSE)}}))
    print(head(.x2))
    .x2
  })
  .x <- do.call('rbind', .x)
  
  ## I believe MGLM is the one we want
  
  ## example data
  exposures_cancertype_obj=(objects_sigs_per_CT_pushedzeros[[1]])
  
  .x_for_mglm <- rbind.data.frame(cbind.data.frame(exposures_cancertype_obj@count_matrices_active[[1]],
                                        id=1:nrow(exposures_cancertype_obj@count_matrices_active[[1]]),
                                        group=1),
                       cbind.data.frame(exposures_cancertype_obj@count_matrices_active[[2]],
                                        id=1:nrow(exposures_cancertype_obj@count_matrices_active[[2]]),
                                        group=2))
  
  require(nnet)
  test <- multinom(Sig ~ group, data = .x)
  z <- summary(test)$coefficients/summary(test)$standard.errors
  # 2-tailed Wald z tests to test significance of coefficients
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  p
  
  ## https://stats.stackexchange.com/questions/63222/getting-p-values-for-multinom-in-r-nnet-package

  ## using MGLM
  require(MGLM)
  gdm.reg <- MGLMreg(as(.x_for_mglm[,1:(ncol(.x_for_mglm)-2)], 'matrix')~as(.x_for_mglm$group, 'matrix'), dist="GDM", LRT=FALSE)
  gdm.reg2 <- MGLMreg(round(as(.x_for_mglm[,1:(ncol(.x_for_mglm)-2)], 'matrix'))~as(.x_for_mglm$group, 'matrix'), dist="GDM", LRT=FALSE)
  
  ## this does give a coefficient for all signatures but one:
  ## alpha_Signature.4 alpha_Signature.5 alpha_Signature.13 alpha_Signature.2 beta_Signature.4 beta_Signature.5 beta_Signature.13 beta_Signature.2
  ## but the p-value is for the whole dataset
  

}

# exposures_cancertype_obj <- objects_sigs_per_CT_pushedzeros[[1]]
# pre_path="../../"
createbarplot_ROOSigs <- function(exposures_cancertype_obj, pre_path="../../", slot="count_matrices_active", ...){
  source(paste0(pre_path, "functions/meretricious/pretty_plots/prettySignatures.R"))
  require(gridExtra)
  
  if(slot == "count_matrices_active"){
    cat('Using only active signatures\n')
  }else if(slot == "count_matrices_all"){
    cat('Using all signatures\n')
  }else{
    stop('Signature slot not found')
  }
      
  mat <- slot(object = exposures_cancertype_obj, name = slot)
  if(length(mat) != 2){
    stop('Not implemented yet: there are more than two categories')
  }
  
  grid.arrange(createBarplot(normalise_rw(mat[[1]]), ...)+ggtitle(exposures_cancertype_obj@id_categories[1])+theme(legend.position = "bottom"),
               createBarplot(normalise_rw(mat[[2]]), ...)+ggtitle(exposures_cancertype_obj@id_categories[2])+theme(legend.position = "bottom"), ncol=2)
  
}

# exposures_cancertype_obj <- objects_sigs_per_CT_pushedzeros[[1]]
# pre_path="../../"
# slot="count_matrices_active"
#' Colour by bool of whether the signature is an anchor signature (1 or 5) or not
createbarplot_colouranchor_ROOSigs <- function(exposures_cancertype_obj, pre_path="../../", slot="count_matrices_active", angle_rotation_axis = 0){
  source(paste0(pre_path, "functions/meretricious/pretty_plots/prettySignatures.R"))
  require(gridExtra)
  
  if(slot == "count_matrices_active"){
    cat('Using only active signatures\n')
  }else if(slot == "count_matrices_all"){
    cat('Using all signatures\n')
  }else{
    stop('Signature slot not found')
  }
  
  mat0 <- slot(object = exposures_cancertype_obj, name = slot)
  if(length(mat0) != 2){
    stop('Not implemented yet: there are more than two categories')
  }
  
  mat <- melt(normalise_rw(mat0[[1]]))
  mat[,'bool'] <- mat$Var2 %in% c('Signature.1', 'Signature.5')
  plt1 <- ggplot(mat, aes(x=Var1, y=value, fill=bool))+geom_bar(stat='identity')+ggtitle(exposures_cancertype_obj@id_categories[1])+theme(legend.position = "bottom")+theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))
  mat <- melt(normalise_rw(mat0[[2]]))
  mat[,'bool'] <- mat$Var2 %in% c('Signature.1', 'Signature.5')
  plt2 <- ggplot(mat, aes(x=Var1, y=value, fill=bool))+geom_bar(stat='identity')+ggtitle(exposures_cancertype_obj@id_categories[2])+theme(legend.position = "bottom")+theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))
  
  grid.arrange(plt1, plt2, ncol=2)
  
}

##' Select only n-k samples. The matrix of active counts and all counts is going to be modified leaving out the same samples.
##' Which samples to remove can also be specified
leave_k_out_ROOSigs <- function(exposures_cancertype_obj, k=NULL, which_rm=NULL, which_slot="count_matrices_active"){
  
  for(slot_it in which_slot){
    if(!all(sapply(slot(exposures_cancertype_obj, slot_it), function(i) is.null(nrow(i))))){
      # if(all(sapply(slot(exposures_cancertype_obj, which_slot), is.na)){
      # }
    # }
    
    if(!is.null(which_rm)){
      which_keep <- which_rm[! (nrow(exposures_cancertype_obj@count_matrices_all) %in% which_rm )]
    }else{
      if(is.null(k)){
        cat('Selecting 80% of samples\n')
        which_keep <- sample(1:nrow(slot(exposures_cancertype_obj, slot_it)[[1]]), round(0.8*nrow(slot(exposures_cancertype_obj, slot_it)[[1]])))
      }else{
        which_keep <- sample(1:nrow(slot(exposures_cancertype_obj, slot_it)[[1]]), k)
      }
      which_keep <- sort(which_keep)
    }
      slot(exposures_cancertype_obj, slot_it) <- lapply(slot(exposures_cancertype_obj, slot_it), function(i) i[which_keep,])
    }
  }
  exposures_cancertype_obj@sample_names <- exposures_cancertype_obj@sample_names[which_keep]
  exposures_cancertype_obj@modification <- paste0('(', exposures_cancertype_obj@modification, '), downsampling')
  exposures_cancertype_obj
}

new_object <- function(){
  
  cat('TO CONTINUE')
  
  ## let's create an object that h
  
  new(Class = "exposures_cancertype",
    cancer_type="simulation",
    type_classification = "simulation",
    number_categories = 2,
    id_categories = c('SimTrunk', 'SimBranch'),
    active_signatures = "character", ## active signatures for this cancer type
    count_matrices_all = "list", ## for each of the categories
    count_matrices_active = "list", ## for each of the categories, only active signatures
    sample_names = "character",
    modification = "",
    is_null_active = FALSE,
    is_empty = "Non-empty"
)
  
}

rbind_matrices <- function(exposures_cancertype_obj, which_slot="count_matrices_active"){
  do.call('rbind', lapply(1:length(slot(exposures_cancertype_obj, which_slot)), function(i){
    x <- slot(exposures_cancertype_obj, which_slot)[[i]]
    rownames(x) <- paste0(rownames(x), '_', exposures_cancertype_obj@id_categories[i])
    x
  }))
}

# exposures_cancertype_obj = objects_sigs_per_CT_features[[1]]
# slot_name = "count_matrices_all"
context_features_to_features <- function(exposures_cancertype_obj, slot_name){
  .x = slot(exposures_cancertype_obj, slot_name)
  new_slot = lapply(1:length(.x), function(i){ ## to each of the lists
    if(is.null(.x[[i]])){
      NA
    }else{
      splt_col = sapply(sapply(colnames(.x[[i]]), function(j) strsplit(j, "\\[|\\]")), function(k) k[2])
      if(nrow(.x[[i]]) == 1){
        ## only one row
        .y = do.call('cbind', lapply(sort(unique(splt_col)), function(l) sum(.x[[i]][,splt_col == l])))
        rownames(.y)= rownames(.x[[i]])
      }else{
        ## multiple rows
        .y = do.call('cbind', lapply(sort(unique(splt_col)), function(l) rowSums(.x[[i]][,splt_col == l])))
      }
      colnames(.y) = sort(unique(splt_col))
      .y
    }
  })
  slot(exposures_cancertype_obj, slot_name) = new_slot
  exposures_cancertype_obj
}

aitchison_perturbation_test = function(object, slot_name){
  .x = (give_perturbation_ROOSigs(object, slot_name=slot_name))
  ## is the centre of the group the identity perturbation? i.e. is this below the zero vector?
  if(length(slot(.x, slot_name)) == 0){
    return(NA)
  }else{
    .x = as(compositions::alr(slot(.x, slot_name)[[1]]), 'matrix')
    .res = try(expr = {Hotelling::hotelling.test(.x, matrix(0, nrow=nrow(.x), ncol=ncol(.x)))})
    if (class(.res) == "try-error") {
      .res = NA
    }
    .res
  }
}

aitchison_perturbation_test_alt = function(object, slot_name){
  .x = give_perturbation_ROOSigs_alt(exposures_cancertype_obj = object, slot_name=slot_name)
  ## is the centre of the group the identity perturbation? i.e. is this below the zero vector?
  if(length(slot(.x, slot_name)) == 0){
    return(NA) ## if there was a slot of active signatures or all signatures, but were removed in the <keep strictly positive> step
  }
  
  nrow_obj = nrow(slot(.x, slot_name)[[1]])
  ncol_obj = ncol(slot(.x, slot_name)[[1]])
  ## Compute ALR
  if(nrow_obj==1){
    # .x = (log(slot(.x, slot_name)[[1]]) - log(slot(.x, slot_name)[[1]][,ncol_obj]))[-ncol_obj]
    .res = NA
  }else{
    logged = log(t(apply(slot(.x, slot_name)[[1]], 1, as.numeric)))
    .x = sweep(logged, 1, logged[,ncol_obj], '-')[,-ncol_obj]
    .res = try(expr = {Hotelling::hotelling.test(.x, matrix(0, nrow=nrow(.x), ncol=ncol(.x)))})
  }
  if (class(.res) == "try-error") {
    .res = NA
  }
  .res
}

give_perturbation_ROOSigs_alt = function(exposures_cancertype_obj, slot_name="count_matrices_active"){
  prev_object <- exposures_cancertype_obj
  if(exposures_cancertype_obj@number_categories != 2){
    stop('Number of categories different from 2 has not been implemented yet')
  }else{
    ## two categories
    ## we will return a matrix with as many columns as signatures, and as many rows as rows originally (e.g. one each sample)
    ## but a single matrix, instead of two
    ## this relies on zero having been removed from the same samples, i.e. if a sample has zeros in
    ## the trunk but not in the branch, the sample will be removed in both cases
    
    if(is.null(slot(exposures_cancertype_obj, slot_name)[[1]])){
      slot(prev_object, slot_name) = list()
    }else{
      if(length(slot(exposures_cancertype_obj, slot_name)[[1]]) == 1){
        if(is.na(slot(exposures_cancertype_obj, slot_name)[[1]])){
          slot(prev_object, slot_name) = list()
        }
      }else{
        if(is.null(nrow(slot(exposures_cancertype_obj, slot_name)[[1]])) | (nrow(slot(exposures_cancertype_obj, slot_name)[[1]]) ==1) ){
          ## only one row
          perturbation = normalise_rw( normalise_rw(t(matrix(apply(slot(exposures_cancertype_obj, slot_name)[[2]], 1, as.numeric)))) /    normalise_rw(t(matrix(apply(slot(exposures_cancertype_obj, slot_name)[[1]], 1, as.numeric)))) )
          colnames(perturbation) = names(slot(exposures_cancertype_obj, slot_name)[[1]])
          rownames(perturbation) = exposures_cancertype_obj@sample_names
        }else{
          perturbation = normalise_rw( normalise_rw(t(apply(slot(exposures_cancertype_obj, slot_name)[[2]], 1, as.numeric))) /    normalise_rw(t(apply(slot(exposures_cancertype_obj, slot_name)[[1]], 1, as.numeric))) )
        }
        ## correct
        # normalise_rw(slot(exposures_cancertype_obj, slot_name)[[2]])[1,]
        # normalise_rw(t(as.matrix(normalise_rw(slot(exposures_cancertype_obj, slot_name)[[1]])[1,]*perturbation[1,])))
        
        slot(prev_object, slot_name) = list(perturbation)
      }
    }
  }
  prev_object@id_categories <- "perturbation 1->2"
  prev_object
}

give_perturbation_ROOSigs_alt_v2 = function(exposures_cancertype_obj, slot_name="count_matrices_active"){
  ## without normalising the perturbations, as otherwise the system is singular
  prev_object <- exposures_cancertype_obj
  if(exposures_cancertype_obj@number_categories != 2){
    stop('Number of categories different from 2 has not been implemented yet')
  }else{
    ## two categories
    ## we will return a matrix with as many columns as signatures, and as many rows as rows originally (e.g. one each sample)
    ## but a single matrix, instead of two
    ## this relies on zero having been removed from the same samples, i.e. if a sample has zeros in
    ## the trunk but not in the branch, the sample will be removed in both cases
    
    if(is.null(slot(exposures_cancertype_obj, slot_name)[[1]])){
      slot(prev_object, slot_name) = list()
    }else{
      if(length(slot(exposures_cancertype_obj, slot_name)[[1]]) == 1){
        if(is.na(slot(exposures_cancertype_obj, slot_name)[[1]])){
          slot(prev_object, slot_name) = list()
        }
      }else{
        if(is.null(nrow(slot(exposures_cancertype_obj, slot_name)[[1]])) | (nrow(slot(exposures_cancertype_obj, slot_name)[[1]]) ==1) ){
          ## only one row
          perturbation = ( normalise_rw(t(matrix(apply(slot(exposures_cancertype_obj, slot_name)[[2]], 1, as.numeric)))) /    normalise_rw(t(matrix(apply(slot(exposures_cancertype_obj, slot_name)[[1]], 1, as.numeric)))) )
          colnames(perturbation) = names(slot(exposures_cancertype_obj, slot_name)[[1]])
          rownames(perturbation) = exposures_cancertype_obj@sample_names
        }else{
          perturbation = ( normalise_rw(t(apply(slot(exposures_cancertype_obj, slot_name)[[2]], 1, as.numeric))) /    normalise_rw(t(apply(slot(exposures_cancertype_obj, slot_name)[[1]], 1, as.numeric))) )
        }
        ## correct
        # normalise_rw(slot(exposures_cancertype_obj, slot_name)[[2]])[1,]
        # normalise_rw(t(as.matrix(normalise_rw(slot(exposures_cancertype_obj, slot_name)[[1]])[1,]*perturbation[1,])))
        
        slot(prev_object, slot_name) = list(perturbation)
      }
    }
  }
  prev_object@id_categories <- "perturbation 1->2"
  prev_object
}



bootstrap_for_statistic = function(object_roo, verbose){
  .object_roo_fun = object_roo
  if(verbose){
    'This function was wrong until 20210920. We were adding up the same first matrix twice'
  }
  ## previous; wrong
  # totals = .object_roo_fun@count_matrices_all[[1]] + .object_roo_fun@count_matrices_all[[1]] ## bootstrap
  ## current version, correct
  totals = .object_roo_fun@count_matrices_all[[1]] + .object_roo_fun@count_matrices_all[[2]] ## bootstrap
  ## sample in two groups
  .x_boot = lapply(1:nrow(totals), function(i){
    multinom1 = t(rmultinom(n=1, size = sum(.object_roo_fun@count_matrices_all[[1]][i,]), prob = totals[i,]/sum(totals[i,])))
    multinom2 = totals[i,]-multinom1
    list(multinom1, multinom2)
  })
  
  .object_roo_fun@count_matrices_all = list(do.call('rbind', lapply(.x_boot, function(i) i[[1]])),
                                            do.call('rbind', lapply(.x_boot, function(i) i[[2]])))
  rownames(.object_roo_fun@count_matrices_all[[1]]) = rownames(.object_roo_fun@count_matrices_all[[2]]) = rownames(object_roo@count_matrices_all[[1]])
  
  give_total_perturbation(.object_roo_fun, add_pseudocounts = TRUE)
  
}

permutation_test_fun_wrapper =  function(objects_sigs_per_CT_features_arg, nbootstraps, verbose=T){
  .x =  permutation_test_fun(objects_sigs_per_CT_features_arg, nbootstraps, verbose)
  return(sum(abs(.x[[1]]) > abs(.x[[2]]))/nbootstraps)
}

permutation_test_fun = function(objects_sigs_per_CT_features_arg, nbootstraps, verbose){
  null_statistics = replicate(n = nbootstraps, expr = bootstrap_for_statistic(objects_sigs_per_CT_features_arg, verbose))
  null_statistics = null_statistics[!is.nan(null_statistics)]
  observed_statistic = give_total_perturbation(objects_sigs_per_CT_features_arg)
  list(null_statistics, observed_statistic)
}

give_total_perturbation = function(object_roo, add_pseudocounts=TRUE, slot_name="count_matrices_all"){
  warning('This function was wrong until 20220208, as it assumed that the zero perturbation was 1/6, i.e. hardcoded for d=6')
  if(add_pseudocounts){
    slot(object_roo, slot_name)[[1]] =     slot(object_roo, slot_name)[[1]] + 1
    slot(object_roo, slot_name)[[2]] =     slot(object_roo, slot_name)[[2]] + 1
  }
  d = ncol(slot(object_roo, slot_name)[[1]])
  # sum(apply(slot(give_perturbation_ROOSigs(object_roo, slot_name = slot_name), name = slot_name)[[1]],
  # 1, function(i) sqrt(sum((i-1/6)**2))))/nrow(slot(object_roo, slot_name)[[1]])
  mean(apply(slot(give_perturbation_ROOSigs(object_roo, slot_name = slot_name), name = slot_name)[[1]],
             1, function(i) sqrt(sum((i-rep(1/d, d))**2))))
}

aitchison_perturbation_test_alt_v2 <- function(object, slot_name, addone=T){
  if(addone)  slot(object, name = slot_name) <- lapply(slot(object, name = slot_name), function(i) i +1) ## to avoid Inf
  .x = give_perturbation_ROOSigs_alt_v2(exposures_cancertype_obj = object, slot_name=slot_name)
  ## is the centre of the group the identity perturbation? i.e. is this below the zero vector?
  if(length(slot(.x, slot_name)) == 0){
    return(NA) ## if there was a slot of active signatures or all signatures, but were removed in the <keep strictly positive> step
  }
  
  nrow_obj = nrow(slot(.x, slot_name)[[1]])
  ncol_obj = ncol(slot(.x, slot_name)[[1]])
  ## Compute ALR
  if(nrow_obj==1){
    # .x = (log(slot(.x, slot_name)[[1]]) - log(slot(.x, slot_name)[[1]][,ncol_obj]))[-ncol_obj]
    .res = NA
  }else{
    .res = try(expr = {Hotelling::hotelling.test(slot(.x, name = slot_name)[[1]],
                                                 matrix(0, nrow=nrow_obj, ncol=ncol_obj))})
  }
  if (class(.res) == "try-error") {
    .res = NA
  }
  .res
}

iterative_chisqrt <- function(S1, S2){
  warning('This function was incorrect until 20220208. It was hardcoded for 96 categories')
  ## copied from forward_variable_selection.R
  ## where S1 and S2 are datasets, in which the columns are the categories
  ## 1. compute 96 independent chi-sqrt tests, for each of the categories
  ncats <- ncol(S1)
  unordered_pvals = lapply(1:ncats, function(mut_cat){
    chisq.test(rbind(c(sum(S1[,mut_cat]), sum(S1) - sum(S1[,mut_cat])),
                     c(sum(S2[,mut_cat]), sum(S2) - sum(S2[,mut_cat]))))
  })
  unordered_pvals = sapply(unordered_pvals, function(res_chisqrt) res_chisqrt$p.value)
  
  ## 2. rank by v-value (lowest to highest)
  rank = order(unordered_pvals, decreasing = FALSE)
  
  ordered_pvals = rep(NA, length(unordered_pvals))
  ordered_pvals[1] = min(unordered_pvals)
  
  ## 3. for mi, from lowest to highest, re-calculate p-value
  S1 = S1[,rank]
  S2 = S2[,rank]
  for(i in 2:(ncats-1)){
    ordered_pvals[i] = chisq.test(rbind(c(sum(S1[,i]), sum(S1[,(i+1):ncats])),
                                        c(sum(S2[,i]), sum(S2[,(i+1):ncats]))))$p.value
  }
  ## here it used to be ordered_pvals[96]
  ordered_pvals[ncats] = chisq.test(rbind(c(sum(S1[,i]), sum(S1[,ncats-1])),
                                       c(sum(S2[,i]), sum(S2[,ncats-1]))))$p.value
  
  # plot(unordered_pvals[rank], ordered_pvals)
  # plot(unordered_pvals, sapply(1:96, function(idx) ordered_pvals[which(rank == idx)]))
  ## 4. Order, again, the p-values
  (sapply(1:ncats, function(idx) ordered_pvals[which(rank == idx)])) ## added padjust in 20220208
}


iterative_chisqrt_wrapper <- function(object, slot_name='count_matrices_all'){
  counts <- slot(object, name = slot_name)
  iterative_chisqrt(counts[[1]],
                    counts[[2]])[1]
}


compute_effect_size_1 <- function(mat1, mat2){
  stopifnot(all(dim(mat1) == dim(mat2)))
  ## the rows in mat1 and mat2 correspond to the same patients
  sum((normalise_rw(mat1) - normalise_rw(mat2))**2)/nrow(mat1)
}

compute_effect_size_2 <- function(mat1, mat2){
  stopifnot(all(dim(mat1) == dim(mat2)))
  ## the rows in mat1 and mat2 correspond to the same patients
  try(sum((normalise_rw(mat1) - normalise_rw(mat2))**2)/(nrow(mat1))*log(mean(rowSums(rbind(mat1,mat2)))))
}


