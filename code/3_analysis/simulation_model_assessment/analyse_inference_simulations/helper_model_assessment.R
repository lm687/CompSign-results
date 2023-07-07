give_accuracies_with_varying_var <- function(var, two_var=F, datasets_arg=datasets,
                                             pvals_data_frame_arg=pvals_data_frame, single_pval=F, CN=F){
  if(two_var){
    do.call('rbind', apply(expand.grid(sapply(var, function(i) unique(unlist(sapply(datasets_arg, `[`, i))))), 1, function(vars_it){
      .pvals <- pvals_data_frame_arg[which(sapply(datasets_arg, '[', var[1]) == vars_it[[1]] & sapply(datasets_arg, '[', var[2]) == vars_it[[2]]),]
      if(single_pval){
        .res_all_subset = summarise_DA_detection(true = .pvals$true, predicted = .pvals$test <= 0.05)
        .res_all_subset <- t(cbind.data.frame(.res_all_subset))
        rownames(.res_all_subset) <- 'test'
      }else{
        if(CN){
          .res_all_subset = give_res_all_CN(.pvals)
        }else{
          .res_all_subset = put_vals_in_table(.pvals)
        }
      }
      .return <- cbind.data.frame(.res_all_subset, VAR1=vars_it[1], VAR2=vars_it[2], model=rownames(.res_all_subset))
      colnames(.return)[(ncol(.return)-2)] <- var[1]
      colnames(.return)[(ncol(.return)-1)] <- var[2]
      return(.return)
    }))    
  }else{
    do.call('rbind', lapply(unique(unlist(sapply(datasets_arg, `[`, var))), function(vars_it){
      .pvals <- pvals_data_frame_arg[which(sapply(datasets_arg, '[', var) == vars_it),]
      if(single_pval){
        .res_all_subset = (summarise_DA_detection(true = .pvals$true, predicted = .pvals$test <= 0.05))
        .res_all_subset <- t(cbind.data.frame(.res_all_subset))
        rownames(.res_all_subset) <- 'test'
      }else{
        if(CN){
          .res_all_subset = give_res_all_CN(.pvals)
        }else{
          .res_all_subset = put_vals_in_table(.pvals)
        }
      }
      .return <- cbind.data.frame(.res_all_subset, d=vars_it, model=rownames(.res_all_subset))
      colnames(.return)[(ncol(.return)-1)] <- var
      return(.return)
    }))
  }
}

give_res_all <- function(pvals_df){
  warning("This was wrong until 20220202. fullM and fullREDMSL were exchanged")
  rbind(fullREM=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$pvals_fullREM <= 0.05),
        fullREDMSL=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$pvals_fullREDMSL <= 0.05),
        diagREDMSL=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$pvals_diagREDMSL <= 0.05),
        diagREDM=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$pvals_diagREDM <= 0.05),
        Harris_chi=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$Harris_chi <= 0.05),
        ttest=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$ttest_props <= 0.05),
        ILR=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$ttest_ilr_adj <= 0.05),
        HMP=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$HMP <= 0.05),
        HMP2=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$HMP2 <= 0.05),
        perturbation=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$perturbation <= 0.05),
        permutation=summarise_DA_detection(true = pvals_df$true, predicted = pvals_df$permutation <= 0.05))
}


put_vals_in_table <- function(.pvals){
  warning('<fullREM> and <fullREDMSL> were interchanged until 20220201')
  rbind(fullREM=summarise_DA_detection(true = .pvals$true, predicted = .pvals$pvals_fullREM <= 0.05),
        fullREDMSL=summarise_DA_detection(true = .pvals$true, predicted = .pvals$pvals_fullREDMSL <= 0.05),
        diagREDMSL=summarise_DA_detection(true = .pvals$true, predicted = .pvals$pvals_diagREDMSL <= 0.05),
        diagREDM=summarise_DA_detection(true = .pvals$true, predicted = .pvals$pvals_diagREDM <= 0.05),
        Harris_chi=summarise_DA_detection(true = .pvals$true, predicted = .pvals$Harris_chi <= 0.05),
        ttest=summarise_DA_detection(true = .pvals$true, predicted = .pvals$ttest_props <= 0.05),
        ILR=summarise_DA_detection(true = .pvals$true, predicted = .pvals$ttest_ilr_adj <= 0.05),
        HMP=summarise_DA_detection(true = .pvals$true, predicted = .pvals$HMP <= 0.05),
        HMP2=summarise_DA_detection(true = .pvals$true, predicted = .pvals$HMP2 <= 0.05),
        perturbation=summarise_DA_detection(true = .pvals$true, predicted = .pvals$perturbation <= 0.05),
        permutation=summarise_DA_detection(true = .pvals$true, predicted = .pvals$permutation <= 0.05))
}

# colours_models <- c(diagREDM= "#943CB4", diagREDMSL= "#194D44", fullREDMSL=  "#C6CF6E",
#                     fullREM= "#5B6DC8", HMP= "#3CA437", HMP2= "#6B244C" ,
#                     ILR= "#6ACDC5", permutation= "#DE1A1A" , perturbation= "#BBB53E",
#                     Harris_chi= "#2A297A", ttest=  "#995533"   )
# colours_models2 <- c(diagREDM= "#943CB4", diagREDMSL= "#194D44", fullREDMSL=  "#C6CF6E",
#                     fullREM= "#5B6DC8", HMP= "#3CA437", HMP2= "#6B244C" ,
#                     ILR= "#6ACDC5", permutation= "#DE1A1A" , perturbation= "#BBB53E",
#                     Harris_chi= "#2A297A", ttest=  "#995533"   )
colours_models <- c(diagREDM= "#943CB4", diagREDMSL= "#194D44", fullREDMSL=  "#C6CF6E",
                    fullREM= "#e69138", HMP= "#3CA437", HMP2= "#6B244C" ,
                    ILR= "#6ACDC5", permutation= "#8fce00" , perturbation= "#BBB53E",
                    Harris_chi= "#2A297A", ttest=  "#995533"   )
colours_models2 <- c(diagREDM= "#943CB4", diagREDMSL= "#194D44", fullREDMSL=  "#C6CF6E",
                     fullREM= "#e69138", HMP= "#3CA437", HMP2= "#6B244C" ,
                     ILR= "#6ACDC5", permutation= "#8fce00" , perturbation= "#BBB53E",
                     Harris_chi= "#2A297A", ttest=  "#995533"   )
