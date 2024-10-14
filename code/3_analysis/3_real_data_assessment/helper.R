give_missclassified_matrix <- function(original_matrix, counts_to_misclassify, missclassify=T){
  
  ## Either missclassify mutations (put them from the clonal group to the subclonal, or vice-versa), -- if missclassify=T
  ## or remove signatures -- if missclassify=F
  
  cumsum_matrix = cumsum(original_matrix$Y) ## this is done by column
  
  new_dataset_obj_trinucleotide = original_matrix
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
    
    if(missclassify){
      ## add it to the row of the same patient but in the other group (clonal/subclonal)
      new_x_in_Y = ifelse(x_in_Y <= nrow(new_dataset_obj_trinucleotide$Y)/2, x_in_Y+nrow(new_dataset_obj_trinucleotide$Y)/2, x_in_Y-nrow(new_dataset_obj_trinucleotide$Y)/2)
      ## add it to the row of the same patient but in the other group (clonal/subclonal). Column remains the same.
      new_dataset_obj_trinucleotide$Y[new_x_in_Y,y_in_Y] = new_dataset_obj_trinucleotide$Y[new_x_in_Y,y_in_Y] + 1
    }
    ## subtract the count from the original row
    new_dataset_obj_trinucleotide$Y[x_in_Y,y_in_Y] = new_dataset_obj_trinucleotide$Y[x_in_Y,y_in_Y] - 1
  }
  return(new_dataset_obj_trinucleotide)
}

give_barplot_agreement_in_DA <- function(diagRE_DM_tests, diagDM_leave_one_out_exposures_tests, ylabel='Number of leave-one-out datasets', colour_version='colourversion1'){
  diagDM_leave_one_out_exposures_tests_accordance_DA = melt(sapply(enough_samples, function(ct) t(t(rep(diagRE_DM_tests[[ct]] <= 0.05, length(diagDM_leave_one_out_exposures_tests[[ct]])) + (diagDM_leave_one_out_exposures_tests[[ct]] <= 0.05))), simplify = F))
  diagDM_leave_one_out_exposures_tests_accordance_DA$accordance = ifelse(diagDM_leave_one_out_exposures_tests_accordance_DA$value == 1, 'Discordance', 'Accordance')
  
  if(colour_version=='colourversion1'){
    res <- ggplot(diagDM_leave_one_out_exposures_tests_accordance_DA, aes(x=L1, #alpha=(accordance== 'Accordance'),
                                                                          fill=ifelse(accordance== 'Accordance',
                                                                                      get_pcawg_name(L1),  'white'), colour='1'))+
      scale_fill_manual(values = c(pcawg_palette, "white"="white"), na.value="cyan")+
      scale_color_manual(values = c('1'='black'))+guides(fill='none', color='none')
  }else  if(colour_version=='colourversion2'){
    res <- ggplot(diagDM_leave_one_out_exposures_tests_accordance_DA, aes(x=L1, #alpha=(accordance== 'Accordance'),
                                                                          fill=ifelse(accordance== 'Accordance',
                                                                                      'Agreement',  'Disagreement'),
                                                                          colour=get_pcawg_name(L1)))+
      scale_fill_manual(values = c("Disagreement"="red", "Agreement"="#eaf7e9"), na.value="#fccfd1")+
      scale_color_manual(values = c(pcawg_palette))+guides(color='none')+labs(fill='Agreement in DA')
  }else{
    stop('Check <colour_version>')
  }
  res <- res + geom_bar()+
    theme(#axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(), legend.position = 'bottom')+
    labs(alpha='Accordance', x='Cancer types', y=ylabel)
  return(res)
}

get_pcawg_name <- function(L1){
  return(tolower(make.names(gsub("[.].*", "", L1))))
}

give_indices_beta <- function(j, items_per_category=2){
  c(1+(j-1)*2, 2+(j-1)*2)
}


give_beta_cor_given_ct <- function(df_from_it, diagRE_DMct, mode, j){
  # print(df_from_it)
  # print(diagRE_DMct)
  # 
  if(is.null(df_from_it)){
    NA
  }else{
    if( (class(df_from_it) == "try-error") | is.na(df_from_it[1])){
      NA
    }else{
      # if( (grepl('leave_one_out', mode) | grepl('add_one', mode) ) & nested & (j == length(df[[ct]]))){
      #   NA ## this is the last signature; it doesn't make any sense for beta coefs to be compared because we are using a different baseline
      #   ## this only applies to leave_one_out and add_one scenarios
      # }else{
        if(mode == 'intact'){ ## do not modify betas: plot as they are
          cor(plot_betas(df_from_it, return_df = T, plot = F)[,'Estimate'],
              plot_betas(diagRE_DMct, return_df = T, plot = F)[,'Estimate'])
        }else if(mode == 'intact_beta_slope'){ ## do not modify betas: plot as they are
          ## select only even coefficients
          cor(plot_betas(df_from_it, return_df = T, plot = F)[,'Estimate'][c(F,T)],
              plot_betas(diagRE_DMct, return_df = T, plot = F)[,'Estimate'][c(F,T)])
        
        }else if(mode == 'splitting_first_signature_all_betas'){
          ## remove the first two signatures from df, and the first signature from the original dataset diagRE_DMct
          stopifnot(give_indices_beta(c(1), 2) == c(1,2))
          cor(plot_betas(df_from_it, return_df = T, plot = F)[,'Estimate'][-give_indices_beta(c(1, 2), 2)],
              plot_betas(diagRE_DMct, return_df = T, plot = F)[,'Estimate'][-give_indices_beta(c(1), 2)])
        }else if(mode == 'splitting_first_signature_beta_slope'){
          ## select only even coefficients
          cor(plot_betas(df_from_it, return_df = T, plot = F)[,'Estimate'][-give_indices_beta(c(1, 2), 2)][c(F,T)],
              plot_betas(diagRE_DMct, return_df = T, plot = F)[,'Estimate'][-give_indices_beta(c(1), 2)][c(F,T)])
        
        }else if(mode == 'leave_one_out_all_betas'){
          cor(plot_betas(df_from_it, return_df = T, plot = F)[,'Estimate'],
              plot_betas(diagRE_DMct, return_df = T, plot = F)[,'Estimate'][-give_indices_beta(j, 2)])
        }else if(mode == 'leave_one_out_beta_slope'){
          ## select only even coefficients
          cor(plot_betas(df_from_it, return_df = T, plot = F)[,'Estimate'][c(F,T)],
              plot_betas(diagRE_DMct, return_df = T, plot = F)[,'Estimate'][-give_indices_beta(j, 2)][c(F,T)])
        }else if(mode == 'add_one'){
          cor(plot_betas(df_from_it, return_df = T, plot = F)[,'Estimate'][-c(1:2)], ## removing the first (newly added) signature
              plot_betas(diagRE_DMct, return_df = T, plot = F)[,'Estimate'])
        }else if(mode == 'add_one_beta_slope'){
          cor(plot_betas(df_from_it, return_df = T, plot = F)[,'Estimate'][-c(1:2)][c(F,T)], ## removing the first (newly added) signature
              plot_betas(diagRE_DMct, return_df = T, plot = F)[,'Estimate'][c(F,T)])
        }else{
          stop('<mode> not found')
        }
      # }
    }
  }
}

give_beta_cor <- function(df, diagRE_DM, mode='leave_one_out_all_betas', df_is_nested=F, diagRE_DM_is_nested=F){
  # df_is_nested: if there are replicates, df_is_nested. df should have this structure df[[ct]][[covariate]][[repl]]
  
  if( grepl('leave_one_out', mode) | grepl('add_one', mode) | grepl('intact', mode) | grepl('splitting_first_signature', mode)){
    ## correlation can only be computed for the first n-1 signatures, as if we leave out the last category we don't have a shared baseline
    
    if(diagRE_DM_is_nested & (!df_is_nested)){
      stop('For diagRE_DM_is_nested=T df_is_nested should be T\n')
    }
    ## this is with both betas, intercept and slope
    if(!df_is_nested){
      ## by default, df has two levels
      sapply(enough_samples, function(ct) sapply(1:length(df[[ct]]), function(j){
        if( (grepl('leave_one_out', mode) | grepl('add_one', mode) ) & (j == length(df[[ct]]))){
          #   NA ## this is the last signature; it doesn't make any sense for beta coefs to be compared because we are using a different baseline
          #   ## this only applies to leave_one_out and add_one scenarios
          NA
        }else{
          give_beta_cor_given_ct(df_from_it=df[[ct]][[j]],
                                 diagRE_DMct=diagRE_DM[[ct]], mode=mode, j=j)
        }
      }, simplify = F), simplify = F)
    }else{
      ## if df has three levels
      if(diagRE_DM_is_nested){
        sapply(enough_samples, function(ct) sapply(1:length(df[[ct]]), function(j)
          sapply(1:length(df[[ct]][[j]]), function(repl)
            give_beta_cor_given_ct(df_from_it=df[[ct]][[j]][[repl]],
                                   diagRE_DMct=diagRE_DM[[ct]][[repl]], mode=mode, j=j), simplify = F, USE.NAMES = T),
          simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
        
      }else{
        sapply(enough_samples, function(ct) sapply(1:length(df[[ct]]), function(j)
          sapply(1:length(df[[ct]][[j]]), function(repl)
            give_beta_cor_given_ct(df_from_it=df[[ct]][[j]][[repl]],
                                 diagRE_DMct=diagRE_DM[[ct]], mode=mode, j=j), simplify = F, USE.NAMES = T),
          simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
      }
      }
  }else{
    stop('Not implemented')
  }
}


give_beta_cor_wrapper_replicates <- function(repl, ...){
  sapply(df, )
}


give_plot_betacors_nested <- function(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP, diagRE_DM, mode, title_y, fold_decrease_nmuts, diagRE_DM_is_nested=F, level2_is_fold_decrease=T, gsub_from_x="", ncol_arg=8){
  
  cat('Getting correlations...\n')
  beta_cors = give_beta_cor(df = diagREDM_varying_nmuts_with_replicates_signature_exposuresQP,
                            diagRE_DM = diagRE_DM, mode = mode, df_is_nested = T, diagRE_DM_is_nested=diagRE_DM_is_nested)
  
  cat('Correlations have been computed\n')
  
  ## adding names
  for(ct in enough_samples){
    names(beta_cors[[ct]]) <- names(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP[[ct]])
    for(j in 1:length(fold_decrease_nmuts)){
      names(beta_cors[[ct]][[j]]) <- paste0('Repl', 1:nreplicates)
    }
  }
  
  cat('Plotting...\n')
  
  beta_cors_melt <- melt(beta_cors)
  beta_cors_melt$L2 <- gsub(gsub_from_x, "", beta_cors_melt$L2)
  
  if(level2_is_fold_decrease){
    beta_cors_melt$L2 <- signif(1/as.numeric(gsub("_fold_decrease_nmuts","" , beta_cors_melt$L2)), digits = 2)
    beta_cors_melt$L2 <- factor(beta_cors_melt$L2, levels = sort(unique(beta_cors_melt$L2), decreasing = T))
  # head(beta_cors_melt)
  }
  
  ggplot(beta_cors_melt, aes(x=L2, y=value, col=L1, group=interaction(L1,L2)))+geom_boxplot()+guides(col='none')+
    facet_wrap(.~L1, ncol=ncol_arg)+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )+
    scale_color_manual(values = pcawg_palette_2)+
    labs(x='Fraction of mutations kept',
         y=title_y,
         remove_col_legend=T, additional_title = ' in\nleave-one-out (D1A)')
  
}


give_plot_betacors_notnested <- function(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP, diagRE_DM, mode, title_y, fold_decrease_nmuts, level2_is_fold_decrease=T, gsub_from_x="", ncol_arg=8){
  
  cat('Getting correlations...\n')
  beta_cors = give_beta_cor(df = diagREDM_varying_nmuts_with_replicates_signature_exposuresQP,
                            diagRE_DM = diagRE_DM, mode = mode, df_is_nested = F)
  ## adding names
  for(ct in enough_samples){
    names(beta_cors[[ct]]) <- names(diagREDM_varying_nmuts_with_replicates_signature_exposuresQP[[ct]])
  }
  
  cat('Plotting...\n')
  
  beta_cors_melt <- melt(beta_cors)
  beta_cors_melt$L2 <- gsub(gsub_from_x, "", beta_cors_melt$L2)
  
  if(level2_is_fold_decrease){
    beta_cors_melt$L2 <- signif(1/as.numeric(gsub("_fold_decrease_nmuts","" , beta_cors_melt$L2)), digits = 2)
    beta_cors_melt$L2 <- factor(beta_cors_melt$L2, levels = sort(unique(beta_cors_melt$L2), decreasing = T))
  }
  head(beta_cors_melt)
  
  
  ggplot(beta_cors_melt, aes(x=L2, y=value, col=L1, group=interaction(L1,L2)))+geom_boxplot()+guides(col='none')+
    facet_wrap(.~L1, ncol=ncol_arg)+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )+
    scale_color_manual(values = pcawg_palette_2)+
    labs(x='Fraction of mutations kept',
         y=title_y,
         remove_col_legend=T, additional_title = ' in\nleave-one-out (D1A)')
  
}
