plotcontour <- function(patient_idx, group_idx, model_name,  add=TRUE, true_contour, ...){
  if(group_idx == 1){
    patient_idx = patient_idx
  }else if(group_idx == 2){
    patient_idx = patient_idx +  npatientsx2/2
  }
  
  pcomp_idx = prcomp_all[[model_name]]
  groupsplit = splits_df[[model_name]]
  df_subset = (scale(groupsplit[[patient_idx]],
                     center = TRUE, scale = FALSE) %*% pcomp_idx$rotation)
  
  title = paste0(model_name, ' ', c('Early', 'Late')[group_idx])
  #paste0(c('Dirichlet-Multinomial', 'Multinomial')[model_name], ' ', c('Early', 'Late')[group_idx])
  if(true_contour){
    ## plot contour plot
    contour(kde2d(df_subset[,1],
                  df_subset[,2], n=50), drawlabels=FALSE, nlevels=10, add=add,  col='#5c5c5c', ...,
            main=title)
  }else{
    ## plot simulated points instead
    if(!add){
      plot(df_subset[,1],
           df_subset[,2], main=title,  pch=19, col=alpha('black', 0.02), ...)
    }else{
      points(df_subset[,1],
             df_subset[,2], main=title,  pch=19, col=alpha('black', 0.02),  ...)
    }
  }
  
}

plot_whole_contour = function(group_idx, model_name, true_contour=TRUE){
  projected_observed_idx = projected_observed[[model_name]]
  
  if(length(projected_observed_idx) == 1 & is.na(projected_observed_idx)){
    plot(0, 0)
  }else{
    lims = apply(rbind(prcomp_res[[model_name]], projected_observed_idx), 2, function(i) c(min(i, na.rm = TRUE), max(i, na.rm = TRUE)))
    
    plotcontour(patient_idx = 1, group_idx = group_idx, model_name = model_name, add = FALSE,
                xlim = c(lims[,1]), ylim = c(lims[,2]), true_contour = true_contour)
    if(npatientsx2 > 2){
      for(patient_idx in 2:(npatientsx2/2)){
        plotcontour(patient_idx = patient_idx, group_idx = group_idx, model_name = model_name, add = TRUE, true_contour)
      }
    }
    
    for(patient_idx in 1:(npatientsx2/2)){
      if(group_idx == 1){
        points(projected_observed_idx[patient_idx,1], projected_observed_idx[patient_idx,2], col='red', pch=19)
      }else if(group_idx == 2){
        points(projected_observed_idx[patient_idx +  npatientsx2/2,1], projected_observed_idx[patient_idx +  npatientsx2/2,2],
               col='red', pch=19)
      }
    }
  }
}