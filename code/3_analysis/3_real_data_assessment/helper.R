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
