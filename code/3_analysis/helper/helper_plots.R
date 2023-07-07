createBarplot <- function(matrix_exposures, angle_rotation_axis = 0, order_labels=NULL,
                          remove_labels=FALSE, levels_signatures=NULL, includeMelt=NULL, Melt=NULL, verbose=TRUE){
  #' error due to it not being a matrix:  
  #'    No id variables; using all as measure variables
  #'    Rerun with Debug
  #'    Error in `[.data.frame`(.mat, , "Var1") : undefined columns selected \\
  #' (use tomatrix())
  
  
  if(is.null(colnames(matrix_exposures))) stop('columns must have names')
  require(reshape2)
  require(ggplot2)
  library(RColorBrewer)
  if(verbose){
    cat(paste0('Creating plot... it might take some time if the data are large. Number of samples: ', nrow(matrix_exposures), '\n'))
  }
  
  if( (!is.null(order_labels) & typeof(order_labels) == "logical")){if(!order_labels){cat('WARNING: Order labels is either a vector with desired order or NULL, not bool')}}
  
  if(!is.null(levels_signatures)){
    library(ggthemes)
    ggthemes_data$economist
  }else{
    levels_signatures <- colnames(matrix_exposures)
  }
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  
  myColors <- col_vector[1:length(levels_signatures)]
  names(myColors) <- levels_signatures
  if(is.null(order_labels)) order_labels = rownames(matrix_exposures)
  if(!is.null(includeMelt)){
    cat('For whatever reason sometimes the melt does not work. Here it is passed as argument.')
    .mat <- Melt
  }else{
    .mat <- melt(matrix_exposures)
  }
  .mat[,'Var1'] <- factor(.mat[,'Var1'], levels=order_labels)
  .mat[,'Var2'] <- factor(.mat[,'Var2'], levels=levels_signatures)
  ###rownames(.mat) <- rownames(matrix_exposures) ### new
  
  if(!remove_labels){
    if(!is.null(levels_signatures)){
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        #theme(axis.title.x=element_blank(),
        #      axis.text.x=element_blank(),
        #      axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))+
        scale_fill_manual(name = "grp",values = myColors)
    }else{
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        #theme(axis.title.x=element_blank(),
        #      axis.text.x=element_blank(),
        #      axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))
    }
  }else{
    if(!is.null(levels_signatures)){
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))+
        scale_fill_manual(name = "grp",values = myColors)
    }else{
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))
    }
  }
}

is_empty_active = function(x){
  ## given an ROO object, determine whether the active signature matrix is empty
  sum(sapply(attr(x, "count_matrices_active"), length)) == 0
}

wrapper_barplot_ROO = function(x, bool_normalised=F, ...){
  if(is_empty_active(x)){
    matrices = attr(x, "count_matrices_all")
  }else{
    matrices = attr(x, "count_matrices_active")
  }
  if(bool_normalised){
    matrices = lapply(matrices, function(i) sweep(x = i, MARGIN = 1, STATS = rowSums(i), FUN = '/'))
  }
  lapply(matrices, createBarplot, ...)
}

