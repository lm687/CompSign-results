#########################################################
################### Work in progress ####################
#########################################################

objects_sigs_per_CT = lapply(it_features, function(type_data){
  .x = lapply(list_CT, function(cancer_type){
    savefolder_features = "../data/roo/"
    objects_sigs_per_CT_features <- list()
    fle <- paste0(savefolder_features, cancer_type, '_', it_features, '_ROO.RDS')
    print(fle)
    if(file.exists(fle)){
      objects_sigs_per_CT_features <- readRDS(fle)
    }else{
      cat('Object for cancer type', cancer_type, 'was not found')
    }
    
    if(type_data == "nucleotidesubstitution1"){
      objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
    }else if(type_data == "nucleotidesubstitution3"){
      objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
    }else if(type_data == "signatures"){
      if(is.null(attr(objects_sigs_per_CT_features,"count_matrices_active")[[1]]) | length(attr(objects_sigs_per_CT_features,"count_matrices_active")[[1]]) == 0){
        ## no active signatures
        if(! (is.null(attr(objects_sigs_per_CT_features,"count_matrices_all")[[1]]) | length(attr(objects_sigs_per_CT_features,"count_matrices_all")[[1]]) == 0)){
          objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
        }else{
          ## no type of exposure, either all or active
          return(NULL)
        }
      }else{
        objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_active")
      }
      if(ncol(objects_sigs_per_CT_features[[1]]) == 0){
        cat(cancer_type, 'does not have any mutations in one of the two groups. NULL is returned.\n')
        return(NULL) ## zero exposures in some group
      }
      objects_sigs_per_CT_features = lapply(objects_sigs_per_CT_features, function(i){
        rwn = rownames(i)
        if(nrow(i) > 1){
          .x = apply(i, 2, as.numeric)
        }else{
          cln = colnames(i)
          .x = t(matrix(as.numeric(i)))
          colnames(.x) = cln
        }
        rownames(.x) = rwn
        round(.x)
      })
    }
    objects_sigs_per_CT_features
  })
  names(.x) = list_CT
  .x
})
names(objects_sigs_per_CT) = it_features

cat('List object objects_sigs_per_CT loaded\n')
