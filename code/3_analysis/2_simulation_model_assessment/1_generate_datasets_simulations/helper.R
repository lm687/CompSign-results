simulate_mutations_from_signatures <- function(exposure_mat){
  
  muts_mat <- t(apply(exposure_mat, 1, function(j){ 
    rowSums(sapply(1:ncol(exposure_mat), function(sig_it){
      table(factor(sapply(j[sig_it], function(i) sample(x = 1:nrow(signature_definitions), size = i, replace = T,
                                                        prob = signature_definitions[,colnames(exposure_mat)[sig_it]]) ),
                   levels=1:nrow(signature_definitions)))
    }))
  }))
  colnames(muts_mat) <- rownames(signature_definitions)
  if(!is.null(rownames(exposure_mat))){
    rownames(muts_mat) <- rownames(exposure_mat)
  }
  
  return(muts_mat)
}