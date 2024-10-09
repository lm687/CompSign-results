##' functions:
##' [7] QPsig
##' [33] complementbase
##' [48] change_muts
##' [59] createRDS_ROOSigs_object
##' [172] createROO_ROOSigs


QPsig<-function(tumour.ref = NA,signatures.ref){
  ## modified from Lynch, 2016
  # we normalize the observations so that they sum to 1
  W <- table(factor(tumour.ref, levels = rownames(signatures.ref)))

  obs<-as.numeric(W/sum(W))
  # to allow use of the deconstructSigs objects we convert to matrices
  signatures.ref<-t(as.matrix(signatures.ref))
  # we use the factorized version of solve.QP -
  # see the helpfile of that function for details of the required form
  # otherwise we would set Dmat = signatures.ref %*% t(signatures.ref) as indicated # in the article
  Rinverse <- backsolve(chol(signatures.ref %*% t(signatures.ref)),
                        diag(dim(signatures.ref)[1]))
  # we also need to define the linear part of the function that we are minimizing
  dvec <- (obs) %*% t(signatures.ref)
  # we have one constraint that the sum of weights is equal to 1 # we have more constraints that each weight is positive
  Amat <- cbind(rep(1,dim(Rinverse)[1]), diag(dim(Rinverse)[1]))
  bvec <- c(1,rep(0,dim(Rinverse)[1]))
  # we now call the solve.QP function from the quadprog library
  myQP<-quadprog::solve.QP(Dmat = Rinverse, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1, factorized = TRUE)
  result <- myQP$solution
  names(result) <- rownames(signatures.ref)
  return(result)
}

## change if the mutations are G or A
complementbase <- function(mut){
  if(mut == 'A'){
    'T'
  }else if(mut == 'T'){
    'A'
  }else if(mut == 'C'){
    'G'
  }else if(mut == 'G'){
    'C'
  }else{
    NA
  }
}

## change if the mutations are G or A
change_muts <- function(mut){
  .x <- strsplit(mut, '[|]|>')
  .x <- strsplit(mut, '\\[|\\]|>')
  if(.x[[1]][2] %in% c('G', 'A')){
    paste0(complementbase(.x[[1]][4]), '[', complementbase(.x[[1]][2]), '>', complementbase(.x[[1]][3]), ']',  complementbase(.x[[1]][1]) )
  }else{
    mut
  }
}

bases = c('T', 'C', 'G', 'A')

## dataframe with mutation equivalencies (for the purposes of speeding up code)
grid_trinucleotides = expand.grid(data.frame(replicate(3, list('A', 'C', 'G', 'T')) ))
grid_mutations = sort(as.vector(apply(grid_trinucleotides, 1, function(i){sapply(bases[i[2] != bases], function(j){ return(paste0(c(i[1], '[', i[2], '>', j, ']', i[3]), collapse = ""))})})))
grid_mutations_df = cbind(grid_mutations, sapply(grid_mutations, change_muts))

createRDS_ROOSigs_object <- function(pre_path="",
                                     vcf_path="../data/restricted/pcawg/pcawg_restricted_snv_counts/",
                                     ccf_threshold=NULL,
                                     type_features=NULL,
                                     active_sigs_version="active_signatures_transposed_clonesig2",
                                     outfolder=NULL,
                                     cancer_type=NULL,
                                     in_dataframe=NULL,
                                     type_signatures="QP_COSMIC",
                                     check_size_exposures=T){

  if(is.null(cancer_type)){stop('Cancer type cannot be null')}

  j = outfolder
  cancer_types = cancer_type
  
  ## we need to know:
  ## 1. name of file
  ## 2. cancer type it belongs to for active signatures
  
  if((type_features %in% c('signatures', 'signaturesPCAWG')) & (type_signatures=="QP_COSMIC") ){
    sigs_cosmic <- read.table(paste0("../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
                              stringsAsFactors = FALSE, sep = ',', header = TRUE)
    # sigs_cosmic <- read.table(paste0(pre_path, "../../data/cosmic/signatures_probabilities.txt"),
    #                           stringsAsFactors = FALSE, sep = '\t', header = TRUE)
    # sigs_cosmic <- sigs_cosmic[,!(apply(sigs_cosmic, 2, function(x) sum(is.na(x))/length(x)) == 1)] ## remove empty columns which for whatever reason are there
    # rownames(sigs_cosmic) <- sigs_cosmic[,3]
    rownames(sigs_cosmic) <- apply(sigs_cosmic, 1, function(i) paste0(substr(i[2], 1, 1), '[', i[1], ']', substr(i[2], 3, 3)))
    # sigs_cosmic <- sigs_cosmic[,-c(1:3)]
    sigs_cosmic <- sigs_cosmic[,-c(1:2)]
    
    ## read in active signatures
    active <- read.table(paste0("../data/cosmic/", active_sigs_version),
                         stringsAsFactors = FALSE, sep = '\t', header = TRUE)
    subset_active <- active[,-c(1, 2, ncol(active))]
    
    if(! ("id2" %in% colnames(active))){
      stop('The dataframe of active signatures should have a column called <id2> with the cancer type identifier')
    }
    # present_active <- sapply(active$acronym, function(i) if(grepl(',', i)){strsplit(i, ', ')}else{i}); present_active[present_active == ""] <- NULL; present_active <- as.character(unlist(present_active))
    present_active <- sapply(active$id2, function(i) if(grepl(',', i)){strsplit(i, ', ')}else{i})
    present_active[present_active == ""] <- NULL; present_active <- as.character(unlist(present_active))
    
  }
  
  ## below: file which should have active signatures file but does not
  
  try({
    # name=gsub('/out_', '', gsub('.consensus.20160830.somatic.snv_mnv.vcf_merged', '', gsub(outfolder, '', j)))
    
    
    if(is.null(in_dataframe)){stop('<in_dataframe> should not be empty - this is your data!')}
    x = in_dataframe
    
    if(type_features != 'signaturesmutSigExtractor'){

      x$mutation <- gsub('/', '>', x$mutation)
      x$mutation <- grid_mutations_df[match(x$mutation, grid_mutations_df[,1]),2]
      
      if(type_features == "nucleotidesubstitution1"){
        x$mutation <- sapply(x$mutation, function(i){strsplit(i, "\\[|\\]")[[1]][2]})
      }
    }else{
      # is signaturesmutSigExtractor
      x$mutation = NA
    }
    ## bin by clonal/subclonal
    num_muts = nrow(x)
    Clonal <- x %>% filter(x$bool_group_1)
    Subclonal <- x %>% filter(!x$bool_group_1)
    
    ## infer signatures
    if(type_features %in% c('signatures', 'signaturesPCAWG')){
      sigsClonal_QP <- length(Clonal$mutation)*QPsig(tumour.ref = Clonal$mutation, signatures.ref = sigs_cosmic)
      sigsSubclonal_QP <- length(Subclonal$mutation)*QPsig(Subclonal$mutation, signatures.ref = sigs_cosmic)
    }else if(type_features == 'signaturesmutSigExtractor'){
      library(mutSigExtractor)
      library(BSgenome.Hsapiens.UCSC.hg19)
      
      ## Clonal
      contexts_snv_clonal <- extractSigsSnv(df=Clonal[,c('chrom','pos','ref','alt')], output='contexts',
                                     ref.genome=BSgenome.Hsapiens.UCSC.hg19)
      sigs_snv_clonal <- fitToSignatures(
        mut.context.counts=contexts_snv_clonal[,1], 
        signature.profiles=SBS_SIGNATURE_PROFILES_V3
      )
      ## Subclonal
      contexts_snv_subclonal <- extractSigsSnv(df=Subclonal[,c('chrom','pos','ref','alt')], output='contexts',
                                            ref.genome=BSgenome.Hsapiens.UCSC.hg19)
      sigs_snv_subclonal <- fitToSignatures(
        mut.context.counts=contexts_snv_subclonal[,1], 
        signature.profiles=SBS_SIGNATURE_PROFILES_V3
      )
      
      sigsClonal_QP <- sigs_snv_clonal
      sigsSubclonal_QP <- sigs_snv_subclonal
    }else{
      if(type_features == "nucleotidesubstitution1"){
        lvls <- sort(unlist( sapply(c('C', 'T'), function(base){ paste0(base, '>', c('A', 'C', 'G', 'T')[! (c('A', 'C', 'G', 'T') == base)]) })))
      }else if(type_features == "nucleotidesubstitution3"){
        lvls <- sort(unlist( sapply(c('C', 'T'), function(base) sapply(sapply(sapply(c('A', 'C', 'G', 'T'), paste0, '[', base), paste0, '>', c('A', 'C', 'G', 'T')[! (c('A', 'C', 'G', 'T') == base)]), paste0, ']', c('A', 'C', 'G', 'T')))))
      }else{
        stop('Incorrect type_features')
      }
      sigsClonal_QP <- table(factor(Clonal$mutation, levels = lvls))
      sigsSubclonal_QP <- table(factor(Subclonal$mutation, levels = lvls))
    }      

    return_object <- list(sigsClonal_QP, sigsSubclonal_QP)
    
    if(check_size_exposures){
      if( ! all.equal((sum(sigsClonal_QP) + sum(sigsSubclonal_QP)), num_muts) ){
        stop('The number of allocated mutations is not the number of original mutations')
      }
    }
    
    ## are there any active signatures?
    if(type_features %in% c('signatures', 'signaturesPCAWG')){
      ## that depends on whether we are reading from the file or not

      cancer_type_abbrev = toupper(strsplit(cancer_type, '[.]')[[1]][1])
      if(cancer_type_abbrev %in% present_active){
        active_j <- colnames(subset_active)[which(subset_active[grepl(cancer_type_abbrev, active$id2),] == 1)]
        sigsClonal_QP_active <- length(Clonal$mutation)*QPsig(Clonal$mutation, signatures.ref = sigs_cosmic[,active_j])
        sigsSubclonal_QP_active <- length(Subclonal$mutation)*QPsig(Subclonal$mutation, signatures.ref = sigs_cosmic[,active_j])
        return_object <- list(return_object, list(sigsClonal_QP_active, sigsSubclonal_QP_active))
      }else{
        cat(paste0('Cancer type', cancer_type_abbrev, ' does not have active signatures in file <', active_sigs_version, '>'))
        return_object <- list(return_object, list(list(), list()))
      }
    }else if(type_features == "signaturesmutSigExtractor"){
      # only populate the active signatures slot
      return_object <- list(return_object, list(list(), list()))
    }
  })
  if(!exists("return_object")){
    ## it doesn't exist because it did not appear in the cancer types
  }else{
    return_object
  }
}



# ## ROO for PCAWG data
# createROO_ROOSigs <- function(savefolder, type_features='nucleotidesubstitution3', read_in=TRUE, all_objs_activesigs=NULL, save_bool){
#   # rm(list = ls())
#   # setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#   # Sys.setenv(LANG='en')
#   
#   library(dplyr)
#   library(reshape2)
#   source(paste0(pre_path, "roo_functions.R"))
#   
#   outfolder="../data/roo/"
#   getname <- function(j)  gsub("−", "-", gsub('out_', '', gsub('.consensus.20160830.somatic.snv_mnv.vcf_merged', '', gsub(outfolder, '', j))))
#   
#   ## Cancer types
#   
#   ## iterate over cancer types
#   ##   iterate over samples from this cancer type
#   
#   objects_sigs_per_CT <- list()
#   for(i in unique(cancer_types)){
#     i_matrix_all <- list(do.call('rbind', lapply(all_objs[names(subset(cancer_types, cancer_types==i))], function(k) k[[1]])),
#                          do.call('rbind', lapply(all_objs[names(subset(cancer_types, cancer_types==i))], function(k) k[[2]])))
#     i_matrix_active <- list(do.call('rbind', lapply(all_objs_activesigs[paste0(names(subset(cancer_types, cancer_types==i)), '_active')], function(k) k[[1]])),
#                             do.call('rbind', lapply(all_objs_activesigs[paste0(names(subset(cancer_types, cancer_types==i)), '_active')], function(k) k[[2]])))
#     nulls_active <- (sapply(i_matrix_active, is.null))
#     if(sum(nulls_active) == 2){
#       is_null_active <- TRUE
#     }else if(sum(nulls_active) == 1){
#       stop('You should not have active signatures for only one of the categories')
#     }else{
#       is_null_active <- FALSE
#     }
#     objects_sigs_per_CT[[i]] <- new("exposures_cancertype",
#                                     cancer_type=i,
#                                     type_classification = "trunk_vs_branch",
#                                     number_categories = 2,
#                                     id_categories = c('trunk', 'branch'),
#                                     active_signatures = "character", ## active signatures for this cancer type
#                                     count_matrices_all = i_matrix_all,
#                                     count_matrices_active = i_matrix_active,
#                                     fitted_distrib_count_matrices_all = vector("list", 2),
#                                     fitted_distrib_count_matrices_active = vector("list", 2),
#                                     sample_names = names(subset(cancer_types, cancer_types==i)),
#                                     modification = "none",
#                                     is_null_active = is_null_active
#     )
#   }
#   
#   ## save
#   if(save_bool){
#     for(fle_it in 1:length(objects_sigs_per_CT)){
#       saveRDS(objects_sigs_per_CT[[fle_it]], file = paste0(savefolder, names(objects_sigs_per_CT)[fle_it], '_ROOSigs.RDS'))
#     }
#   }else{
#     objects_sigs_per_CT
#   }
#   
# }


createROO_ROOSigs_object <- function(type_features='nucleotidesubstitution3',
                                     all_objs_activesigs=NULL, save_bool, outfiles, all_objs, pre_path="",
                                     cancer_type_given=NULL, file_name_given=NULL){
  
  # source(paste0(pre_path, "../../../../other_repos/Pseudo-Ordering-Signatures/code/pcawg/helper_functions_pcawg.R"))
  library(dplyr)
  library(reshape2)

  type_features=opt$feature_type
  all_objs=rds_object
  all_objs_activesigs=all_objs_activesigs
  save_bool=FALSE
  outfiles=names(rds_object)
  pre_path="1_create_ROO/"
  cancer_type_given = rep(opt$cancer_type, length(rds_object))
  file_name_given=NULL
  
  source(paste0(pre_path, "roo_functions.R"))
  
  
  outfolder <- "../../data/roo/"
  getname <- function(j)  gsub("−", "-", gsub('out_', '', gsub('.consensus.20160830.somatic.snv_mnv.vcf_merged', '', gsub(outfolder, '', j))))
  
  ## Cancer types
  if(!is.null(cancer_type_given)){
    ## We are telling it what the cancer type is
    cancer_types = cancer_type_given
    outfiles_name <- basename(outfiles)
  }else{
    stop()
  }
  names(cancer_types) <- outfiles_name
  stopifnot(length(unique(cancer_types)) == 1)
  ## iterate over cancer types
  ##   iterate over samples from this cancer type
  
  objects_sigs_per_CT <- list()
  
  i_matrix_all <- list(do.call('rbind', lapply(all_objs, function(k) k[[1]])),
                       do.call('rbind', lapply(all_objs, function(k) k[[2]])))
  i_matrix_active <- list(do.call('rbind', lapply(all_objs_activesigs, function(k) k[[1]])),
                          do.call('rbind', lapply(all_objs_activesigs, function(k) k[[2]])))
  if(sum(sapply(all_objs_activesigs, function(i) length(i[[1]]) + length(i[[2]]))) == 0){
    ## no active signatures
    i_matrix_active <- list(list(), list())
  }
  nulls_active <- (sapply(i_matrix_active, is.null))
  if(sum(nulls_active) == 2){
    is_null_active <- TRUE
  }else if(sum(nulls_active) == 1){
    stop('You should not have active signatures for only one of the categories')
  }else{
    is_null_active <- FALSE
  }
  objects_sigs_per_CT <- new("exposures_cancertype",
                                  cancer_type=cancer_types,
                                  type_classification = "trunk_vs_branch",
                                  number_categories = 2,
                                  id_categories = c('trunk', 'branch'),
                                  active_signatures = "character", ## active signatures for this cancer type
                                  count_matrices_all = i_matrix_all,
                                  count_matrices_active = i_matrix_active,
                                  sample_names = names(cancer_types),
                                  modification = "none",
                                  is_null_active = is_null_active,
                                  is_empty="Non-empty"
  )

  
  ## save
  if(save_bool){
      saveRDS(objects_sigs_per_CT, file = paste0(savefolder, names(objects_sigs_per_CT)[fle_it], '_ROOSigs.RDS'))
  }else{
    objects_sigs_per_CT
  }
  
}

give_dummy_row_names = function(df, prefix='Sample'){
  rownames(df) = paste0(prefix, 1:nrow(df))
  df
}

give_dummy_col_names = function(df, prefix='Feature'){
  colnames(df) = paste0(prefix, 1:ncol(df))
  df
}
