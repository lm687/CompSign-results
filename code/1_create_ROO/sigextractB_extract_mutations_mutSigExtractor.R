##' Note: file helper_1_create_ROO.R also has a function to create mutsigextractor signatures,
##' even though they were not implemented using
##' the subset fo active signatures

debug_bool <- F

library(optparse)
if(debug_bool){
  rm(list = ls())
  setwd("~/Documents/PhD/GlobalDA/code/")
  opt=list();
  # opt$input_files='../data/restricted/pcawg/pcawg_restricted_snv_counts/00b9d0e6-69dc-4345-bffd-ce32880c8eef ../data/restricted/pcawg/pcawg_restricted_snv_counts/02917220-6a7a-46a1-8656-907e96bef88e ../data/restricted/pcawg/pcawg_restricted_snv_counts/03ad38a6-0902-4aaa-84a3-91ea88fa9883 ../data/restricted/pcawg/pcawg_restricted_snv_counts/05616329-e7ba-4efd-87b1-d79cd0f7af3d ../data/restricted/pcawg/pcawg_restricted_snv_counts/068f4f69-d2fe-4f25-912e-ca7d4623efb6 ../data/restricted/pcawg/pcawg_restricted_snv_counts/0e7f46ca-6f5c-4538-b6d6-00af65d57fcf ebe0ed67-2d3f-45cd-8f9b-4912595b16a0 ../data/restricted/pcawg/pcawg_restricted_snv_counts/fa676301-902f-473f-8313-5bff34ae549a ../data/restricted/pcawg/pcawg_restricted_snv_counts/fce8d8c6-f2a0-43a8-9a7a-b9c519a6686c'
  # opt$input_files='../data/restricted/pcawg/pcawg_restricted_snv/00b9d0e6-69dc-4345-bffd-ce32880c8eef.consensus.20160830.somatic.snv_mnv.vcf.gz ../data/restricted/pcawg/pcawg_restricted_snv/02917220-6a7a-46a1-8656-907e96bef88e.consensus.20160830.somatic.snv_mnv.vcf.gz'
  opt$cancer_type="Lymph-BNHL";
  opt$output = "../data/roo/Lymph-BNHL_signaturesMSE_ROO.RDS"
}else{
  option_list = list(
    make_option(c("--cancer_type"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--output"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--input_files"), type="character", default=NA,
                help="", metavar="character")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
}

## get all the vcf files of this cancer type
meta_pcawg <- read.table("../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt", header = T)

meta_pcawg <- meta_pcawg[meta_pcawg$histology_detailed == opt$cancer_type,]

opt$input_files <- paste0("../data/restricted/pcawg/pcawg_restricted_snv/", meta_pcawg$samplename, ".consensus.20160830.somatic.snv_mnv.vcf.gz")

# opt$input_files <- strsplit(opt$input_files, ' ')[[1]]

source("2_inference_TMB/helper_TMB.R")
source("3_analysis/validation/helper_load_vcf_and_clonal.R")
source("1_create_ROO/roo_functions.R") # exposures_cancertype class defintion

# devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')
library(mutSigExtractor)
library(BSgenome.Hsapiens.UCSC.hg19)
library(vcfR)

flder_clonal = "../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/"
flder_vcf = "../data/restricted/pcawg/pcawg_restricted_snv_counts/"
fles_cancer_type = paste0(basename(strsplit(opt$input_files, " ")[[1]]), "_cluster_assignments.txt.gz")

#' https://github.com/UMCUGenetics/mutSigExtractor

##' Problem: the vcf from PCAWG is from their own aasembly hs37d5.
##' Looking at genome.fa, it looks as though it's GRCh37 (i.e. hg19)
##' mutSigExtractor uses as a default genome hg19, which we specify here

## read ct and information

## read vcf
# list.files("../../../data/restricted/pcawg/pcawg_restricted_snv/")

contexts_snv_listCLONAL <- list()
contexts_snv_listSUBCLONAL <- list()

for(vcf_snv in  opt$input_files){
  
  if(file.exists(vcf_snv)){
    sample_id <- gsub(".consensus.20160830.somatic.snv_mnv.vcf.gz", "", basename(vcf_snv))
    
    # filecluster <-  "../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/0a6be23a-d5a0-4e95-ada2-a61b2b5d9485_cluster_assignments.txt.gz"
    # fileccf <- "../data/restricted/pcawg/consensus_subclonal_reconstruction_mutccf_20170325/0a6be23a-d5a0-4e95-ada2-a61b2b5d9485_mutation_ccf.txt"
    filecluster <-  paste0("../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/",
                           sample_id, "_cluster_assignments.txt.gz")
    filecluster
  
    # clonal <- (read.table(filecluster, header=T))
    # ccf <- (read.table(fileccf, header=T))
    # clonal
    
    # dim(ccf)
    # dim(clonal)
    
    ## read information about clonal/subclonal
    unified_vcf_and_clonal <- give_unified_vcf_and_ccf(vcf_file_path = vcf_snv,
                                                       subclonality_df_path = filecluster, 
                                                       type_clonal_file = 'cluster')
    
    cluster_colnames <- colnames(unified_vcf_and_clonal$subclonality_df)[grepl('cluster_', colnames(unified_vcf_and_clonal$subclonality_df))]
    
    if(length(cluster_colnames) == 1){
      cat('Only one group\n')
      next
    }
    
    clusters_MAP0 <- apply(unified_vcf_and_clonal$subclonality_df[,cluster_colnames], 1, which.max)
    if(any(sapply(clusters_MAP0, typeof) != "integer")){
      clusters_MAP0[sapply(clusters_MAP0, typeof) != "integer"]
      which(sapply(clusters_MAP0, typeof) != "integer")
    }
    clusters_MAP0_check <- sapply(clusters_MAP0, function(i) i %in% c(1,2,3))
    clusters_MAP0_check <- sapply(clusters_MAP0, length)
    clusters_MAP0[clusters_MAP0_check == 0] <- NA
    clusters_MAP0 <- unlist(clusters_MAP0)
    
    compare_ccf_and_cluster = F
    if(compare_ccf_and_cluster){
      ## read ccf
      fileccf <- paste0("../data/restricted/pcawg/consensus_subclonal_reconstruction_mutccf_20170325/",
                        sample_id, "_mutation_ccf.txt")
      
      unified_vcf_and_ccf <- give_unified_vcf_and_ccf(vcf_file_path = vcf_snv,
                                                      subclonality_df_path = fileccf,
                                                      type_clonal_file = 'ccf')
      
      image(as(unified_vcf_and_clonal$subclonality_df[,c('cluster_1', 'cluster_2', 'cluster_3')], 'matrix'))
      
      nrow((unified_vcf_and_clonal$subclonality_df[,c('cluster_1', 'cluster_2', 'cluster_3')])[unified_vcf_and_clonal$subclonality_df[,c('cluster_1')] > 0.5,])
      table(clusters_MAP0)
      
      plot(clusters_MAP0, unified_vcf_and_ccf$subclonality_df$ccf)
      nrow(unified_vcf_and_clonal$vcf)
      length(clusters_MAP0)
    }
    
    ## split vcf into clonal and subclonal mutations
    unified_vcf_and_clonal$vcf[unlist(clusters_MAP0) == 1,] ## clonal
    unified_vcf_and_clonal$vcf[unlist(clusters_MAP0) != 1,] ## subclonal
    filenameclonal <- gsub(".vcf", "_CLONAL.vcf", gsub("/pcawg_restricted_snv/", "/pcawg_restricted_snv_split/" , vcf_snv))
    filenamesubclonal <- gsub(".vcf", "_SUBCLONAL.vcf", gsub("/pcawg_restricted_snv/", "/pcawg_restricted_snv_split/" , vcf_snv))
    unified_vcf_and_clonalCLONAL <- unified_vcf_and_clonal$vcf[unlist(clusters_MAP0) == 1,]
    unified_vcf_and_clonalSUBCLONAL <- unified_vcf_and_clonal$vcf[unlist(clusters_MAP0) != 1,]
    
    unified_vcf_and_clonalCLONAL <- unified_vcf_and_clonalCLONAL[!is.na(unified_vcf_and_clonalCLONAL@fix[,1]),]
    unified_vcf_and_clonalSUBCLONAL <- unified_vcf_and_clonalSUBCLONAL[!is.na(unified_vcf_and_clonalSUBCLONAL@fix[,1]),]
    
    vcfR::write.vcf(unified_vcf_and_clonalCLONAL,
                    file = filenameclonal)
    vcfR::write.vcf(unified_vcf_and_clonalSUBCLONAL,
                    file = filenamesubclonal)
    
    readclonal <- vcfR::read.vcfR(filenameclonal)
    # readclonal[80:82,]@fix
    
    contexts_snv_listCLONAL[[vcf_snv]] <- extractSigsSnv(vcf.file=filenameclonal, output='contexts',
                                                         ref.genome=BSgenome.Hsapiens.UCSC.hg19, verbose = T)
    contexts_snv_listSUBCLONAL[[vcf_snv]] <- extractSigsSnv(vcf.file=filenamesubclonal, output='contexts',
                                                            ref.genome=BSgenome.Hsapiens.UCSC.hg19, verbose = T)
    
    system(paste0("rm ", filenameclonal)) ## heavy
    system(paste0("rm ", filenamesubclonal)) ## heavy
    rm(contexts_snv)
  }else{
    cat('File ', vcf_snv, ' does not exist\n')
  }
}


## add counts to ROO
contexts_snv_listCLONAL <- do.call('cbind', contexts_snv_listCLONAL)
contexts_snv_listSUBCLONAL <- do.call('cbind', contexts_snv_listSUBCLONAL)

## Save 96 contexts
objects_sigs_per_CT <- new("exposures_cancertype",
                           cancer_type=opt$cancer_type,
                           type_classification = "cluster1_vs_other",
                           number_categories = 2,
                           id_categories = c('cluster1', 'other'),
                           active_signatures = "character", ## active signatures for this cancer type
                           count_matrices_all = list(contexts_snv_listCLONAL, contexts_snv_listSUBCLONAL),
                           count_matrices_active = list(contexts_snv_listCLONAL, contexts_snv_listSUBCLONAL),
                           sample_names = colnames(contexts_snv_listCLONAL),
                           modification = "none",
                           is_null_active = T,
                           is_empty="Non-empty")
saveRDS(object = objects_sigs_per_CT, file = paste0("../data/roo/", opt$cancer_type, "_nucleotidesubstitution3MSE_ROO2.RDS"))
saveRDS(object = objects_sigs_per_CT, file = paste0("../data/roo/", opt$cancer_type, "_nucleotidesubstitution3MSE_ROO.RDS"))

debug <- F
if(debug){
  
  opt$cancer_type <- 'Kidney-RCC.clearcell'
  objects_sigs_per_CT <- readRDS("../data/roo/Kidney-RCC.clearcell_nucleotidesubstitution3MSE_ROO2.RDS")
  contexts_snv_listCLONAL <- objects_sigs_per_CT@count_matrices_all[[1]]
  contexts_snv_listSUBCLONAL <- objects_sigs_per_CT@count_matrices_all[[2]]
}

sigs_snvCLONAL <- fitToSignatures(
  mut.context.counts=t(contexts_snv_listCLONAL), 
  signature.profiles=SBS_SIGNATURE_PROFILES_V3
)
sigs_snvSUBCLONAL <- fitToSignatures(
  mut.context.counts=t(contexts_snv_listSUBCLONAL), 
  signature.profiles=SBS_SIGNATURE_PROFILES_V3
)

active_sigs <- read.table("../data/cosmic/active_signatures_PCAWGpaper.txt")
ct_changed <- strsplit(toupper(opt$cancer_type), '[.]')[[1]][1]
active_sigs_vec <- names(active_sigs[ct_changed,-1])[active_sigs[ct_changed,-1] == 1]
remove_na <- function(i) i[!is.na(i)]
SBS_SIGNATURE_PROFILES_V3_subset <- SBS_SIGNATURE_PROFILES_V3[,remove_na(match(active_sigs_vec, colnames(SBS_SIGNATURE_PROFILES_V3)))]

sigs_snvCLONALsubset <- fitToSignatures(
  mut.context.counts=t(contexts_snv_listCLONAL), 
  signature.profiles=SBS_SIGNATURE_PROFILES_V3_subset
)
sigs_snvSUBCLONALsubset <- fitToSignatures(
  mut.context.counts=t(contexts_snv_listSUBCLONAL), 
  signature.profiles=SBS_SIGNATURE_PROFILES_V3_subset
)


## Save signature exposures

objects_sigs_per_CT <- new("exposures_cancertype",
                           cancer_type=opt$cancer_type,
                           type_classification = "cluster1_vs_other",
                           number_categories = 2,
                           id_categories = c('cluster1', 'other'),
                           active_signatures = "character", ## active signatures for this cancer type
                           count_matrices_all = list(sigs_snvCLONAL, sigs_snvSUBCLONAL),
                           count_matrices_active = list(sigs_snvCLONALsubset, sigs_snvSUBCLONALsubset),
                           sample_names = rownames(sigs_snvCLONALsubset),
                           modification = "none",
                           is_null_active = F,
                           is_empty="Non-empty")

saveRDS(object = objects_sigs_per_CT, file = paste0("../data/roo/", opt$cancer_type, "_signaturesMSE_ROO2.RDS"))
saveRDS(object = objects_sigs_per_CT, file = paste0("../data/roo/", opt$cancer_type, "_signaturesMSE_ROO.RDS"))

# # 0a6be23a-d5a0-4e95-ada2-a61b2b5d9485	Prost-AdenoCA	TRUE	TRUE
# 
# my_ct <- load_PCAWG(ct = "Prost-AdenoCA", typedata = "signatures", path_to_data = "../data/")
# my_ct <- colSums(my_ct$Y[(grepl("0a6be23a-d5a0-4e95-ada2-a61b2b5d9485", rownames(my_ct$Y))),])
# 
# comparison <- cbind(my_ct, sigs_snv[match(names(my_ct), names(sigs_snv))])
# plot(comparison)
# 
# sigs_snv
# 
# my_ct

