
rm(list = ls())

library(TMB)
library(optparse)
library(CompSign)
library(BSgenome.Hsapiens.UCSC.hg19) ## for TrackSig
library(parallel)

source("3_analysis/2_simulation_model_assessment/3_Comparison_other_models/helper_TrackSig.R") ## functions and loading ex_trinucleotides: locations of each trinucleotide in hg19, as an example

debug <- F
if(debug){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("../../../") ## from where snakemake is run
  
  opt <- list()
  opt$input = ("../data/assessing_models_simulation/datasets/multiple_GenerationMixturewithCCFb_68_NA_100_96_0.05_NA_NA_NA/multiple_GenerationMixturewithCCFb_68_NA_100_96_0.05_NA_NA_NA_dataset1.RDS")
  opt$model <- 'SigFreq'
  opt$datasetgeneration = 'GenerationMixturewithCCFb'
}
option_list = list(
  make_option(c("--model"), type="character", default='SigFreq',
              help="Which model to use for inference", metavar="character"),
  make_option(c("--datasetgeneration"), type="character", default=NULL,
              help="How was the data generated? This is used to get the number of active signatures", metavar="character"),
  make_option(c("--input"), type="character", default=NA,
              help="Input file with dataset (RDS)", metavar="character"),
  make_option(c("--output"), type="character", default=NA,
              help="Output file in which to write the results of the inference (RDS file)", metavar="character"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Get active signatures for the cancer type we are simulating from

if(opt$datasetgeneration == 'GenerationMixturewithCCFa'){
  ## Lymph_BNHL
  ct <- "Lymph-BNHL"
}else if(opt$datasetgeneration == 'GenerationMixturewithCCFb'){
  ct <- "Lung-SCC"
}else{
  stop()
}

cancer_type_active_sigs <- load_PCAWG(ct, typedata = 'signaturesPCAWG', path_to_data = "../data/", override_warning_X_Z = T)
active_sigs <- colnames(cancer_type_active_sigs$Y)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
##' Including all signatures, even with artefacts (SBS_SIGNATURE_PROFILES_V2 and SBS_SIGNATURE_PROFILES_V3 from
##' mutSigExtractor do not include artefactual signatures)
sigs_cosmic0 <- read.table(paste0( "../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
                           stringsAsFactors = FALSE, sep = ',', header = TRUE)
rownames(sigs_cosmic0) <- paste0(substr(sigs_cosmic0$SubType, 1, 1),'[',
                                 sigs_cosmic0$Type, ']', substr(sigs_cosmic0$SubType, 3, 3))
sigs_cosmic0 <- sigs_cosmic0[-c(1,2)];
sigs_cosmic <- colnames(sigs_cosmic0)
## rownames should be in this format:
# C_A_GCT 0.00584 2.04e-03 0.00345 0.00421 1.14e-03 6.73e-03 7.38e-04 1.57e-02 1.95e-03 0.000350 4.99e-01 1.33e-04
# C_A_TCA 0.03290 1.90e-02 0.03570 0.46400 2.67e-03 6.69e-04 1.74e-03 2.19e-01 1.43e-02 0.007370 2.37e-05 2.02e-03

give_tracksig_mutation_format <- function(o){
  paste0(substr(o, 3, 3), '_', substr(o, 5, 5), '_', substr(o, 1, 1), substr(o, 3, 3), substr(o, 7, 7))
}
rownames(sigs_cosmic0) <- sapply(rownames(sigs_cosmic0), give_tracksig_mutation_format)
##-----------------------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------------------##
system(paste0('mkdir -p ', dirname(opt$output)))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
input_obj = readRDS(opt$input)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## split in patients
sample_names <- names(input_obj$simulated_mutations)
# input_obj$exposures_bins_QP_list$Y

## in the input object input_obj, (chromosome  position) are just placeholders and do not correspond to the genomic change
## of the mutation. E.g. in position 21 of 24861963 there might not even be a C, let alone a ACG context.
# chromosome  position      ccf mutation
#          21  24861963 1.023263  A[C>T]G
#           8  90911235 1.023263  T[C>G]T
# for this, we use the genome used in TrackSig
# refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
## and for each trinucleotide context (e.g. GCA for G[C>A]A) get a position
## this is done manually. Check object ex_trinucleotides, which is created in:
# 3_analysis/2_simulation_model_assessment/3_Comparison_other_models/helper_TrackSig.R

scaling_factor <- 100000

temp_files_list <- sapply(sample_names, function(sample_it){
  ## prepare the file for TrackSig, for each sample independently
  all_muts = do.call('rbind', input_obj$simulated_mutations[[sample_it]])
  all_muts$trinucleotide_context <- sapply(all_muts$mutation, get_trinucleotide_context)
  matched_trinucleotides_contexts <- do.call('rbind', ex_trinucleotides[match(all_muts$trinucleotide_context, names(ex_trinucleotides))])
  colnames(matched_trinucleotides_contexts) <- c('chr', 'position')
  
  all_muts <- cbind(all_muts, mut_with_trinucleotide_context=matched_trinucleotides_contexts)
  
  all_muts_to_write = data.frame(all_muts$mut_with_trinucleotide_context.chr,
                                 all_muts$mut_with_trinucleotide_context.position, '.',
             sapply(all_muts$mutation, grab_ref), sapply(all_muts$mutation, grab_alt),
             '.', '.', paste0('t_ref_count=', scaling_factor, ';t_alt_count=', round(all_muts$ccf*scaling_factor)))
  
  uuid_session <- uuid::UUIDgenerate()
  system("mkdir -p tmp/")
  temp_file = paste0('tmp/', uuid_session, sample_it)
  writeLines('##fileformat=VCFv4.1\n##INFO=<ID=t_alt_count,Number=1,Type=Integer,Description="Tumour alt count simulated from ccf">\n##INFO=<ID=t_ref_count,Number=1,Type=Integer,Description="Tumour reference count simulated from ccf">\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO', con = temp_file)
  write.table(all_muts_to_write, file = temp_file, sep = '\t', quote = F, col.names = F, row.names = F, append = T)
  ## create a file for TrackSig. The file is tab-separated and should look like that:
  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
  # chr3	134920639	.	C	G	.	.	t_ref_count=55;t_alt_count=54
  # chr20	11380899	.	C	T	.	.	t_ref_count=40;t_alt_count=56
  # chr2	213943163	.	C	T	.	.	t_ref_count=42;t_alt_count=59
  # chr14	89865291	.	C	T	.	.	t_ref_count=64;t_alt_count=53
  # chr10	77313198	.	C	G	.	.	t_ref_count=45;t_alt_count=60
  # chr20	3706848	.	C	A	.	.	t_ref_count=50;t_alt_count=48
  return(temp_file)
})

# TrackSig:::TrackSig(vcfFile = , activeInSample = active_sigs, purity = , sampleID = , referenceSignatures = , binSize = 100, scoreMethod = )
# TrackSig_TrackSigmod(vcfFile = temp_file, activeInSample = active_sigs, purity = 1,
#                      sampleID = sample_it, referenceSignatures = sigs_cosmic0, binSize = 100,
#                      scoreMethod = opt$model)


##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
time1=Sys.time()
## this has to be run for each sample
traj <- parallel::mclapply(1:length(sample_names), function(it){
  cat('Running TrackSig\n')
  TrackSig::TrackSig(vcfFile = temp_files_list[[it]], activeInSample = active_sigs, purity = 1,
                             sampleID = sample_names[[it]], referenceSignatures = sigs_cosmic0, binSize = 100,
                             scoreMethod = opt$model)
})
time2=Sys.time() 

for(temp_file_it in temp_files_list){
  system(paste0("rm ", temp_file_it))
}

## traj is a list of TrackSig outputs for all samples in a cancer type
saveRDS(object = traj, file = opt$output)
