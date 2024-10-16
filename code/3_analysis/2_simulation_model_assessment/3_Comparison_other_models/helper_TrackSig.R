# BSgenome.Hsapiens.UCSC.hg19$chr1[13003:13005]
ex_trinucleotides=list("AAA" = 13004,
                       "CAA" = 130092,
                       "GAA" = 13019,
                       "TAA" = 10002,
                       "ACA" = 130099,
                       "CCA" = 130091,
                       "GCA" = 13027,
                       "TCA" = 20032,
                       "AGA" = 13018,
                       "CGA" = 40024,
                       "GGA" = 13016,
                       "TGA" = 13030,
                       "ATA" = 13050,
                       "CTA" = 14029,
                       "GTA" = 130024,
                       "TTA" = 130060,
                       "AAC" = 13020,
                       "CAC" = 130034,
                       "GAC" = 130098,
                       "TAC" = 10537,
                       "ACC" = 13021,
                       "CCC" = 130105,
                       "GCC" = 130104,
                       "TCC" = 20001,
                       "AGC" = 130209,
                       "CGC" = 20128,
                       "GGC" = 13007,
                       "TGC" = 13026,
                       "ATC" = 20031,
                       "CTC" = 13009,
                       "GTC" = 20105,
                       "TTC"  = 20079,
                       "AAG" = 13032,
                       "CAG" = 130204,
                       "GAG" = 13017,
                       "TAG" = 20118,
                       "ACG" = 10563,
                       "CCG" = 130234,
                       "GCG"  = 10579,
                       "TCG" = 10577,
                       "AGG" = 13006,
                       "CGG" = 10580,
                       "GGG"  = 13076,
                       "TGG" = 13012,
                       "ATG" = 13029,
                       "CTG" = 13011,
                       "GTG" = 13014,
                       "TTG" = 130217,
                       "AAT" = 41017,
                       "CAT" = 13028,
                       "GAT" = 40015,
                       "TAT" = 41009,
                       "ACT" = 41024,
                       "CCT" = 13022,
                       "GCT" = 13008,
                       "TCT" = 13010,
                       "AGT" = 130221,
                       "CGT" = 10621,
                       "GGT" = 13013,
                       "TGT" = 13024,
                       "ATT" = 20077,
                       "CTT" = 10529,
                       "GTT" = 130216,
                       "TTT" = 20078)

# which(!sapply(names(ex_trinucleotides), function(i){
#   as.character(BSgenome.Hsapiens.UCSC.hg19$chr1[(ex_trinucleotides[[i]]-1):(ex_trinucleotides[[i]]+1)]) == i
# }))

stopifnot(all(sapply(names(ex_trinucleotides), function(i){
  as.character(BSgenome.Hsapiens.UCSC.hg19$chr1[(ex_trinucleotides[[i]]-1):(ex_trinucleotides[[i]]+1)]) == i
})))

for(i in 1:length(ex_trinucleotides)){
  ex_trinucleotides[[i]] <- c('chr1', ex_trinucleotides[[i]])
}

grab_ref <- function(o){
  # e.g. o=A[C>T]G
  substr(o, 3, 3)
}
grab_alt <- function(o){
  # e.g. o=A[C>T]G
  substr(o, 5, 5)
}

get_trinucleotide_context <- function(o){
  paste0(substr(o, 1, 1), substr(o, 3, 3), substr(o, 7, 7))
}


TrackSig_TrackSigmod <- function (vcfFile, activeInSample, cnaFile = NULL, purity = 1, 
          sampleID = NULL, referenceSignatures = TrackSig:::alex_merged, 
          scoreMethod = "SigFreq", binSize = 100, nCutoff = 10000, 
          desiredMinSegLen = 1, refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){
  
  #############################################################################
  ## change: replacing TrackSig:::vcfToCounts() by TrackSig_vcfToCountsmod() ##
  #############################################################################
  
  if (missing(activeInSample)) {
    assertthat::assert_that(scoreMethod == "Frequency", msg = "When scoreMethod is not equal to \"Frequency\", activeInSample must be provided.")
  }
  else {
    assertthat::assert_that(all(activeInSample %in% colnames(referenceSignatures)))
  }
  assertthat::assert_that(is.numeric(purity) & (0 < purity) & 
                            (purity <= 1), msg = "Purity should be a proportion between 0 and 1\n")
  if (is.null(sampleID)) {
    if (grepl(".txt$", vcfFile)) {
      sampleID <- strsplit(unlist(strsplit(vcfFile, "/"))[length(strsplit(vcfFile, 
                                                                          "/")[[1]])], ".txt")[[1]]
    }
    if (grepl(".vcf$", vcfFile)) {
      sampleID <- strsplit(unlist(strsplit(vcfFile, "/"))[length(strsplit(vcfFile, 
                                                                          "/")[[1]])], ".vcf")[[1]]
    }
    else {
      stop("Failed setting sampleID. Please check input vcf file.")
    }
  }
  context <- TrackSig:::generateContext(c("CG", "TA"))
  vcaf <- countsPerBin <- NULL
  print('here')
  list[vcaf, countsPerBin] <- TrackSig_vcfToCountsmod(vcfFile = vcfFile, 
                                          cnaFile = cnaFile, purity = purity, binSize = binSize, 
                                          nCutoff = nCutoff, context = context, refGenome = refGenome)
  assertthat::assert_that(all(rownames(countsPerBin) %in% rownames(referenceSignatures)), 
                          msg = "Mutation type counts failed.")
  print('here2')
  countsPerBin <- countsPerBin[rownames(referenceSignatures), 
                               , drop = FALSE]
  referenceSignatures <- referenceSignatures[activeInSample]
  if (any(rowSums(countsPerBin)[rowSums(referenceSignatures) == 
                                0] != 0)) {
    print(sprintf("Error in sample %s: Some mutation types have probability 0 under the model, but their count is non-zero. This count vector is impossible under the model.", 
                  sampleID))
  }
  mixtures <- changepoints <- NULL
  list[mixtures, changepoints] <- getChangepointsPELT(vcaf = vcaf, 
                                                      countsPerBin = countsPerBin, referenceSignatures = referenceSignatures, 
                                                      scoreMethod = scoreMethod, binSize = binSize, desiredMinSegLen = desiredMinSegLen)
  if (is.null(changepoints)) {
    vcaf$clust = 1
  }
  else {
    clustIdx = rep(1:(length(changepoints) + 1), times = c(changepoints, 
                                                           max(vcaf$bin)) - c(0, changepoints))
    vcaf$clust = clustIdx[vcaf$bin]
  }
  return(list(mixtures = mixtures, changepoints = changepoints, 
              sampleID = sampleID, binData = vcaf))
}

TrackSig_vcfToCountsmod <- function (vcfFile, cnaFile = NULL, purity = 1, binSize = 100, 
          context = generateContext(c("CG", "TA")), nCutoff = 10000, 
          verbose = F, refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){
  ###############################################################################
  ## NO: replacing TrackSig:::getTrinuc() by TrackSig_getTrinuc_mod()      ##
  ## NO: replacing TrackSig:::parseVcfFile() by TrackSig_parseVcfFilemod() ##
  ###############################################################################
  
  print('here (TrackSig_vcfToCountsmod)')
  vcfFile <- path.expand(vcfFile)
  stopifnot(file.exists(vcfFile))
  if (!is.null(cnaFile)) {
    cnaFile <- path.expand(cnaFile)
    stopifnot(file.exists(cnaFile))
  }
  print('here 2 (TrackSig_vcfToCountsmod)')
  # vcf <- TrackSig:::parseVcfFile(vcfFile, nCutoff, refGenome)
  vcf <- TrackSig_parseVcfFilemod(vcfFile, nCutoff, refGenome)
  # print(head(vcf))
  print('here 3 (TrackSig_vcfToCountsmod)')
  cna <- TrackSig:::parseCnaFile(cnaFile)
  print('here 4 (TrackSig_vcfToCountsmod)')
  vcaf <- TrackSig:::getVcaf(vcf, purity, cna, refGenome, verbose)
  print('here 5 (TrackSig_vcfToCountsmod)')
  
  vcaf <- TrackSig:::getTrinuc(vcaf, refGenome, verbose)
  # vcaf <- TrackSig_getTrinuc_mod(vcaf, refGenome, verbose)
  # print(head(vcaf))
  countsPerBin <- NULL
  list[vcaf, countsPerBin] <- TrackSig:::getBinCounts(vcaf, binSize, context, 
                                           verbose)
  vcaf <- vcaf[, c("chr", "pos", "cn", "mutType", "alt", "phi", 
                   "qi", "bin")]
  rownames(vcaf) <- NULL
  return(list(vcaf = vcaf, countsPerBin = countsPerBin))
}

TrackSig_getTrinuc_mod = function (vcaf, refGenome, verbose = F) {
  if (verbose) {
    print("Making mutation types...")
  }
  print(head(vcaf))
  # assertthat::assert_that(class(refGenome) == "BSgenome")
  # mutRanges <- GenomicRanges::GRanges(paste0("chr", vcaf$chr, 
  #                                            ":", vcaf$pos - 1, "-", vcaf$pos + 1, ":+"))
  # triNuc <- Biostrings::getSeq(refGenome, mutRanges)
  # vcaf$mutType <- as.character(triNuc)
  # mismatchedRef <- which(!(vcaf$ref == substr(vcaf$mutType, 
  #                                             2, 2)))
  # if (length(mismatchedRef) > 0) {
  #   warning(sprintf("%s (of %s) mutations have vcf refrence allele mismatch with the selected reference genome\n", 
  #                   length(mismatchedRef), dim(vcaf)[1]))
  #   substr(vcaf$mutType, 2, 2) <- vcaf$ref
  # }
  # rmSet <- !sapply(triNuc, FUN = Biostrings::hasOnlyBaseLetters)
  # if (sum(rmSet) > 0) {
  #   warning(sprintf("%s (of %s) mutations dropped for uncertain identity in reference genome\n", 
  #                   sum(rmSet), dim(vcaf)[1]))
  #   vcaf <- vcaf[!rmSet, ]
  # }
  # complementSel <- (vcaf$ref == "G" | vcaf$ref == "A")
  # vcaf$mutType[complementSel] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(vcaf$mutType))[complementSel])
  # vcaf$alt[vcaf$ref == "G"] <- as.character(Biostrings::complement(Biostrings::DNAStringSet(vcaf$alt[vcaf$ref == 
  #                                                                                                      "G"])))
  # vcaf$alt[vcaf$ref == "A"] <- as.character(Biostrings::complement(Biostrings::DNAStringSet(vcaf$alt[vcaf$ref == 
  #                                                                                                      "A"])))
  # vcaf$ref[vcaf$ref == "G"] <- "C"
  # vcaf$ref[vcaf$ref == "A"] <- "T"
  return(vcaf)
}

# TrackSig_parseVcfFilemod <- function (vcfFile, cutoff = 10000, refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) 
# {
#  read.table(vcfFile, sep = '\t') 
#   
# }

TrackSig_parseVcfFilemod <- function (vcfFile, cutoff = 10000, refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) 
{
  require(VariantAnnotation)
  print('In TrackSig_parseVcfFilemod')
  vcf <- VariantAnnotation::readVcf(vcfFile, genome = GenomeInfoDb::bsgenomeName(refGenome))
  print(dim(vcf))
  print(VariantAnnotation::info(vcf))
  GenomeInfoDb::seqlevelsStyle(vcf) <- "UCSC"
  print(VariantAnnotation::info(vcf)$t_alt_count)
  # vcf <- vcf[!is.na(VariantAnnotation::info(vcf)$t_alt_count) & 
  #              !is.na(VariantAnnotation::info(vcf)$t_ref_count)]
  print(is.na(VariantAnnotation::info(vcf)$t_alt_count))
  print(dim(vcf))
  assertthat::assert_that(dim(vcf)[1] > 0, msg = "No variants with sufficient t_alt_count or t_ref_count info in VCF")
  assertthat::assert_that("t_alt_count" %in% rownames(VariantAnnotation::info(VariantAnnotation::header(vcf))), 
                          msg = "Tumor alternate variant count \"t_alt_count\" was not found in the vcf header. Please check formatting\n")
  assertthat::assert_that("t_ref_count" %in% rownames(VariantAnnotation::info(VariantAnnotation::header(vcf))), 
                          msg = "Tumor refrence variant count \"t_ref_count\" was not found in the vcf header. Please check formatting\n")
  if (dim(vcf)[1] > cutoff) {
    vcf <- sample(vcf, cutoff)
  }
  if (dim(vcf)[1] <= 400) {
    warning("For high resolution segmentation, recommend >400 mutations total\n")
  }
  return(vcf)
}
