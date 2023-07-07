# install.packages("devtools")
# devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')
library(mutSigExtractor)
library(BSgenome)

#' https://github.com/UMCUGenetics/mutSigExtractor

##' Problem: the vcf from PCAWG is from their own aasembly hs37d5.
##' Looking at genome.fa, it looks as though it's GRCh37 (i.e. hg19)
##' mutSigExtractor uses as a default genome hg19, or we can also use hg38, but that throws an error

library('BSgenome.Hsapiens.UCSC.hg19')
library('BSgenome.Hsapiens.UCSC.hg38')

# hs37d5 <- BSgenome::

file <- ("/Users/morril01/Documents/PhD/CDA_in_Cancer/data/genome.fa.gz")
fasta.seqlengths(file)
  
vcf_snv <- "../../../data/restricted/pcawg/pcawg_restricted_snv/0a6be23a-d5a0-4e95-ada2-a61b2b5d9485.consensus.20160830.somatic.snv_mnv.vcf"

contexts_snv_hg38 <- extractSigsSnv(vcf.file=vcf_snv, output='contexts',
                               ref.genome=BSgenome.Hsapiens.UCSC.hg38)
contexts_snv <- extractSigsSnv(vcf.file=vcf_snv, output='contexts',
                               ref.genome=BSgenome.Hsapiens.UCSC.hg19)
head(contexts_snv)

sigs_snv <- fitToSignatures(
  mut.context.counts=contexts_snv[,1], 
  signature.profiles=SBS_SIGNATURE_PROFILES_V3
)
head(sigs_snv)

# 0a6be23a-d5a0-4e95-ada2-a61b2b5d9485	Prost-AdenoCA	TRUE	TRUE

my_ct <- load_PCAWG(ct = "Prost-AdenoCA", typedata = "signatures", path_to_data = "../../../data/")
my_ct <- colSums(my_ct$Y[(grepl("0a6be23a-d5a0-4e95-ada2-a61b2b5d9485", rownames(my_ct$Y))),])

comparison <- cbind(my_ct, sigs_snv[match(names(my_ct), names(sigs_snv))])
plot(comparison)

sum(my_ct)
sum(sigs_snv)

