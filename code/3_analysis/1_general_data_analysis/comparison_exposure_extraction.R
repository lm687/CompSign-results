## Comparison of exposure extraction methods

#' The comparison to values obtained by considering mutSigExtractor [19], as well as to values reported
#' in [20] and extracted using sigProfiler and without splitting the mutations into two groups, yield
#' comparable results.

##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
# source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
# source("../3_analysis/recovery_COSMIC_signatures/recover_COSMIC_signatures.R")

library(TMB)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(dplyr)
library(jcolors)
library(viridis)
# library(mutSigExtractor)
library(reshape2)
theme_set(theme_bw())

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../../data/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
enough_samples
# nonexogenous = read.table("../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t",
#                           comment.char = "#", fill = F)
# nonexogenouswSBS1SBS5 = read.table("../../data/cosmic/exogenous_signatures_SBS_withSBS1SBS5.txt", sep = "\t",
#                                    comment.char = "#", fill = F)
# sigs_cosmic0 <- read.table(paste0( "../../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
#                            stringsAsFactors = FALSE, sep = ',', header = TRUE)
# rownames(sigs_cosmic0) <- paste0(substr(sigs_cosmic0$SubType, 1, 1),'[',
#                                  sigs_cosmic0$Type, ']', substr(sigs_cosmic0$SubType, 3, 3))
# sigs_cosmic0 <- sigs_cosmic0[-c(1,2)];
# sigs_cosmic <- colnames(sigs_cosmic0)
# 
# nucleotide_colours_logR <- c('C$>$A/T$>$G'= '#3cb371', 'C$>$G/T$>$G'= '#90ee90', 'C$>$T/T$>$G'= '#66cdaa',
#                              'T$>$A/T$>$G'= '#cd5c5c', 'T$>$C/T$>$G'= '#f4a460')
# nucleotide_colours <- c('C>A' = '#3cb371', 'C>G'= '#90ee90', 'C>T'= '#66cdaa',
#                         'T>A'= '#cd5c5c', 'T>C'= '#f4a460', 'T>G'='red')
# nucleotide_colours_dollar <- c('C$>$A' = '#3cb371', 'C$>$G'= '#90ee90', 'C$>$T'= '#66cdaa',
#                                'T$>$A'= '#cd5c5c', 'T$>$C'= '#f4a460', 'T$>$G'='red')
# 
# pcawg_meta <- read.table("../../data/restricted/pcawg/repository_1567600367.tsv", sep = "\t", head=T)
# metastatic_samples <- pcawg_meta[grepl('Metastatic', pcawg_meta$Specimen.Type),]
# pcawg_meta2 <- read.table("../../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt", sep = "\t", head=T)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
read_info <- function(ct){
  .x <- list(#dataset_all_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../data/", load_all_sigs = T, override_warning_X_Z = T),
             dataset_active_sigs_SP = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../data/", load_all_sigs = F, override_warning_X_Z = T),
             dataset_active_sigs_MSE = load_PCAWG(ct = ct, typedata = "signaturesMSE", path_to_data = "../../../data/", load_all_sigs = F, override_warning_X_Z = T)
             )
  # SaA
   return(.x)
}
##-----------------------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------------------##
read_info_list <- lapply(enough_samples, read_info)
names(read_info_list) <- enough_samples
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## MSE vs QP with SP active signatures
# idx <- 1
# plot(as.vector(read_info_list[[idx]]$dataset_active_sigs_SP$Y),
# as.vector(read_info_list[[idx]]$dataset_active_sigs_MSE$Y))
# 
# MSE_QP_df <- lapply(1:length(read_info_list), function(idx) {
#   cat(idx, '\n')
#   cbind(as.vector(read_info_list[[idx]]$dataset_active_sigs_SP$Y),
# as.vector(read_info_list[[idx]]$dataset_active_sigs_MSE$Y))
#   })
lapply(read_info_list, function(j) sapply(lapply(j, `[[`, 'Y'), ncol))
lapply(read_info_list, function(j) sapply(lapply(j, `[[`, 'Y'), nrow))


## if a signature is inactive in one of the two, add zeros
for(idx in 1:length(read_info_list)){
  mtch1 <- match(colnames(read_info_list[[idx]]$dataset_active_sigs_SP$Y), colnames(read_info_list[[idx]]$dataset_active_sigs_MSE$Y))

  if(!all(!is.na(mtch1))){
    read_info_list[[idx]]$dataset_active_sigs_MSE$Y <- read_info_list[[idx]]$dataset_active_sigs_MSE$Y[,mtch1]
    read_info_list[[idx]]$dataset_active_sigs_MSE$Y[is.na(read_info_list[[idx]]$dataset_active_sigs_MSE$Y)] <- 0
    colnames(read_info_list[[idx]]$dataset_active_sigs_MSE$Y) <- colnames(read_info_list[[idx]]$dataset_active_sigs_SP$Y)
  }

  mtch2 <- match(colnames(read_info_list[[idx]]$dataset_active_sigs_MSE$Y), colnames(read_info_list[[idx]]$dataset_active_sigs_SP$Y))
  
  if(!all(!is.na(mtch2))){
    read_info_list[[idx]]$dataset_active_sigs_SP$Y <- read_info_list[[idx]]$dataset_active_sigs_SP$Y[,mtch2]
    read_info_list[[idx]]$dataset_active_sigs_SP$Y[is.na(read_info_list[[idx]]$dataset_active_sigs_SP$Y)] <- 0
    colnames(read_info_list[[idx]]$dataset_active_sigs_SP$Y) <- colnames(read_info_list[[idx]]$dataset_active_sigs_MSE$Y)
    
  }
  
  stopifnot(colnames(read_info_list[[idx]]$dataset_active_sigs_SP$Y) == colnames(read_info_list[[idx]]$dataset_active_sigs_SP$Y))
  
  intersect_rownames <- intersect(rownames(read_info_list[[idx]]$dataset_active_sigs_SP$Y),  
                                  gsub("[.].*", "", rownames(read_info_list[[idx]]$dataset_active_sigs_MSE$Y)))
  read_info_list[[idx]]$dataset_active_sigs_SP$Y <- read_info_list[[idx]]$dataset_active_sigs_SP$Y[rownames(read_info_list[[idx]]$dataset_active_sigs_SP$Y) %in% intersect_rownames,]
  read_info_list[[idx]]$dataset_active_sigs_MSE$Y <- read_info_list[[idx]]$dataset_active_sigs_MSE$Y[gsub("[.].*", "", rownames(read_info_list[[idx]]$dataset_active_sigs_MSE$Y)) %in% intersect_rownames,]
  stopifnot(nrow(read_info_list[[idx]]$dataset_active_sigs_SP$Y) == nrow(read_info_list[[idx]]$dataset_active_sigs_MSE$Y))
}


read_info_list_df <- melt(lapply(read_info_list, function(j) {
  .x <- sapply(lapply(j, `[[`, 'Y'), as.vector)
  cbind.data.frame(.x, idx=1:nrow(.x))
}), id.vars='idx')
head(read_info_list_df)
read_info_list_df <- (dcast(read_info_list_df, L1+idx~variable, value.var = 'value'))
comparison1_plot <- ggplot(read_info_list_df, aes(x=dataset_active_sigs_SP+1, y=dataset_active_sigs_MSE+1))+
  geom_point()+labs(x='Exposures from QP (group-specific)', y='Exposures from MSE (group-specific)')+
  scale_x_log10()+scale_y_log10()+geom_abline(slope = 1, intercept = 0, lty='dashed', col='blue')

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Exposures from SP directly downloaded from Synapse
SP_PCAWG <- read.table("../../../data/pcawg/PCAWG_sigProfiler_SBS_signatures_in_samples.csv", sep = ',', header = T)
SP_PCAWG <- (melt(SP_PCAWG, id.vars=c('Cancer.Types', 'Sample.Names', 'Accuracy')))
SP_PCAWG <- SP_PCAWG %>% group_by( Cancer.Types, variable ) %>% summarise(sum_value=sum(value))
## These data only have the SP id, and not the patient id
SP_PCAWG

SP_QP_summary <- melt(lapply(read_info_list, function(j) {
  .x <- colSums(Reduce("+", split_matrix_in_half(j$dataset_active_sigs_SP$Y)))
  data.frame(sig=names(.x), value=.x)
  }), value='value')
head(SP_QP_summary)
SP_QP_summary$variable <- NULL

SP_QP_summary <- SP_QP_summary %>%  rename('Cancer.Types' = L1,
                               'variable' = sig,
                               'sum_value' = value)

head(SP_QP_summary)
head(SP_PCAWG)

comparison2 <- (left_join(SP_QP_summary, SP_PCAWG, by=c('Cancer.Types', 'variable')))
comparison2 <- comparison2 %>% group_by(Cancer.Types) %>% mutate(frac_value.x = sum_value.x/sum(sum_value.x),
                                                  frac_value.y = sum_value.y/sum(sum_value.y))
ggplot(comparison2, aes(x=sum_value.x, y=sum_value.y))+geom_point()
comparison2_plot <- ggplot(comparison2, aes(x=frac_value.x, y=frac_value.y))+geom_point()+
  labs(x='Exposures from SP (sum of group-specific)', y='Exposures from QP (cancer type-specific, downloaded)')+
  geom_abline(slope = 1, intercept = 0, lty='dashed', col='blue')

cowplot::plot_grid(comparison1_plot, comparison2_plot)
ggsave("../../../results/figures_paper/signatureextraction_comparison.pdf", height = 4.5, width = 7.5)
##-----------------------------------------------------------------------------------------------------##

#' # install.packages("devtools")
#' # devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')
#' library(mutSigExtractor)
#' library(BSgenome)
#' 
#' #' https://github.com/UMCUGenetics/mutSigExtractor
#' 
#' ##' The vcf from PCAWG is from their own assembly hs37d5.
#' ##' Looking at genome.fa, it looks as though it's GRCh37 (i.e. hg19)
#' ##' mutSigExtractor uses as a default genome hg19, or we can also use hg38, but that throws an error
#' 
#' library('BSgenome.Hsapiens.UCSC.hg19')
#' library('BSgenome.Hsapiens.UCSC.hg38')
#' 
#' # hs37d5 <- BSgenome::
#' 
#' file <- ("/Users/morril01/Documents/PhD/CDA_in_Cancer/data/genome.fa.gz")
#' fasta.seqlengths(file)
#'   
#' vcf_snv <- "../../../data/restricted/pcawg/pcawg_restricted_snv/0a6be23a-d5a0-4e95-ada2-a61b2b5d9485.consensus.20160830.somatic.snv_mnv.vcf"
#' 
#' contexts_snv_hg38 <- extractSigsSnv(vcf.file=vcf_snv, output='contexts',
#'                                ref.genome=BSgenome.Hsapiens.UCSC.hg38)
#' contexts_snv <- extractSigsSnv(vcf.file=vcf_snv, output='contexts',
#'                                ref.genome=BSgenome.Hsapiens.UCSC.hg19)
#' head(contexts_snv)
#' 
#' sigs_snv <- fitToSignatures(
#'   mut.context.counts=contexts_snv[,1], 
#'   signature.profiles=SBS_SIGNATURE_PROFILES_V3
#' )
#' head(sigs_snv)
#' 
#' # 0a6be23a-d5a0-4e95-ada2-a61b2b5d9485	Prost-AdenoCA	TRUE	TRUE
#' 
#' my_ct <- load_PCAWG(ct = "Prost-AdenoCA", typedata = "signatures", path_to_data = "../../../data/")
#' my_ct <- colSums(my_ct$Y[(grepl("0a6be23a-d5a0-4e95-ada2-a61b2b5d9485", rownames(my_ct$Y))),])
#' 
#' comparison <- cbind(my_ct, sigs_snv[match(names(my_ct), names(sigs_snv))])
#' plot(comparison)
#' 
#' sum(my_ct)
#' sum(sigs_snv)

