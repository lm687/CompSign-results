rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(CompSign)
# data(package='CompSign')

simplified_object
## Looking at the differential abundance of mutational signatures in chromosomes
data(ProstAdenoCA_chrom, package='CompSign')

matrices_exposures_active_sigs <- sapply(ProstAdenoCA_chrom$all_exposures_ct_active, t)

chroms <- rownames(matrices_exposures_active_sigs[[1]])

## we will unify chromosomes: if some is absent, we will add exposures of z
matrices_exposures_active_sigs <- sapply(matrices_exposures_active_sigs, function(i){
  rwni <- rownames(i)
  n_missing <- chroms %in% rwni
  if(any(!n_missing)){
    i <- rbind(i, t(sapply(1:sum(!n_missing), function(j) rep(0, ncol(i)))))
    rownames(i) <- c(rwni, chroms[!(chroms %in% rwni)])
  }
  i[chroms,]
}, simplify = F)

do.call('rbind', matrices_exposures_active_sigs)


## random effects:
## there are 24 observations per patient (one for each chromosome). We use patient-specific intercepts

object_for_DA <- list(x=cbind(1, do.call('rbind', lapply(matrices_exposures_active_sigs, function(x) diag(nrow(x))))[,-1]), ## removing the first column: the first chromosome will become the baseline
                      z=do.call('rbind', lapply(1:length(chroms), function(unused) diag(length(matrices_exposures_active_sigs)))),
                      Y=do.call('rbind', matrices_exposures_active_sigs))
sapply(object_for_DA, dim)


## to debug
object_for_DA$x <- cbind(rep(1, nrow(object_for_DA$x)), sample(c(0,1), size = nrow(object_for_DA$x), replace = T))
sapply(object_for_DA, dim)

re_run <- F
if(re_run){
  chrom_diagREDM <- wrapper_run_TMB(model = 'diagRE_DM', object = object_for_DA, smart_init_vals = F)
  saveRDS(chrom_diagREDM, "../../../../data/additional_use_cases/chromosomes.RDS")
}

model = 'diagRE_DM'
object = object_for_DA

diagDM_no_small_sigs <- wrapper_run_TMB(object = simplified_object,
                                        model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagDM_no_small_sigs
