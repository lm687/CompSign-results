rm(list = ls())
setwd("/Users/morril01/Documents/PhD/GlobalDA/code")
source("2_inference_TMB/helper_TMB.R")

library(ggplot2)

##' differential precision: do full model vs precision model, chisqrt with 2 degrees of freedom.
##' I can use the results from the multinomial simulation of PCAWG data for that.
##' 
##' Same for differential abundance
##' 
##' Run bias pipeline in non-DA dataset, to get the distribution of p-values
##' 
##' check distribution if it's chi sqrt or not. simulate under equality (under the null)
##' 
##' check if, for the differential precision analysis, and under the null, the statistics follow
##' a chi-sqrt distribution (which would be the normal test to use in a LRT of full vs restricted
##' model for differential precision)
##' 

## get the results for the no-DA model simulating using the multinomial, using PCAWG

fles <- list.files("../data/assessing_models_simulation/inference_results/TMB/nlminb/", full.names = T)
flesPCAWG <- fles[grepl('PCAWG', fles)]
flesPCAWG <- flesPCAWG[grepl('-999_', flesPCAWG)] ## selecting no DA

readPCAWG <- lapply(flesPCAWG, readRDS)
pvalsPCAWG <- sapply(readPCAWG, wald_TMB_wrapper)
ggplot(data.frame(x=pvalsPCAWG), aes(x=x))+geom_histogram()+theme_bw()+
  ggtitle(paste0('pvals no-DA all ', generation))+labs(x='p-value')
ggsave(paste0("../results/models_explanatory/uniformity_pvals_null/", 'uniformity_pvals_null_all_generations_1', '.pdf'))
    
df_pvalsPCAWG <- data.frame(name=gsub("_.*", "", gsub("multiple_", "", basename(flesPCAWG))),
                            pvalsPCAWG, model=sapply(basename(flesPCAWG), function(i) strsplit(i, '_')[[1]][8]))
pdf(paste0("../results/models_explanatory/uniformity_pvals_null/", 'uniformity_pvals_null_all_generations', '.pdf'),
    width=5, height = 5)
for(n in unique(df_pvalsPCAWG$name)){
  cat(n,'\n')
  df <- df_pvalsPCAWG[(df_pvalsPCAWG$name == n) & !is.na(df_pvalsPCAWG$pvalsPCAWG),]
  if(nrow(df) > 0){
    print(ggplot(df, aes(x=pvalsPCAWG))+geom_histogram()+theme_bw()+
            lims(x=c(0,1))+ggtitle(n)+facet_wrap(.~model))
  }
  rm(df)
}
dev.off()


# for(generation in c('GenerationMixturefewersignaturespairedBoneOsteosarcPCAWG',
#                       'GenerationMixturefewersignaturespairedKidneyRCCPCAWG')){
#   
#   flesgeneration <- flesPCAWG[grepl(generation, flesPCAWG)]
#   flesgenerationdiagREDM <- flesgeneration[grepl('diagREDM_', flesgeneration)]
#   readgeneration <- lapply(flesgeneration, readRDS)
#   readgenerationdiagREDM <- lapply(flesgenerationdiagREDM, readRDS)
# 
#   pvalsNODA_all <- sapply(readgeneration, wald_TMB_wrapper)
#   pvalsNODA <- sapply(readgenerationdiagREDM, wald_TMB_wrapper)
#   
#   hist(pvalsNODA_all)
#   hist(pvalsNODA_all, breaks=100)
#   
#   # ggplot(data.frame(x=pvalsNODA_all), aes(x=x))+geom_histogram()+theme_bw()+
#   #   ggtitle(paste0('pvals no-DA all ', generation))+labs(x='p-value')
#   # ggsave(paste0("../results/models_explanatory/uniformity_pvals_null/", 'uniformity_pvals_null_', generation, '_all.pdf'))
#   
#   pvalsNODA_all
#   
#   df_pvals_all <- data.frame(name=gsub("_.*", "", gsub("multiple_", "", basename(flesgeneration))),
#                         pvalsNODA_all, model=sapply(basename(flesgeneration), function(i) strsplit(i, '_')[[1]][8]))
#   pdf(paste0("../results/models_explanatory/uniformity_pvals_null/", 'uniformity_pvals_null_', generation, '.pdf'),
#       width=3, height = 3)
#   for(n in unique(df_pvals_all$name)){
#     print(ggplot(df_pvals_all[df_pvals_all$name == n,], aes(x=pvalsNODA_all))+geom_histogram()+theme_bw()+
#             lims(x=c(0,1))+ggtitle(n)+facet_wrap(.~model))
#   }
#   dev.off()
# 
# }

LRT_full_reduced <- list.files("../data/assessing_models_simulation/inference_results/TMB_LRT/nlminb/", full.names = T)
LRT_full_reduced <- LRT_full_reduced[!grepl('_NC.RDS', LRT_full_reduced)]
LRT_full_reduced_diagREDM <- LRT_full_reduced[grepl('diagREDM_', LRT_full_reduced)]
LRT_full_reduced_diagREDMsinglelambda <- LRT_full_reduced[grepl('diagREDMsinglelambda_', LRT_full_reduced)]
LRT_full_reduced_diagREDMreducedmodel <- LRT_full_reduced[grepl('diagREDMreducedmodel_', LRT_full_reduced)]

runs_LRT <- list(diagREDM=lapply(LRT_full_reduced_diagREDM, readRDS),
                 diagREDMsinglelambda=lapply(LRT_full_reduced_diagREDMsinglelambda, readRDS),
                 diagREDMreducedmodel=lapply(LRT_full_reduced_diagREDMreducedmodel, readRDS))

runs_LRT_summaries <- lapply(runs_LRT, function(i) sapply(i, function(j) j[2]))
runs_LRT_pdHess <- lapply(runs_LRT_summaries, function(j) sapply(j, `[`, 'pdHess'))
sapply(runs_LRT_pdHess, function(i) sum(unlist(i) == T, na.rm = T))/length(runs_LRT_pdHess[[1]])

## distribution of LRT statistic for the differential abundance test

## distribution of LRT statistic for the differential precision test
## for each dataset, compute the differential precision LRT
all(gsub('diagREDM_', '', LRT_full_reduced_diagREDM) == gsub('diagREDMsinglelambda_', '', LRT_full_reduced_diagREDMsinglelambda))
all(gsub('diagREDM_', '', LRT_full_reduced_diagREDM) == gsub('diagREDMreducedmodel_', '', LRT_full_reduced_diagREDMreducedmodel))

runs_LRT_opt <- lapply(runs_LRT, function(i) sapply(i, function(j) j[1]))

## diagREDM vs diagREDMsinglelambda
lrts <- as.numeric(sapply(runs_LRT_opt$diagREDM, function(i) try(i$objective)))-as.numeric(sapply(runs_LRT_opt$diagREDMsinglelambda, function(i) try(i$objective)))
hist(lrts, breaks = 10000)

## diagREDM vs diagREDMreducedmodel
lrt2 <- as.numeric(sapply(runs_LRT_opt$diagREDM, function(i) try(i$objective)))-as.numeric(sapply(runs_LRT_opt$diagREDMreducedmodel, function(i) try(i$objective)))
hist(lrt2)
                                                                                  

LLs <- data.frame(LLdiagREDM= as.numeric(sapply(runs_LRT_opt$diagREDM, function(i) try(i$objective))),
   LLdiagREDMsinglelambda = as.numeric(sapply(runs_LRT_opt$diagREDMsinglelambda, function(i) try(i$objective))),
   LLdiagREDMreducedmodel = as.numeric(sapply(runs_LRT_opt$diagREDMreducedmodel, function(i) try(i$objective))))

pairs(LLs)

remove_na <- function(i) i[!is.na(i)]

noDA <- LLs[grep('-999', LRT_full_reduced_diagREDM),]
hist(noDA$LLdiagREDM - noDA$LLdiagREDMreducedmodel)
plot(density(remove_na(noDA$LLdiagREDM - noDA$LLdiagREDMreducedmodel)))

plot(density(remove_na(noDA$LLdiagREDM - noDA$LLdiagREDMsinglelambda)))
