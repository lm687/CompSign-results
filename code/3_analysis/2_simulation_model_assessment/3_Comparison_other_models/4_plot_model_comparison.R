##------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(CompSign)
# library(HiLDA)
library(ggplot2)
# library(lsa) ## cosine similarity
library(reshape2)
library(dplyr)
# library(GGally)
library(gridExtra)
theme_set(theme_bw())
source("../../../2_inference_TMB/helper_TMB.R")
source("helper.R")

##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
# all_in_datasets <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/datasets/")
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
names_out_TMB <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/", return_filenames = T)
names_out_HiLDA <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/HiLDA/", remove_HiLDAGlobal = T, return_filenames = T)
names_out_TCSM <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/TCSM/", return_filenames = T)
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
## TMB
all_out_TMB <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/")
# 1283  files found
all_out_TMB <- all_out_TMB[-grep('_NC', names(all_out_TMB))]
names(all_out_TMB) <- sapply(names(all_out_TMB), rename_datasets_fun)

names(all_out_TMB)
pvals <- sapply(all_out_TMB, wald_TMB_wrapper)
max(pvals, na.rm = T)
min(pvals, na.rm = T)

df_pvals <- data.frame(pvals_TMB = pvals, DA_TMB=pvals <= 0.05, 
                       pi_softmax=sapply(as.numeric(sapply(names(pvals), function(i) strsplit(i, '_')[[1]][6])), function(i) softmax(c(i, 0))[1]),
                       pi=as.numeric(sapply(names(pvals), function(i) strsplit(i, '_')[[1]][6])),
                       numsamples=sapply(names(pvals), function(i) strsplit(i, '_')[[1]][2]),
                       numsmuts=sapply(names(pvals), function(i) strsplit(i, '_')[[1]][3]),
                       datasetgeneration=sapply(names(pvals), function(i) strsplit(i, '_')[[1]][1]))
ggplot(df_pvals, aes(x=factor(pi, levels=sort(unique(df_pvals$pi))), fill=factor(DA_TMB)))+geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(.~interaction(numsamples, numsmuts, datasetgeneration))

ggplot(df_pvals, aes(x=factor(round(pi_softmax, 3), levels=sort(unique(round(df_pvals$pi_softmax, 3)))), fill=factor(DA_TMB)))+geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(.~interaction(numsamples, numsmuts, datasetgeneration))

ggplot(df_pvals %>% group_by(pi_softmax, pi, numsamples, numsmuts, datasetgeneration) %>% 
         summarise(frac_DA=mean(DA_TMB, na.rm = T)),
       aes(x=factor(round(pi_softmax, 3), levels=sort(unique(round(df_pvals$pi_softmax, 3)))),
           y = frac_DA, group=interaction( numsamples, numsmuts, datasetgeneration),
           col = interaction( paste0('n=', numsamples), paste0('T=', numsmuts), sep = ', '),
           shape = interaction( paste0('n=', numsamples), paste0('T=', numsmuts), sep = ', ')))+
  geom_vline(xintercept = '0.007', lty=4, alpha=0.4)+
  geom_line()+geom_point()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  facet_wrap(.~interaction(datasetgeneration))+
  labs(x=latex2exp::TeX('$\\pi$'), y='Fraction of DA runs')+
  labs(col='Parameter combination', shape='Parameter combination')+
  theme(legend.position = 'bottom')+
  scale_color_manual(values = c('#004D40', '#C86518', '#3270BB', '#A7BBAD', '#7D0B1B'))
ggsave("../../../../results/figures_paper/comparison_methods/TMB_fraction_DA.pdf", height = 3.5, width = 8)
# geom_rect(data=df_pvals[1,], aes(xmin=-Inf, xmax='0.007', ymin=-Inf, ymax=Inf, fill='red'), alpha=0.02)

df_pvals$datasetnotDA = df_pvals$pi_softmax < 0.01

df_pvals_parameter_expand = expand.grid(apply(df_pvals[,c('numsamples', 'numsmuts', 'datasetgeneration')], 2, unique, simplify = F))

pval_thresholds
give_AUC_df <- function(df_pvals, df_pvals_parameter_expand, column_pvalues_arg, pval_thresholds_arg=pval_thresholds, remove_NA=T, ...){
  aucs_TMB <- apply(df_pvals_parameter_expand, 1, function(parameter_combination){
    .x <- sapply(pval_thresholds_arg, function(pval_threshold_it){
      input_df = df_pvals %>% filter(numsamples == parameter_combination[1],
                          numsmuts == parameter_combination[2],
                          datasetgeneration == parameter_combination[3])
      if(nrow(input_df) == 0){
        return(c(NA, NA))
      }else{
        get_sensivity_specificity(pvals_df = input_df,
                                  model = NULL,
                                  pval_threshold = pval_threshold_it,
                                  not_subset_df = T, column_pvalues = column_pvalues_arg, ...)
      }
      
    })
    colnames(.x) <- pval_thresholds_arg
    .x
  }, simplify = F)
  names(aucs_TMB) <- apply(df_pvals_parameter_expand, 1, paste0, collapse='_')
  aucs_TMB = dcast(melt(aucs_TMB), L1+Var2~Var1, value.var = 'value')
  aucs_TMB$numsamples =  sapply(aucs_TMB$L1, function(i) strsplit(i, '_')[[1]][1])
  aucs_TMB$numsmuts =  sapply(aucs_TMB$L1, function(i) strsplit(i, '_')[[1]][2])
  aucs_TMB$datasetgeneration =  sapply(aucs_TMB$L1, function(i) strsplit(i, '_')[[1]][3])
  if(remove_NA) aucs_TMB <- aucs_TMB[!is.na(aucs_TMB$FPR),]
  return(aucs_TMB)
}

give_AUC_plot <- function(df_AUC){
  ggplot(df_AUC,
         aes(y=TPR, x=FPR,
             col = interaction( paste0('n=', numsamples), paste0('T=', numsmuts), sep = ', '),
             shape = interaction( paste0('n=', numsamples), paste0('T=', numsmuts), sep = ', ')
         ))+geom_line()+geom_point()+
    labs(y='TPR', x='FPR')+
    labs(col='Parameter combination', shape='Parameter combination')+
    theme(legend.position = 'bottom')+
    scale_color_manual(values = c('#004D40', '#C86518', '#3270BB', '#A7BBAD', '#7D0B1B'))+
    facet_wrap(.~datasetgeneration)
}
aucs_TMB <- give_AUC_df(df_pvals, df_pvals_parameter_expand, column_pvalues_arg = 'pvals_TMB')
give_AUC_plot(aucs_TMB %>% filter(datasetgeneration != 'C2'))
ggsave("../../../../results/figures_paper/comparison_methods/TMB_ROC.pdf", height = 3.5, width = 8)

##------------------------------------------------------------------------##
## HiLDA
# all_out_HiLDA <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/HiLDA/",
#                                      remove_HiLDAGlobal = T)
# names(all_out_HiLDA) <- modify_names_HiLDA(names(all_out_HiLDA) )
all_out_HiLDA_betas <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/HiLDA/",
                                     remove_HiLDAGlobal = T, HiLDA_return='betas')
# 800  files found
names(all_out_HiLDA_betas) <- modify_names_HiLDA(names(all_out_HiLDA_betas))

## globaltest for hildaglobal runs
df_pvals_HiLDA <- sapply(all_out_HiLDA, hildaGlobal_from_beta)
df_pvals_HiLDA_alpha <- sapply(all_out_HiLDA, hildaGlobal_from_alpha)

table(df_pvals_HiLDA, df_pvals_HiLDA_alpha)

head(all_out_HiLDA$`TwoCT_100_50_NA_NA_1_HiLDAGlobal_NA_NA_NA_dataset0`$BUGSoutput$summary)

# ggplot(melt(all_out_HiLDA$`TwoCT_100_50_NA_NA_0_HiLDAGlobal_NA_NA_NA_dataset0`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()
# ggplot(melt(all_out_HiLDA$`TwoCT_100_50_NA_NA_1_HiLDAGlobal_NA_NA_NA_dataset0`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()
# ggplot(melt(all_out_HiLDA$`TwoCT_100_50_NA_NA_-1_HiLDAGlobal_NA_NA_NA_dataset0`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()
# ggplot(melt(all_out_HiLDA$`TwoCT_100_50_NA_NA_-8_HiLDAGlobal_NA_NA_NA_dataset0`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()
# ggplot(melt(all_out_HiLDA$`TwoCT_100_50_NA_NA_-8_HiLDAGlobal_NA_NA_NA_dataset1`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()
# ggplot(melt(all_out_HiLDA$`TwoCT_100_50_NA_NA_-8_HiLDAGlobal_NA_NA_NA_dataset2`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()
# ggplot(melt(all_out_HiLDA$`TwoCT_100_50_NA_NA_-8_HiLDAGlobal_NA_NA_NA_dataset3`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()
# 
# ggplot(melt(all_out_HiLDA$`v4_100_50_NA_NA_-8_HiLDA_NA_NA_NA_dataset0`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()
# ggplot(melt(all_out_HiLDA$`v4_100_50_NA_NA_-1_HiLDA_NA_NA_NA_dataset0`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()
# ggplot(melt(all_out_HiLDA$`v4_100_50_NA_NA_1_HiLDA_NA_NA_NA_dataset0`$BUGSoutput$sims.list$alpha[,,1]),
#        aes(x=value, col=Var2, group=Var2))+geom_density()

table(df_pvals_HiLDA)
table(df_pvals_HiLDA_alpha)
df_pvals$pvals_HiLDA <- df_pvals_HiLDA[match(names(pvals), names(df_pvals_HiLDA))]
df_pvals$DA_HiLDA <- (df_pvals$pvals_HiLDA <= 0.05)
aucs_HiLDA <- give_AUC_df(df_pvals, df_pvals_parameter_expand, column_pvalues_arg = 'pvals_HiLDA')

##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
## TCSM
all_out_TCSM <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/TCSM/")
# 1280  files found
names(all_out_TCSM) <- sapply(names(all_out_TCSM), rename_datasets_fun)
all_out_TCSM[[1]]$effect
all_out_TCSM[[1]]$sigma
## one value for each signature
as.numeric(all_out_TCSM[[1]]$EstSignificance$Mann.Whitney.U.p.value)
pvals_TCSM <- sapply(all_out_TCSM, function(j) min(as.numeric(j$EstSignificance$Mann.Whitney.U.p.value)))
df_pvals$pvals_TCSM <- pvals_TCSM[match(names(pvals), names(pvals_TCSM))]
df_pvals$DA_TCSM <- df_pvals$pvals_TCSM  <= 0.05
ggplot(df_pvals, aes(x=factor(round(pi_softmax, 3), levels=sort(unique(round(pi_softmax, 3)))), 
                     fill=factor(DA_TCSM)))+geom_bar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  facet_wrap(.~interaction(numsamples, numsmuts, datasetgeneration))

aucs_TCSM <- give_AUC_df(df_pvals, df_pvals_parameter_expand, column_pvalues_arg = 'pvals_TCSM')

##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
## Comparison

## Comparison of p-values or DA, as mixing fraction increases
plot(pvals[match(names(pvals_TCSM), names(pvals))], pvals_TCSM)
plot(pvals[match(names(pvals_TCSM), names(pvals))], pvals_TCSM)
table(TMB=df_pvals$DA_TMB, TCSM=df_pvals$DA_TCSM)
table(TMB=df_pvals$DA_TMB, TCSM=df_pvals$DA_TCSM, HiLDA=df_pvals$DA_HiLDA)

## FP in TCSM
df_pvals$DA_TMB[grepl('-8', rownames(df_pvals))]
df_pvals$DA_TCSM[grepl('-8', rownames(df_pvals))]

accuracies_all <- lapply(c('pvals_TCSM', 'pvals_HiLDA', 'pvals_TMB'), function(pval_name){
  give_AUC_df(df_pvals, df_pvals_parameter_expand, column_pvalues_arg = pval_name, 
                              give_accuracy=T, pval_thresholds_arg=0.05, remove_NA = FALSE)})
names(accuracies_all) <- c('TCSM', 'HiLDA', 'diagREDM')
accuracies_all_df <- melt(accuracies_all, id.vars=colnames(accuracies_all[[1]]))
accuracies_all_df <- accuracies_all_df %>% arrange(datasetgeneration, as.numeric(numsamples),
                                                   factor(numsmuts, levels=c("50","1773", "Obs" )),
                                                   L1)
accuracies_all_df <- accuracies_all_df[,-c(3,4)]
colnames(accuracies_all_df) <- c('Model', 'Significance level', 'Accuracy', 'Balanced accuracy', 'FPR', 'TPR', 'N', 'T', 'Dataset')
accuracies_all_df <- accuracies_all_df[is.nan(accuracies_all_df$Accuracy) | !is.na(accuracies_all_df$Accuracy),]
print(xtable::xtable(accuracies_all_df),include.rownames=FALSE)
# View(accuracies_all_df)
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
## Compare extracted signatures
## This is only possible between TCSM and TMB, as HiLDA uses Shiraishi signatures

## using all_out_TCSM instead of all_out_TCSM_matched, because we want to show the results for all the runs
## which were successful in TCSM, and not exclude those from HiLDA for which we were not able to get results
## (those of highest mutation number)

## actual signatures
## using v4 and v7, which use the signatures from Lymph-CLL
signature_definitions <- read.table("../../../../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv", sep=',', header = T)
signature_definitions_cats <- paste0(substr(signature_definitions$SubType, 1, 1), '[', signature_definitions$Type, ']', substr(signature_definitions$SubType, 3, 3))
signature_definitions <- apply(signature_definitions[,-c(1,2)], 2, function(i) i)
rownames(signature_definitions) <- signature_definitions_cats
stopifnot(colnames(all_out_TCSM[[1]]$signatures) == gsub("\\[|\\]|>", ".", rownames(signature_definitions)))

true_signatures <- load_PCAWG("Lymph-CLL", typedata = "signaturesPCAWG", simulation = F, path_to_data="../../../../data/", override_warning_X_Z=T)
true_signatures <- signature_definitions[,colnames(true_signatures$Y)]

## cosine similarity for each TCSM signature to the true signatures
## higher: higher similarity
get_set_cosines_from_TCSM <- function(TCSM_object, true_signatures_arg){
  TCSM_object$signatures <- as.matrix(TCSM_object$signatures)
  if(ncol(TCSM_object$signatures) != length(rownames(signature_definitions))){
    TCSM_object$signatures <- TCSM_object$signatures[,match(gsub("\\[|\\]|>", ".", rownames(signature_definitions)), colnames(TCSM_object$signatures))]
    colnames(TCSM_object$signatures) = gsub("\\[|\\]|>", ".", rownames(signature_definitions))
    TCSM_object$signatures[is.na(TCSM_object$signatures)] <- 0
  }
  stopifnot(colnames(TCSM_object$signatures) == gsub("\\[|\\]|>", ".", rownames(signature_definitions)))
  
  all_signature_cosine <- apply(TCSM_object$signatures, 1, function(i){
    sapply(1:ncol(true_signatures_arg), function(j){
      lsa::cosine(x = unlist(i), y = as.vector(true_signatures_arg[,j]))
    })
  })
  
  find_signatures_cosines_from_matrix <- function(all_signature_cosine, true_signatures){
    ## find best set of matches. Start with the highest similarities
    all_signature_cosine0 <- all_signature_cosine
    order_vec_cosine <- order(apply(all_signature_cosine, 2, max), decreasing = T)
    res <- rep(NA, ncol(true_signatures))
    for(i in 1:(ncol(true_signatures)-1)){
      res[i] <- max(all_signature_cosine0[,order_vec_cosine[1]])
      all_signature_cosine0 <- all_signature_cosine0[,-order_vec_cosine[1]]
      if(i != (ncol(true_signatures)-1)){
        order_vec_cosine <- order(apply(all_signature_cosine0, 2, max), decreasing = T)
      }
    }
    res[ncol(true_signatures)] = max(all_signature_cosine0)
    return(res)
  }
  
  return(find_signatures_cosines_from_matrix(all_signature_cosine, true_signatures_arg))
}

cosines_from_TCSM <- lapply(all_out_TCSM, get_set_cosines_from_TCSM, true_signatures)
## find these cosine similarities compared to random sigs
signatures_random <- lapply(1:length(all_out_TCSM), function(i) sample(colnames(signature_definitions), 4, replace = F) )
cosines_from_TCSM_randomsigs <- lapply(1:length(all_out_TCSM), function(j) get_set_cosines_from_TCSM(TCSM_object = all_out_TCSM[[j]], true_signatures_arg = signature_definitions[,signatures_random[[j]]]))
names(cosines_from_TCSM_randomsigs) <- names(cosines_from_TCSM)
cosines_from_TCSM_df <- melt(list(true_signatures=cosines_from_TCSM,
                                  random_signatures=cosines_from_TCSM_randomsigs))
cosines_from_TCSM_df$pi_softmax=sapply(as.numeric(sapply(cosines_from_TCSM_df$L2, function(i) strsplit(i, '_')[[1]][6])), function(i) softmax(c(i, 0))[1])
cosines_from_TCSM_df$numsamples=sapply(cosines_from_TCSM_df$L2, function(i) strsplit(i, '_')[[1]][2])
cosines_from_TCSM_df$numsmuts=sapply(cosines_from_TCSM_df$L2, function(i) strsplit(i, '_')[[1]][3])
cosines_from_TCSM_df$datasetgeneration=sapply(cosines_from_TCSM_df$L2, function(i) strsplit(i, '_')[[1]][1])
cosines_from_TCSM_df$L1[cosines_from_TCSM_df$L1 == 'true_signatures'] <- 'True signatures'
cosines_from_TCSM_df$L1[cosines_from_TCSM_df$L1 == 'random_signatures'] <- 'Random signatures'
cosines_from_TCSM_df$L1 <- factor(cosines_from_TCSM_df$L1, levels=c('True signatures',  'Random signatures'))
cosines_from_TCSM_df$x = interaction( paste0('n=', cosines_from_TCSM_df$numsamples), paste0('T=', cosines_from_TCSM_df$numsmuts), sep = ', ')
cosines_from_TCSM_df$x <- factor(cosines_from_TCSM_df$x, levels=c('n=20, T=50', 'n=100, T=50', 'n=100, T=1773', 'n=100, T=Obs', 'n=150, T=Obs'))
ggplot(cosines_from_TCSM_df, aes(x= x,
                                 group=interaction(numsamples, numsmuts, pi_softmax), y=value,
                                 col=factor(pi_softmax)))+
  geom_boxplot()+
  # facet_wrap(.~datasetgeneration+L1, nrow=1, scales = 'free_x')+
  facet_grid(.~datasetgeneration+L1, scales='free', space='free_x')+
  scale_colour_manual(values = colorRampPalette(c("cyan", "yellow"))(length(unique(cosines_from_TCSM_df$pi_softmax))))+
  guides(col='none')+labs(x='Parameter combination', y='Cosine similarity')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../../../results/figures_paper/comparison_methods/recovery_signatures_TCSM.pdf", height = 3.4, width = 10)

# scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0'))
## [possibly plot]

## here we can also add the additional TwoCT with large # mutations for which we do not have results from HiLDA

##------------------------------------------------------------------------##


##------------------------------------------------------------------------##
## Match runs
names_runs <- table(c(names(all_out_TCSM), names(all_out_TMB), names(all_out_HiLDA)))
names_runs <- names(names_runs[names_runs==3])
all_out_TCSM_matched = all_out_TCSM[names_runs]
all_out_TMB_matched = all_out_TMB[names_runs]
all_out_HiLDA_matched = all_out_HiLDA[names_runs]
sapply(list(all_out_TCSM_matched, all_out_TMB_matched, all_out_HiLDA_matched), length)
stopifnot(names(all_out_TCSM_matched) == names(all_out_TMB_matched))
stopifnot(names(all_out_TCSM_matched) == names(all_out_HiLDA_matched))
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
## Comparison of beta coefficients (sorting them, as we don't have their correspondence)

# output_to_beta_matrix(all_out_TCSM_matched[[idx]], 'TCSM') ## d
# output_to_beta_matrix(all_out_TMB_matched[[idx]], 'TMB') ## d-1
# output_to_beta_matrix(all_out_HiLDA_matched[[idx]], 'HiLDA') ## d 

out_betas <- list( TCSM = lapply(all_out_TCSM_matched, output_to_beta_matrix, model='TCSM'), ## d
                   TMB = lapply(all_out_TMB_matched, output_to_beta_matrix, model='TMB'), ## d-1
                   HiLDA = lapply(all_out_HiLDA_matched, output_to_beta_matrix, model='HiLDA')) ## d
out_betas_without_HiLDA <- list( TCSM = lapply(all_out_TCSM[!(sapply(all_out_TMB, length) == 1)], output_to_beta_matrix, model='TCSM'), ## d
                   TMB = lapply(all_out_TMB[!(sapply(all_out_TMB, length) == 1)], output_to_beta_matrix, model='TMB')) ## d-1
##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
##' as we do not know which de-novo signatures correspond to the original signatures, we sort each of the 
##' results based on increasing values of the intercept/group 1 abundance (as this is likely to be more
##' stable than the slope). We do this for the data from all models. 
out_betas$HiLDA_sorted <- lapply(out_betas$HiLDA, function(i) try(i[order(i[,1]),] ))
out_betas$TCSM_sorted <- lapply(out_betas$TCSM, function(i) i[order(i[,1]),] )

## for TMB runs, first convert to softmax, so that coefficients can be compared (i.e. d coefficients in all cases)
out_betas$TMB_softmax <- lapply(out_betas$TMB, function(i) apply(rbind(i, c(0,0)), 2, softmax) )
out_betas$TMB_softmax_sorted <- lapply(out_betas$TMB_softmax, function(i) i[order(i[,1]),] )

## good agreement in intercepts, TCSM seems to recover signatures well [plot]
plot(unlist(lapply(out_betas$TCSM_sorted, function(i) i[,1])),
     unlist(lapply(out_betas$TMB_softmax_sorted, function(i) i[,1])))

plot(unlist(lapply(out_betas$TCSM_sorted, function(i) i[,2])),
     unlist(lapply(out_betas$TMB_softmax_sorted, function(i) i[,2])))

## For HiLDA, the alphas without normalising can be confounded by the dispersion parameter. 
## Hence, computing a normalised alpha, and sorting the results by this
## Moreover, the second column is the abundance compared to the first group. To get the actual
## alpha parameter we need to add them up
out_betas$HiLDA0 <- out_betas$HiLDA
out_betas$HiLDA0_sorted <- lapply(out_betas$HiLDA0, function(i) i[order(i[,1]),] )
out_betas$HiLDA <- lapply(out_betas$HiLDA0, function(i) cbind(i[,1], i[,2]+i[,1]))
out_betas$HiLDA0_sorted_SE <- lapply(out_betas$HiLDA0_sorted, function(i)  cbind(NA, ((i[,1]/sum(i)-i[,2]/sum(i))**2)))
out_betas$HiLDA_normalised <- lapply(out_betas$HiLDA, function(i) t(normalise_cl(i)))
# stopifnot(max(unlist(out_betas$HiLDA_normalised), na.rm = T) <= 1)
## there might be negative values in the second column, but not the first
out_betas$HiLDA_normalised_sorted <-  lapply(out_betas$HiLDA_normalised, function(i) try(i[order(i[,1]),] ))

list_sorted_intercepts <- list(TCSM=unlist(lapply(out_betas$TCSM_sorted, function(i) i[,1])),
                               diagREDM=unlist(lapply(out_betas$TMB_softmax_sorted, function(i) i[,1])),
                               HiLDA_normalised=unlist(lapply(out_betas$HiLDA_normalised_sorted, function(i) i[,1])),
                               HiLDA=unlist(lapply(out_betas$HiLDA, function(i) i[,1])))
pairs(list_sorted_intercepts)
## we remove the non-normalised version of HiLDA, which is indeed confounded
list_sorted_intercepts$HiLDA <- list_sorted_intercepts$HiLDA_normalised
list_sorted_intercepts$HiLDA_normalised <- NULL

list_sorted_intercepts_df <- do.call('cbind.data.frame', list_sorted_intercepts)
list_sorted_intercepts_df <- add_metadata(list_sorted_intercepts_df, rownames(list_sorted_intercepts_df))
list_sorted_intercepts_df$params <- NULL
GGally::ggpairs(list_sorted_intercepts_df, columns=1:3, 
                aes(col=factor(model), shape=n, alpha=0.2))+
  labs(x='Signature abundance in clonal group', y='Signature abundance in subclonal group')+
  scale_color_manual(values = c('#377eb8', '#ff7f00', 'red'))
ggsave("../../../../results/figures_paper/comparison_methods/intercept_ggpairs.pdf", height = 4, width = 4)
## [plot] Note that these data are compositional so a good correlation is to some extent imposed!

plt_intercept1 <- GGally::ggpairs(list_sorted_intercepts_df %>% filter(model=='C1'), columns=1:3, 
                aes(col=factor(signif(pi_softmax, 2)),  shape=n, alpha=0.2))+
  labs(x='Signature abundance in clonal group', y='Signature abundance in subclonal group')+
  # scale_colour_brewer(palette = 'BuGn')+scale_fill_brewer(palette = 'BuGn')
  scale_colour_manual(values = colorRampPalette(c("cyan", "yellow"))(length(unique(list_sorted_intercepts_df$pi))))

plt_intercept2 <- GGally::ggpairs(list_sorted_intercepts_df %>% filter(model=='C2'), columns=1:3, 
                aes(col=factor(signif(pi_softmax, 2)),  shape=n, alpha=0.2))+
  labs(x='Signature abundance in clonal group', y='Signature abundance in subclonal group')+
  # scale_colour_brewer(palette = 'BuGn')+scale_fill_brewer(palette = 'BuGn')
  scale_colour_manual(values = colorRampPalette(c("cyan", "yellow"))(length(unique(list_sorted_intercepts_df$pi))))

plt_intercept3 <- GGally::ggpairs(list_sorted_intercepts_df %>% filter(model=='C3'), columns=1:3, 
                aes(col=factor(signif(pi_softmax, 2)),  shape=n, alpha=0.2))+
  labs(x='Signature abundance in clonal group', y='Signature abundance in subclonal group')+
  scale_colour_manual(values = colorRampPalette(c("cyan", "yellow"))(length(unique(list_sorted_intercepts_df$pi))))
  #scale_colour_brewer(palette = 'BuGn')#+scale_fill_brewer(palette = 'BuGn')
plt_intercept_legend <- GGally::ggpairs(list_sorted_intercepts_df %>% filter(model=='C3'), columns=1:3, 
                                  aes(col=factor(pi),  shape=n, alpha=0.2))+
  labs(x='Signature abundance in clonal group', y='Signature abundance in subclonal group')+
  scale_colour_manual(values = colorRampPalette(c("cyan", "yellow"))(length(unique(list_sorted_intercepts_df$pi))))
plt_intercept_legend_v2 <- GGally::ggpairs(list_sorted_intercepts_df %>% filter(model=='C3'), columns=1:3, 
                                        aes(col=factor(signif(pi_softmax, 2)),  shape=n, alpha=0.2))+
  labs(x='Signature abundance in clonal group', y='Signature abundance in subclonal group')+
  scale_colour_manual(values = colorRampPalette(c("cyan", "yellow"))(length(unique(list_sorted_intercepts_df$pi))))
  
plot_new_ggpairs <- function(ggpair_object, title){
  do.call('grid.arrange', c(list(ggpair_object$plots[[4]]+geom_abline(slope = 1, intercept = 0, lty='dashed')+  guides(color='none', alpha='none', shape='none'),
                               ggpair_object$plots[[7]]+geom_abline(slope = 1, intercept = 0, lty='dashed')+  guides(color='none', alpha='none', shape='none'),
                               ggpair_object$plots[[8]]+geom_abline(slope = 1, intercept = 0, lty='dashed')+  guides(color='none', alpha='none', shape='none')), nrow=1, top=title))
}

legend_pi <- cowplot::get_legend(plt_intercept_legend_v2$plots[[4]]+theme(legend.position = 'bottom')+
                      guides(alpha='none', shape='none', color=guide_legend(nrow=2,byrow=TRUE))+
                      labs(col=latex2exp::TeX('$\\pi$')))
legend_pi_v2 <- cowplot::get_legend(plt_intercept_legend$plots[[4]]+theme(legend.position = 'bottom')+
                                   guides(alpha='none', shape='none', color=guide_legend(nrow=1,byrow=TRUE))+
                                   labs(col=latex2exp::TeX('$\\pi$')))
legend_pi_v3 <- cowplot::get_legend(plt_intercept_legend_v2$plots[[4]]+theme(legend.position = 'bottom')+
                                      guides(alpha='none', shape='none', color=guide_legend(nrow=1,byrow=TRUE))+
                                      labs(col=latex2exp::TeX('$\\pi$')))
beta0_correlations_softmax <- list(plot_new_ggpairs(plt_intercept1, title = 'C1') ,
     plot_new_ggpairs(plt_intercept2, title = 'C2'),
     plot_new_ggpairs(plt_intercept3, title = 'C3'), 
     legend_pi)

pdf("../../../../results/figures_paper/comparison_methods/beta0_correlations_softmax.pdf", height = 8, width = 8, onefile = T)
cowplot::plot_grid(plotlist = beta0_correlations_softmax,
                          nrow=4, rel_heights=c(1,1,1,0.4))
dev.off()

pdf("../../../../results/figures_paper/comparison_methods/legend_pi.pdf", height = 1, width = 8, onefile = T)
cowplot::plot_grid(legend_pi)
dev.off()

pdf("../../../../results/figures_paper/comparison_methods/legend_pi_v2.pdf", height = 1, width = 10, onefile = T)
cowplot::plot_grid(legend_pi_v2)
dev.off()

pdf("../../../../results/figures_paper/comparison_methods/legend_pi_v3.pdf", height = 1, width = 12, onefile = T)
cowplot::plot_grid(legend_pi_v3)
dev.off()

##------------------------------------------------------------------------##

##------------------------------------------------------------------------##

## Beta_1, or slope
## In TMB, if the beta1 are low in absolute value, there is no DA
## In TCSM it is the same
## In HiLDA (the unchanged coefficients, stored in HiLDA0_sorted), that is too the case, 
## although it is confounded by the dispersion
list_sorted_slopes <- list(TCSM=unlist(lapply(out_betas$TCSM_sorted, function(i) i[,2])),
                           diagREDM=unlist(lapply(out_betas$TMB_softmax_sorted, function(i) i[,2])),
                           HiLDA=unlist(lapply(out_betas$HiLDA0_sorted, function(i) i[,2])))
list_sorted_slopes_df <- do.call('cbind.data.frame', list_sorted_slopes)
list_sorted_slopes_df <- add_metadata(list_sorted_slopes_df, rownames(list_sorted_slopes_df))
pairs(list_sorted_slopes) ## no correlation - possibly confouded by the absolute values of the coefficients

ggplot(list_sorted_slopes_df, aes(x=pi, y=diagREDM))+geom_point()+
  geom_violin(aes(group=pi))

## "Importance" of beta_2 compared to beta_1: we divide beta_2 by beta_1.
## Relative importance of the group compared to the abundance captured by the intercept
## Note that for HiLDA we must use out_betas$HiLDA and not out_betas$HiLDA0, as we want the alpha value
list_sorted_slopes_importance <- list(TCSM=unlist(lapply(out_betas$TCSM, function(j) sum(abs(j[,2]))/sum(abs(j[,1])))),
                                      diagREDM=unlist(lapply(out_betas$TMB, function(j) sum(abs(j[,2]))/sum(abs(j[,1])))),
                           HiLDA=unlist(lapply(out_betas$HiLDA, function(j) sum(abs(j[,2]))/sum(abs(j[,1])))))
list_sorted_slopes_importance2 <- list(TCSM=unlist(lapply(out_betas$TCSM, function(j) sum(abs((j[,2]+j[,1])/j[,1])))),
                                       diagREDM=unlist(lapply(out_betas$TMB, function(j) sum(abs((j[,2]+j[,1])/j[,1])))),
                                      HiLDA=unlist(lapply(out_betas$HiLDA, function(j) sum(abs((j[,2]+j[,1])/j[,1])))))
plot(list_sorted_slopes_importance$diagREDM,
     list_sorted_slopes_importance$TCSM)
plot(list_sorted_slopes_importance$diagREDM,
     list_sorted_slopes_importance$HiLDA)

pairs(list_sorted_slopes_importance) ## better correlation
pairs(list_sorted_slopes_importance2) ## bad correlation
## [plot]

plot(log(list_sorted_slopes_importance2$diagREDM),
     log(list_sorted_slopes_importance2$TCSM))

pairs(list_sorted_slopes_importance2, xlim=c(0,4), ylim=c(0,4))



##' Do this for the non-softmaxed and
##' the softmaxed results. In HiLDA and TMB, convert coefficients to ALR, too
# out_betas$HiLDA_ALR <- lapply(out_betas$HiLDA, function(i) try(log(sweep(i, 2, i[nrow(i),], '/') )[-nrow(i),]))
# out_betas$HiLDA_ALR[sapply(out_betas$HiLDA_ALR, typeof) == 'character'] <- NA ## failed runs
## ALR not possible for TCSM, as group coefficients can be positive or negative


## normalise alphas in HiLDA
# out_betas$HiLDA_sorted <- lapply(out_betas$HiLDA, function(i) try(i[order(i[,1]),] ))
# out_betas$HiLDA_sorted[sapply(out_betas$HiLDA_sorted, typeof) == 'character'] <- NA ## failed runs

## also good agreement in HiLDA
# plot(unlist(lapply(out_betas$TCSM_sorted[-which(is.na(out_betas$HiLDA_normalised_sorted))], function(i) i[,1])),
#      unlist(lapply(out_betas$HiLDA_normalised_sorted[-which(is.na(out_betas$HiLDA_normalised_sorted))], function(i) i[,1])))


out_betas$TCSM_sorted[[1]][,2]
plot(unlist(lapply(out_betas$TMB_softmax_sorted, function(j) as.vector(scale(j[,2], center = T)))),
     unlist(lapply(out_betas$TCSM_sorted, function(j) as.vector(scale(j[,2], center = T)))))


unlist(lapply(out_betas$TMB_sorted, function(j) as.vector(scale(j[,2], center = T))))
GGally::ggpairs(do.call('cbind.data.frame', list_sorted_slopes))+
  labs(x='Signature abundance in clonal group', y='Signature abundance in subclonal group')

out_betas$TCSM[[1]]
out_betas$TMB[[1]]
out_betas$HiLDA[[1]]


## unify
for(model in names(out_betas)){
  for(i in 1:length(out_betas[[model]])){
    if(!all(is.na(out_betas[[model]][[i]]))){
      colnames(out_betas[[model]][[i]]) <- c('1', '2')
      out_betas[[model]][[i]] <- matrix(unlist(out_betas[[model]][[i]]), nrow=nrow(out_betas[[model]][[i]]))
    }
  }
} 

for(model in names(out_betas_without_HiLDA)){
  for(i in 1:length(out_betas_without_HiLDA[[model]])){
    if(!all(is.na(out_betas_without_HiLDA[[model]][[i]]))){
      colnames(out_betas_without_HiLDA[[model]][[i]]) <- c('1', '2')
      out_betas_without_HiLDA[[model]][[i]] <- matrix(unlist(out_betas_without_HiLDA[[model]][[i]]), 
                                                      nrow=nrow(out_betas_without_HiLDA[[model]][[i]]))
    }
  }
} 

out_betas$TCSM[[1]]
out_betas$TMB[[1]]
out_betas$HiLDA[[1]]

# table(out_betas_df$Var1)
# table(out_betas_df$Var2)
# 
# head(out_betas_df)
# table(out_betas_df$Var2)

# head(out_betas_df)
# 
# table(out_betas_df$dataset)

##------------------------------------------------------------------------##

##------------------------------------------------------------------------##
give_out_betas_df <- function(list_betas){
  out_betas_df <- melt(list_betas)
  out_betas_df <- out_betas_df[!is.na(out_betas_df$value),]
  out_betas_df$Var2[out_betas_df$Var2 == '1'] <- 'Intercept'
  out_betas_df$Var2[out_betas_df$Var2 == '2'] <- 'Slope'
  out_betas_df <- dcast(out_betas_df, L2+Var1+Var2~interaction(L1), value='value')
  # out_betas_df
  out_betas_df$dataset <- gsub("_.*", "", out_betas_df$L2)
  out_betas_df2 <- melt(out_betas_df, id.vars=c('L2', 'Var1', 'Var2', 'dataset'))
  out_betas_df2 <- add_metadata(df = out_betas_df2, out_betas_df2$L2)
  out_betas_df2
}
out_betas_df2 <- give_out_betas_df(out_betas)
out_betas_without_HiLDA_df2 <- give_out_betas_df(out_betas_without_HiLDA)
out_betas_df2$variable <- as.character(out_betas_df2$variable)
out_betas_df2$variable <- gsub("TMB", "diagREDM", out_betas_df2$variable)

# mse <- sapply(out_betas$HiLDA, function(i) mean((i[,1]-i[,2])**2) )
# mse <- sapply(out_betas$HiLDA, function(i) mean((i[,1]/sum(i)-i[,2]/sum(i))**2) )
# mse <- data.frame(mse)
# mse <- add_metadata(mse, rownames(mse))
# mse$L2 <- rownames(mse)
# mse$Var1 <- 'NA' ## idx for category: here only one value per run. as we group by that variable, add this
# mse$Var2 <- 'Slope' ## Intercept or slope: here it does not apply, but we are subsetting by slope below
# mse$dataset <- mse$model
# mse$value <- mse$mse
# mse$mse <- NULL
# mse$variable <- 'HiLDA_MSE'
# sort(colnames(out_betas_df2))
# sort(colnames(mse))
# out_betas_df2 <- rbind.data.frame(out_betas_df2, mse[,colnames(out_betas_df2)])


# head(out_betas_df2)
# out_betas_df2$Var1 %>% unique
# ggplot(out_betas_df2, aes(x=factor(pi), y=value, col=Var1, group=interaction(pi, Var1)))+geom_boxplot()+
#   facet_wrap(.~interaction(Var1, dataset, model, n, nlambda), ncol=8, scales = 'free_y')
# 
# ggplot(out_betas_df2 %>% filter(model == 'HiLDA'), aes(x=factor(pi), y=value, col=Var1, group=interaction(pi, Var1)))+geom_boxplot()+
#   facet_wrap(.~interaction(Var1, dataset, model, n, nlambda), ncol=8, scales = 'free_y')

ggplot(out_betas_df2 %>% dplyr::filter(dataset == 'C1'), aes(x=factor(pi), y=value, col=Var1, group=interaction(pi)))+geom_boxplot()+
  facet_wrap(.~interaction(variable, model, n, nlambda, drop = T), ncol=8, scales = 'free_y')

ggplot(out_betas_df2 %>% dplyr::filter(dataset == 'C1', Var2 == 'Slope'), 
       aes(x=interaction(pi), y=value, col=Var1, group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(variable, dataset, n, nlambda, drop = T), ncol=8, scales = 'free_y')

## the clearest (as it should be, given it is the mixture of two cancer types)

ggplot(out_betas_df2 %>% 
         dplyr::filter(dataset == 'C1', Var2 == 'Slope',
                       variable %in% c('diagREDM_softmax_sorted', 'TCSM_sorted', 'HiLDA0_sorted_SE')),  ## HiLDA_MSE
       aes(x=interaction(round(pi_softmax, 3)), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(gsub("HiLDA0 sorted SE", "HiLDA SE sorted", gsub("_", " ", variable)),
                           ' ', dataset, ' n=', n, ', T=', nlambda, drop = T, sep = ''), nrow=1, scales = 'free_y')+
  scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0', 'black'))+
  labs(x=latex2exp::TeX('$\\pi$'), y='(Transformed) Coefficient', col='Signature')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'bottom')+
  guides(color=guide_legend(nrow=1))
ggsave("../../../../results/figures_paper/comparison_methods/beta1_boxplotc1.pdf", height = 4, width = 10)

ggplot(out_betas_df2 %>% 
         dplyr::filter(dataset == 'C2', Var2 == 'Slope',
                       variable %in% c('diagREDM_softmax_sorted', 'TCSM_sorted', 'HiLDA0_sorted_SE')),  ## HiLDA_MSE
       aes(x=interaction(round(pi_softmax, 3)), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(gsub("HiLDA0 sorted SE", "HiLDA SE sorted", gsub("_", " ", variable)), dataset,
                           ' ', dataset, ' n=', n, ', T=', nlambda, drop = T, sep = ''), nrow=2, scales = 'free_y')+
  scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0', 'black'))+
  labs(x=latex2exp::TeX('$\\pi$'), y='(Transformed) Coefficient', col='Signature')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'bottom')+
  guides(color=guide_legend(nrow=1))
ggsave("../../../../results/figures_paper/comparison_methods/beta1_boxplotc2.pdf", height = 4, width = 10)

ggplot(out_betas_df2 %>% 
         dplyr::filter(dataset == 'C3', Var2 == 'Slope',
                       variable %in% c('diagREDM_softmax_sorted', 'TCSM_sorted', 'HiLDA0_sorted_SE')),  ## HiLDA_MSE
       aes(x=interaction(round(pi_softmax, 3)), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(gsub("HiLDA0 sorted SE", "HiLDA SE sorted", gsub("_", " ", variable)),
                           ' ', dataset, ' n=', n, ', T=', nlambda, drop = T, sep = ''), nrow=2, scales = 'free_y')+
  scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0', 'black'))+
  labs(x=latex2exp::TeX('$\\pi$'), y='(Transformed) Coefficient', col='Signature')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'bottom')+
  guides(color=guide_legend(nrow=1))
ggsave("../../../../results/figures_paper/comparison_methods/beta1_boxplotc3.pdf", height = 4, width = 10)

ggplot(out_betas_without_HiLDA_df2 %>% 
         dplyr::filter(dataset == 'C1', Var2 == 'Slope',
                       variable %in% c('TMB_softmax_sorted', 'TCSM_sorted')), 
       aes(x=interaction(round(pi_softmax, 3)), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(variable, dataset, n, nlambda, drop = T), nrow=1, scales = 'free_y')+
  scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0'))+
  labs(x=latex2exp::TeX('$\\pi$'), y='(Transformed) Coefficient', col='Signature')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


unique(out_betas_df2$dataset)
## with patient pairing
# ggplot(out_betas_df2 %>% 
#          dplyr::filter(dataset == 'C3', Var2 == 'Slope',
#                        variable %in% c('TMB_softmax_sorted', 'TCSM_sorted', 'HiLDA_normalised_sorted')), 
#        aes(x=interaction(round(pi_softmax, 3)), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
#   facet_wrap(.~interaction(variable, dataset, n, nlambda, drop = T), nrow=2, scales = 'free_y')+
#   scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0'))+
#   labs(x=latex2exp::TeX('$\\pi$'), y='(Transformed) Coefficient', col='Signature')+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../../../../results/figures_paper/comparison_methods/beta1_boxplot2.pdf", height = 5, width = 12)

## without patient pairing
ggplot(out_betas_df2 %>% 
         dplyr::filter(dataset == 'C2', Var2 == 'Slope',
                       variable %in% c('TMB_softmax_sorted', 'TCSM_sorted', 'HiLDA_normalised_sorted')), 
       aes(x=interaction(pi), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(variable, dataset, n, nlambda, drop = T), nrow=2, scales = 'free_y')+
  scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0'))

out_betas$HiLDA_normalised[[1]]

## No apparent different in TMB when we simulate using patient-specific intercepts or not
ggplot(out_betas_df2 %>% 
         dplyr::filter(dataset %in% c('C2', 'C3'), Var2 == 'Slope',
                       variable %in% c('TMB_softmax_sorted')), 
       aes(x=interaction(pi), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(variable, dataset, n, nlambda, drop = T), nrow=2, scales = 'free_y')+
  scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0'))

## Analysing HiLDA results
ggplot(out_betas_df2 %>% 
         dplyr::filter(variable %in% c('HiLDA')), 
       aes(x=interaction(pi), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(variable, dataset, n, nlambda, Var2, drop = T), nrow=2, scales = 'free_y')+
  scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0'))
ggplot(out_betas_df2 %>% 
         dplyr::filter(dataset == 'C1', variable %in% c('HiLDA')), 
       aes(x=interaction(pi), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(variable, dataset, n, nlambda, Var2, drop = T), nrow=2, scales = 'free_y')+
  scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0'))

## HiLDA0: the slope is the subtracted alpha between the two conditions
ggplot(out_betas_df2 %>% 
         dplyr::filter(variable %in% c('HiLDA0')), 
       aes(x=interaction(pi), y=value, col=as.character(Var1), group=interaction(pi, Var1)))+geom_boxplot()+
  facet_wrap(.~interaction(variable, dataset, n, nlambda, Var2, drop = T), nrow=2, scales = 'free_y')+
  scale_color_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#572AC7', '#51D0E0'))

example_hilda <- cbind.data.frame(DA=out_betas$HiLDA$`C1_100_50_NA_NA_1_0`,
      midDA=out_betas$HiLDA$`C1_100_50_NA_NA_-1_0`,
      noDA=out_betas$HiLDA$`C1_100_50_NA_NA_-8_0`)
sweep(example_hilda, 2, colSums(example_hilda), '/')

ggplot(mse %>% filter(model == 'C1', n==100, nlambda==50), aes(x=factor(pi_softmax), y=mse, group=pi_softmax))+
  facet_wrap(.~model+n+nlambda)+geom_boxplot()
ggplot(mse, aes(x=factor(pi_softmax), y=mse, group=pi_softmax))+
  facet_wrap(.~model+n+nlambda, scales = "free")+geom_boxplot()

rhats <- get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/HiLDA/",
                                           remove_HiLDAGlobal = T, HiLDA_return='rhats')

# rhats <- lapply(all_out_HiLDA, function(i) try(python_like_select_rownames(i$BUGSoutput$summary, 'alpha|beta')[,'Rhat']))
plot(density(unlist(rhats)))
min(unlist(rhats))

##------------------------------------------------------------------------##


##------------------------------------------------------------------------##
## Comparison of runtime

runtime_list <- list(diagREDM = get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/", runtime = T),
                     HiLDA = get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/HiLDA/", runtime = T, remove_HiLDAGlobal = T), ### HiLDA: only keep one of the two sets of results wherever both local and global tests are considered
                     TCSM = get_inference_files(folder_in = "../../../../data/assessing_models_simulation/inference_results/TCSM/", runtime = T))
sapply(runtime_list, function(i) head(names(i)))

## TCSM: two files of time, as it is a multi-step process. Add them up
runtime_list_TCSM_check = (table(duplicated(gsub(".time.*", "", names(runtime_list$TCSM)))))
stopifnot((length(runtime_list_TCSM_check) == 2) & (runtime_list_TCSM_check[1] == runtime_list_TCSM_check[2]))

for(i in grep(".time$",  names(runtime_list$TCSM))){
  runtime_list$TCSM[[i]] = as.numeric(runtime_list$TCSM[[i]][2,1])+ as.numeric(runtime_list$TCSM[[gsub("time", "time2", names(runtime_list$TCSM)[i])]][2,1])
}
runtime_list$TCSM <- runtime_list$TCSM[grep(".time$",  names(runtime_list$TCSM))]
# runtime_list[[3]] <- unlist(runtime_list[[3]])
tosec <- function(i){
  if(units(i) == "mins"){
    as.numeric(i*60)
  }else  if(units(i) == "secs"){
    as.numeric(i)
  }else  if(units(i) == "hours"){
    as.numeric(i*60*60)
  }else{
    stop
  }
}

runtime_list$HiLDA <- lapply(runtime_list$HiLDA, tosec)
runtime_list$diagREDM <- lapply(runtime_list$diagREDM, tosec)

runtime_list_melt <- (melt(runtime_list))
# View(runtime_list_melt)
runtime_list_melt <- add_metadata(runtime_list_melt, runtime_list_melt$L2)

## dataset from generation GenerationMixtureSimulation -- add the number of mutations as real
## they belong to Generation C3 (they are equivalent to v4, as there is patient-specific information)
## here the number of mutations was not 200 (as stated) but the observed number of mutations
runtime_list_melt$nlambda[runtime_list_melt$model == ''] = "Obs"
runtime_list_melt$model[runtime_list_melt$model == ''] = "C3"
runtime_list_melt$model = sapply(runtime_list_melt$model, rename_datasets_fun)
ggplot(runtime_list_melt, aes(x=interaction(L1), group=interaction(L1, pi_softmax), y=value,
                              col = interaction( paste0('n=', n), paste0('T=', nlambda), sep = ', ')))+geom_boxplot()+
  # scale_y_log10()+
  facet_wrap(.~interaction( paste0(model, ', n=', n, ', T=', nlambda)), nrow=2, scales = 'free')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x=latex2exp::TeX('Model and $\\pi$'), y='Runtime (seconds)', col='Parameter combination')+
  scale_color_manual(values = c('#004D40', '#C86518', '#3270BB', '#A7BBAD', '#7D0B1B'))+
  theme(legend.position = 'bottom')
ggsave("../../../../results/figures_paper/comparison_methods/runtime_boxplots.pdf", height = 5, width = 8)

runtime_list_melt_xtable <- runtime_list_melt %>% group_by(L1, model, n, nlambda) %>% summarise(mean_runtime=mean(value),
                                                                    sd_runtime=sd(value),
                                                                    min_runtime=min(value),
                                                                    max_runtime=max(value))%>%
  arrange(model, as.numeric(n),
          factor(nlambda, levels=c("50","1773", "Obs" )),
          L1)
colnames(runtime_list_melt_xtable) <- c('Model', 'Dataset', 'N', 'T', 'Mean (s)', 'sd (s)', 'Min (s)', 'Max (s)')
print(xtable::xtable(runtime_list_melt_xtable), include.rownames=FALSE)
