##' Compare TMB models where signatures are extracted using the set of active signatures (ground truth)
##' from the simulation, or selecting active signatures based on other strategies

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(CompSign)
library(reshape2)
source("../../../2_inference_TMB/helper_TMB.R")

all_res_TMB_default <- grep('GenerationMixtureSimulationv6', list.dirs("../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/", full.names = T), value = T)
all_files_default <- unlist(sapply(all_res_TMB_default, list.files, full.name=T))
all_files_default <- grep('.RDS', all_files_default, value = T) ## do not read in runtime files

all_res_TMB_with_sigextract <- grep('GenerationMixtureSimulationv6', list.dirs("../../../../data/assessing_models_simulation/inference_results/TMB_with_sigextract/nlminb/", full.names = T), value = T)
all_files_sigextract <- unlist(sapply(all_res_TMB_with_sigextract, list.files, full.name=T))
all_files_sigextract <- grep('.RDS', all_files_sigextract, value = T) ## do not read in runtime files

files_list = list(diagREDMsigextraction1 = grep('diagREDMsigextraction1',
                                                all_files_sigextract, value = T),
                  diagREDMsigextraction2 = grep('diagREDMsigextraction2',
                                                all_files_sigextract, value = T),
                  diagREDMsigextraction3 = grep('diagREDMsigextraction3',
                                                all_files_sigextract, value = T),
                  diagREDM = grep('diagREDM', all_files_default, value = T))

## match files
match1 <- gsub("TMB_with_sigextract", "", gsub("diagREDMsigextraction1", "", files_list$diagREDMsigextraction1))

files_list$diagREDMsigextraction1 ## unchanged
files_list$diagREDMsigextraction2 <- files_list$diagREDMsigextraction2[match(match1, gsub("TMB_with_sigextract", "", gsub("diagREDMsigextraction2", "", files_list$diagREDMsigextraction2)))]
files_list$diagREDMsigextraction3 <- files_list$diagREDMsigextraction3[match(match1, gsub("TMB_with_sigextract", "", gsub("diagREDMsigextraction3", "", files_list$diagREDMsigextraction3)))]
files_list$diagREDM <- as.vector(files_list$diagREDM[match(match1, gsub("TMB", "", gsub("diagREDM", "", files_list$diagREDM)))])


idx_to_remove = rowSums(sapply(files_list, is.na)) > 0
for(i in 1:length(files_list)){
  files_list[[i]] <- files_list[[i]][!idx_to_remove]
}

stopifnot(all.equal(gsub("diagREDMsigextraction1", "", files_list$diagREDMsigextraction1), gsub("diagREDMsigextraction2", "", files_list$diagREDMsigextraction2)))
stopifnot(all.equal(gsub("diagREDMsigextraction1", "", files_list$diagREDMsigextraction1), gsub("diagREDMsigextraction3", "", files_list$diagREDMsigextraction3)))
stopifnot(all.equal(gsub("TMB_with_sigextract", "", gsub("diagREDMsigextraction1", "", files_list$diagREDMsigextraction1)), gsub("TMB", "", gsub("diagREDM", "", files_list$diagREDM))))

files_read <- lapply(files_list, function(i){x <- lapply(i, readRDS); names(x) <- basename(i); x})
pvals <- lapply(files_read, function(i) sapply(i, wald_TMB_wrapper))

pvals_df = do.call('cbind.data.frame', pvals)
pvals_df$pi = as.numeric(sapply(basename(files_list$diagREDM), function(i) strsplit(i, '_')[[1]][7]))
pvals_df$pi_softmax=sapply(pvals_df$pi, function(i) softmax(c(i,0))[1])
empty_to_NA <- function(i){ sapply(i, function(j) if(length(j) == 0){NA}else{j}) }

pvals_df2 <- data.frame(pvals=unlist(pvals), 
                        d = unlist(lapply(files_read, function(i) empty_to_NA(sapply(i, function(j) 1+length(python_like_select_name(j$par.fixed, 'beta'))/2)))),
                        pi = as.numeric(unlist(lapply(files_read, function(j) sapply(basename(names(j)), function(k) strsplit(k, '_')[[1]][7])))),
                        variable=rep(names(pvals), sapply(pvals, length)))
pvals_df2$pi_softmax=sapply(pvals_df2$pi, function(i) softmax(c(i,0))[1])

pvals_df = melt(pvals_df, id.vars = c('pi', 'pi_softmax'))
pvals_df <- pvals_df[!is.na(pvals_df$value),]
pvals_df$datasetnotDA = pvals_df$pi_softmax < 0.01
pvals_df_summary = pvals_df %>% group_by(pi, variable, datasetnotDA) %>% summarise(frac_DA=sum(value <= 0.05)/n())
pvals_df_summary$accuracy = sapply(1:nrow(pvals_df_summary), function(j){
  if(pvals_df_summary$datasetnotDA[j]){
    ## not DA
    1-pvals_df_summary$frac_DA[j]
  }else{
    pvals_df_summary$frac_DA[j]
  }})
ggplot(pvals_df_summary, aes(x=pi, col=variable, y=frac_DA))+geom_line()
ggplot(pvals_df_summary, aes(x=pi, col=variable, y=accuracy))+geom_line()

pvals_df$datasetnotDA

ggplot(pvals_df, aes(x=))

aucs <- lapply(names(files_list), function(model_it){
  .x <- sapply(pval_thresholds, function(pval_threshold_it){
  get_sensivity_specificity(pvals_df, model = model_it,
                            pval_threshold = pval_threshold_it)

    })
  colnames(.x) <- pval_thresholds
  .x
})
names(aucs) <- names(files_list)

aucs = dcast(melt(aucs), L1+Var2~Var1, value.var = 'value')
rename_sigextraction <- function(i){
  res <- sapply(i, function(j){
    if(j == 'diagREDM'){
      return('Simulated signatures')
    }else if(j == 'diagREDMsigextraction1'){
      return('Signatures sum > 80% exposures')
    }else if(j == 'diagREDMsigextraction2'){
      return('Signatures sum > 70% exposures')
    }else if(j == 'diagREDMsigextraction3'){
      return('Signatures active > 80% samples')
    }else{
      stop()
    }
  })
  return(res)
}
aucs$L1 <- rename_sigextraction(aucs$L1)

plt1 <- ggplot(aucs, aes(x=TPR, y=FPR, col=L1))+geom_line()+
  scale_color_manual(values = c('#377eb8', '#ff7f00', '#984ea3', 'red'))+
  labs(col='Signature extraction strategy')
plt1
## betas

all_betas = lapply(1:length(files_read), function(i){
   .x <- lapply(files_read[[i]], function(j) list(plot_betas(j, plot = F, return_df = T, only_slope = T)[,1]))
   names(.x) <- basename(files_list[[i]])
   .x
  })
names(all_betas) <- names(files_read)
all_betas <- (melt(all_betas))
head(all_betas)
all_betas$pi = as.numeric(sapply(all_betas$L2, function(i) strsplit(i, '_')[[1]][7]))

all_betas$L1 <- rename_sigextraction(all_betas$L1)
all_betas$pi_softmax <- sapply(all_betas$pi, function(i) softmax(c(i,0))[1])
plt2 <- ggplot(all_betas, aes(x=factor(round(pi_softmax, 2)), y=value, group=interaction(L1,pi), col=L1))+
  geom_boxplot()+facet_wrap(.~L1, nrow=1, scales = "free_y")+
  scale_color_manual(values = c('#377eb8', '#ff7f00', '#984ea3', 'red'))+
  guides(col='none')+
  geom_hline(yintercept = 0, lty='dashed')+
  labs(x=latex2exp::TeX('$\\pi$'), y=latex2exp::TeX('$\\beta_1$'))+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

plt1
plt2

## number of signatures per run
rownames(pvals_df2) <- NULL
head(pvals_df2)
pvals_df2_summary <- pvals_df2 %>% group_by(pi, variable, pi_softmax) %>% summarise(mean_d=mean(d),
                                                                sd_d = sd(d))
# ggplot(pvals_df2)+
  # geom_jitter(aes(x=factor(round(pi_softmax, 2)), y=d, group=interaction(pi_softmax, variable),
  #                 col=variable))+
pvals_df2_summary$variable <- rename_sigextraction(pvals_df2_summary$variable)
plt3 <- ggplot(pvals_df2_summary, aes(x=factor(round(pi_softmax, 2)), ymin=mean_d-sd_d, ymax=mean_d+sd_d,
                  group = variable, fill=variable))+
  geom_ribbon(alpha=0.2)+
  geom_line(aes(y=mean_d, col=variable))+
  labs(x=latex2exp::TeX('$\\pi$'), y=latex2exp::TeX('Number of signatures (mean$\\pm$) sd'))+
  scale_color_manual(values = c('#377eb8', '#ff7f00', '#984ea3', 'red'))+
  scale_fill_manual(values = c('#377eb8', '#ff7f00', '#984ea3', 'red'))+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
  
cowplot::plot_grid(plt3, plt1, plt2)

plt3+guides(col='none', fill='none')
ggsave("../results/figures_paper/comparison_methods/active_signature_strategies_CompSign/active_signature_strategies_CompSign_numsigs.png",
       width = 3, height=3)
ggsave("../results/figures_paper/comparison_methods/active_signature_strategies_CompSign/active_signature_strategies_CompSign_numsigs.pdf",
       width = 3, height=3)

plt2+guides(col='none', fill='none')
ggsave("../results/figures_paper/comparison_methods/active_signature_strategies_CompSign/active_signature_strategies_CompSign_betas.png",
       width = 10, height=3)
ggsave("../results/figures_paper/comparison_methods/active_signature_strategies_CompSign/active_signature_strategies_CompSign_betas.pdf",
       width = 10, height=3)

plt1
ggsave("../results/figures_paper/comparison_methods/active_signature_strategies_CompSign/active_signature_strategies_CompSign_ROC.png",
       width = 5, height=3)
ggsave("../results/figures_paper/comparison_methods/active_signature_strategies_CompSign/active_signature_strategies_CompSign_ROC.pdf",
       width = 5, height=3)
