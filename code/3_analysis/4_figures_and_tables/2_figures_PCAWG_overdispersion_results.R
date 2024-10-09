##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Replotting overdispersion results
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
library(TMB)
library(ggplot2)
library(latex2exp)
library(tikzDevice)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../../data/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
diagRE_DMDL_SP <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../results/results_TMB/pcawg_robjects/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS")))
}, simplify = F); names(diagRE_DMDL_SP) <- enough_samples
dataset_all_sigs  <- sapply(enough_samples, function(ct){
  load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../data/", load_all_sigs = T, override_warning_X_Z = T)
}, simplify = F)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
ovrdisp <- do.call('rbind.data.frame', lapply(1:length(diagRE_DMDL_SP), try(function(idx){
  if(diagRE_DMDL_SP[[idx]]$pdHess){
    cbind.data.frame( plot_lambdas(diagRE_DMDL_SP[[idx]], return_df=T, plot=F), ct=names(diagRE_DMDL_SP)[idx])
  }else{
    c(NA, NA)
  }
})))
ovrdisp[ovrdisp$name == 'Lambda 1','name'] = 'Clonal'
ovrdisp[ovrdisp$name == 'Lambda 2','name'] = 'Subclonal'

differential_precision <- p.adjust(sapply(diagRE_DMDL_SP, wald_TMB_wrapper_overdisp), method = 'fdr')
names(differential_precision) <- names(diagRE_DMDL_SP)
sort(differential_precision)
table(differential_precision <= 0.05)
differential_precision[(differential_precision <= 0.05)]

ovrdisp$differentially_abundant = ifelse(ovrdisp$ct %in% names(differential_precision[(differential_precision <= 0.05)]), yes = '*', no = '')
ovrdisp$differentially_abundant

differential_precision_2 <- p.adjust(sapply(diagRE_DMDL_SP, ttest_TMB_wrapper_overdisp), method = 'fdr')
names(differential_precision_2) <- names(diagRE_DMDL_SP)
sort(differential_precision_2)
table(differential_precision_2 <= 0.05)
differential_precision_2[(differential_precision_2 <= 0.05)]

ovrdisp$differential_precision_2 = ifelse(ovrdisp$ct %in% names(differential_precision_2[(differential_precision_2 <= 0.05)]), yes = '*', no = '')

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Relationship between differential dispersion and the number of mutations
ovrdisp %>% group_by(ct) %>% summarise(higher_disp_in_subclonal = ifelse())

ovrdisp$higher_disp_in_subclonal <- rep(c(ovrdisp[c(T,F),'Estimate'] > ovrdisp[c(F,T),'Estimate']), each=2) ## higher dispersion in subclonal

ovrdisp$higher_disp_in_subclonal

num_muts_clonal_subclonal = sapply(dataset_all_sigs, function(x) sapply(split_matrix_in_half(x$Y), sum))
num_muts_clonal_subclonal[1,] > num_muts_clonal_subclonal[2,]
higher_nmuts_in_clonal = num_muts_clonal_subclonal[1,]>num_muts_clonal_subclonal[2,]

higher_disp_in_subclonal_vec = ovrdisp$higher_disp_in_subclonal[c(T,F)]
names(higher_disp_in_subclonal_vec) = ovrdisp$ct[c(T,F)]

stopifnot(names(higher_disp_in_subclonal_vec) == colnames(num_muts_clonal_subclonal))
table(higher_nmuts_in_clonal, higher_disp_in_subclonal_vec)

## most ct have more mutations in the clonal group and higher dispersion in the subclonal group
names(which(higher_nmuts_in_clonal & higher_disp_in_subclonal_vec))

## 4 ct have higher mutations in the clonal group and higher dispersion in the clonal group
names(which(higher_nmuts_in_clonal & !higher_disp_in_subclonal_vec)) #"Head-SCC"                "Kidney-RCC.clearcell"    "Kidney-RCC.papillary"    "Skin-Melanoma.cutaneous"

## 3 ct have lower mutations in the clonal group and higher dispersion in the clonal group
names(which(!higher_nmuts_in_clonal & higher_disp_in_subclonal_vec)) #"CNS-GBM"          "CNS-PiloAstro"    "ColoRect-AdenoCA"

## none with lower mutations in the clonal group and lower dispersion in the clonal group
names(which(!higher_nmuts_in_clonal & !higher_disp_in_subclonal_vec))

##-----------------------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------------------##

overdispersion_plot <- ggplot(ovrdisp, aes(x=ct,  y=`Estimate`, group=name, col=name, shape=name))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=`Estimate`-`Std..Error`,
                    ymax=`Estimate`+`Std..Error`), width=.1, position=position_dodge(width=0.5))+
  theme_bw()+
  geom_text(aes(y=Inf, label=differential_precision_2, vjust=1.8), col='black')+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='', 
       y='CHANGE THIS', ## change to  $\ln(\hat{\lambda})$
       # y= TeX(r"($\ln(\hat{\lambda})$)"), 
       col='Group', shape='Group')+
  theme(legend.position = "bottom")+
  theme(        legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(-10,-10,-10,-10),
                plot.margin = unit(c(1,1,1,1), "cm"))+
  scale_color_manual(values=c('#3b4d61', '#ef9d10'))

tikzDevice::tikz(file = "../../../results/figures_paper/overdispersion_params_groups_test_2_updated.tex", height = 3.5, width = 6.5)
overdispersion_plot
dev.off()


##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##

# ovrdisp <- do.call('rbind.data.frame', lapply(1:length(diagRE_DMDL), try(function(idx){
#   if(diagRE_DMDL[[idx]]$pdHess){
#     cbind.data.frame( plot_lambdas(diagRE_DMDL[[idx]], return_df=T, plot=F), ct=names(diagRE_DMDL)[idx])
#   }else{
#     c(NA, NA)
#   }
# })))
# ovrdisp[ovrdisp$name == 'Lambda 1','name'] = 'Clonal'
# ovrdisp[ovrdisp$name == 'Lambda 2','name'] = 'Subclonal'
# # ovrdisp$differentially_abundant = ifelse(ovrdisp$ct %in% names(differential_precision[(differential_precision <= 0.05)]), yes = '*', no = '')
# differential_precision_2 <- p.adjust(sapply(diagRE_DMDL, ttest_TMB_wrapper_overdisp), method = 'fdr')
# names(differential_precision_2) <- names(diagRE_DMDL)
# which(differential_precision_2 <= 0.05) ## which are differentially overdispersed
# length(which(differential_precision_2 <= 0.05))
# cat(names(which(differential_precision_2 <= 0.05)), sep=', ')
# head(dcast(ovrdisp, formula = ct~name, value.var=c('Estimate'))) ## in whih cases subclonal is more dispersed
# table(apply(dcast(ovrdisp, formula = ct~name, value.var=c('Estimate')), 1, function(i) i[2] > i[3]))
# 
# 
# ovrdisp$differential_precision_2 = ifelse(ovrdisp$ct %in% names(differential_precision_2[(differential_precision_2 <= 0.05)]), yes = '*', no = '')
# # tikzDevice::tikz(file = "../../results/figures_paper/overdispersion_params_groups_test_2.tex", height = 3.5, width = 6.5)
# ggplot(ovrdisp, aes(x=ct,  y=`Estimate`, group=name, col=name))+
#   geom_point(position=position_dodge(width=0.5))+
#   geom_errorbar(aes(ymin=`Estimate`-`Std..Error`,
#                     ymax=`Estimate`+`Std..Error`), width=.1, position=position_dodge(width=0.5))+
#   theme_bw()+
#   geom_text(aes(y=Inf, label=differential_precision_2, vjust=1.8), col='black')+
#   theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
#   labs(x='', y='Log lambda', col='Group')+theme(legend.position = "bottom")+
#   theme(        legend.margin=margin(0,0,0,0),
#                 legend.box.margin=margin(-10,-10,-10,-10),
#                 plot.margin = unit(c(1,1,1,1), "cm"))+
#   scale_color_manual(values=c('#3b4d61', '#ef9d10'))
# # dev.off()
## add manually in file: $\log(\widehat{\lambda})$
