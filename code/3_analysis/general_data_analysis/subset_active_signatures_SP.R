
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../")

library(TMB)
library(optparse)
library(gridExtra)
library(ggrepel)

source("1_create_ROO/roo_functions.R")
source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("2_inference_TMB/mm_multinomial/helper_functions.R")
source("2_inference_TMB/helper_TMB.R")
source("../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
source("3_analysis/helper/pcawg.colour.palette.R")

enough_samples = read.table("../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]

pcawg_palette <- pcawg.colour.palette(gsub("\\..*", "", enough_samples),
                                      scheme = "tumour.subtype")
names(pcawg_palette) <- enough_samples


df_all = lapply(enough_samples, function(ct){
  dataset <- load_PCAWG(ct = paste0("../data/roo/", ct, "_signaturesPCAWG_ROO.RDS"),
                     typedata = "signaturesPCAWG",
                     simulation = F,
                     path_to_data = NA, read_directly=T)
  cbind.data.frame(average_norm_exp_nonzero =apply(normalise_rw(dataset$Y), 2, function(i) mean(i[i > 0])),
                   average_nonzero = colMeans(normalise_rw(dataset$Y) > 0),
                   label_signature=colnames(dataset$Y))
})
names(df_all) <- enough_samples

pdf("../results/exploratory/subset_active_signatures_sp.pdf")
for(ct in enough_samples){
  
  print(ggplot(df_all[[ct]])+
    aes(x=average_nonzero, y=average_norm_exp_nonzero, label=label_signature)+
    geom_point()+geom_label_repel()+theme_bw()+ggtitle(ct)+lims(x=c(0,1), y = c(0,1)))
}
dev.off()


df_all_melt <- melt(df_all, id.vars=colnames(df_all[[1]]))
head(df_all_melt)
ggplot(df_all_melt, aes(x=average_nonzero, y=average_norm_exp_nonzero,
                        label=label_signature, col=L1))+
  geom_point()+
  scale_color_manual(values = pcawg_palette)+theme_bw()+
  theme(legend.position = "bottom")+#geom_label()
  facet_wrap(.~label_signature)
ggsave("../results/exploratory/subset_active_signatures_sp_all.pdf", width = 8)

ggplot(df_all_melt, aes(x=average_nonzero, y=average_norm_exp_nonzero,
                        label=label_signature, col=label_signature))+
  geom_point()+
  theme(legend.position = "bottom")+
  theme_bw()
ggsave("../results/exploratory/subset_active_signatures_sp_all2.pdf", width = 8)


ggplot(df_all_melt, aes(x=average_nonzero, y=average_norm_exp_nonzero,
                        label=gsub('SBS', '', label_signature), col=L1))+
  geom_point()+
  theme(legend.position = "bottom")+
  theme_bw()+
  scale_color_manual(values = pcawg_palette)+
  geom_label_repel()
ggsave("../results/exploratory/subset_active_signatures_sp_all3.pdf", width = 8, height = 5)

ggplot(df_all_melt, aes(x=average_nonzero, y=average_norm_exp_nonzero,
                        label=gsub('SBS', '', label_signature), col=L1))+
  geom_point()+
  theme(legend.position = "bottom")+
  theme_bw()+
  scale_color_manual(values = pcawg_palette)+
  geom_label_repel()+
  scale_y_continuous(trans = "log2")
ggsave("../results/exploratory/subset_active_signatures_sp_all4.pdf", width = 8, height = 5)

