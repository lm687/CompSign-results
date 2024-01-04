rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(pheatmap)
library(ComplexHeatmap)

x <- read.table("../../../data/cosmic/active_signatures_PCAWGpaper.txt")
x2 <- read.table("../../../data/pcawg/PCAWG_sigProfiler_SBS_signatures_in_samples.csv", sep = ",", hea=T) ## for the name of the cancer types in correct upper/lower case
x2 <- x2$Cancer.Types

for(i in 1:nrow(x)){
  cat((x2[match(x[i,1], toupper(x2))]), '&', paste0(names(x[i,-1])[which(x[i,-1] == 1)], collapse = ", "), "\\\\\n")
}

rownames(x) <- x2[match(x[,1], toupper(x2))]
x <- x[,-1]

pheatmap::pheatmap(x)
ht <- ComplexHeatmap::Heatmap(as(x, 'matrix'), cluster_columns = F, cluster_rows = F, col=c('#e7feff', '#f08080'),
                        name='Signature active', heatmap_legend_param = list(
                          legend_direction = "horizontal", 
                          legend_width = unit(5, "cm"),
                          nrow=1
                        ))

pdf("../../../results/exploratory/about_cosmic_sigs/active_sigs_in_PCAWG.pdf", width = 10)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

## co-occurrence
co_oc <- outer(1:ncol(x), 1:ncol(x), Vectorize(function(i, j){
  sum(x[,i] & x[,j])/(sum(sum(x[,i] & x[,j]) + sum(x[,i] & !x[,j]) ))
}))
colnames(co_oc) <- rownames(co_oc) <- colnames(x)

x[,c('SBS8', 'SBS39')]

pdf("../../../results/exploratory/about_cosmic_sigs/active_sigs_in_PCAWG_co_occurrence.pdf", width = 10, height = 10)
pheatmap(co_oc)
dev.off()
