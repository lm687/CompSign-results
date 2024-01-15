##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Creating a table for nonexogenous and active signatures
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
library(xtable)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../../data/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
nonexogenous = read.table("../../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t", comment.char = "#", fill = F)

nonexogenous
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
fles <- grep("signaturesPCAWG_ROO", list.files("../../../data/roo/", full.names = T), value = T)
DM_active_sigs <- lapply(fles,
                          function(i){a <- readRDS(i); colnames(attr(a,"count_matrices_active")[[1]])})
names(DM_active_sigs) <- gsub("_signaturesPCAWG_ROO.RDS", "", basename(fles))
DM_active_sigs <- DM_active_sigs[enough_samples]
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
df_table <- data.frame(#names(DM_active_sigs),
                       sapply(DM_active_sigs, paste0, collapse=', '))
xtable::xtable(df_table, include.rownames=FALSE)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
## Change manually:
##' 1. marking nonexogenous signatures with asterisk
##' 2. \begin{tabular}{R{1.5in}L{4in}}
##' 3. remove & 1 & 2 \\ to Cancer type & Active signatures\\
##' 4. Add caption: \caption{Active signatures in each of the cancer types, according to the PCAWG analyses. With asterisks, signatures which are considered exogenous for these analysis: those include signatures known to represent exogenous mutational processes, as well as SBS1 and SBS5. The last signature is used as baseline. Should this signature be exogenous, the last signature without an asterisk is used as baseline. \label{active_and_nonexogenous_signatures}}




# nonexogenous$V1
# [1] "SBS1"  "SBS4"  "SBS5"  "SBS7a" "SBS7b" "SBS7c" "SBS7d" "SBS11" "SBS29" "SBS31"
# [11] "SBS32" "SBS35" "SBS87" "SBS92" "SBS27" "SBS43" "SBS45" "SBS46" "SBS47" "SBS48"
# [21] "SBS49" "SBS50" "SBS51" "SBS52" "SBS53" "SBS54" "SBS55" "SBS56" "SBS57" "SBS58"
# [31] "SBS59" "SBS60"


% latex table generated in R 4.3.1 by xtable 1.8-4 package
% Wed Dec  6 17:16:02 2023
\begin{table}[ht]
\centering
\begin{tabular}{rll}
\hline
& 1 & 2 \\ 
\hline
Bone-Osteosarc & Bone-Osteosarc & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS8, SBS13, SBS17a, SBS17b, SBS30, SBS40 \\ 
Breast-AdenoCA & Breast-AdenoCA & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS8, SBS9, SBS13, SBS17a, SBS17b, SBS18, SBS37, SBS40, SBS41 \\ 
CNS-GBM & CNS-GBM & SBS1$^*$, SBS5$^*$, SBS11$^*$, SBS30, SBS37, SBS40 \\ 
CNS-Medullo & CNS-Medullo & SBS1$^*$, SBS5$^*$, SBS8, SBS18, SBS39, SBS40 \\ 
CNS-PiloAstro & CNS-PiloAstro & SBS1$^*$, SBS5$^*$, SBS19, SBS23, SBS40 \\ 
ColoRect-AdenoCA & ColoRect-AdenoCA & SBS1$^*$, SBS5$^*$, SBS10a, SBS10b, SBS15, SBS17a, SBS17b, SBS18, SBS28, SBS37, SBS40, SBS44, SBS45$^*$ \\ 
Eso-AdenoCA & Eso-AdenoCA & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS13, SBS17a, SBS17b, SBS18, SBS28, SBS40 \\ 
Head-SCC & Head-SCC & SBS1$^*$, SBS2, SBS3, SBS4$^*$, SBS5$^*$, SBS7a$^*$, SBS7b$^*$, SBS7d$^*$, SBS13, SBS16, SBS17a, SBS17b, SBS18, SBS33, SBS40 \\ 
Kidney-ChRCC & Kidney-ChRCC & SBS1$^*$, SBS2, SBS5$^*$, SBS13, SBS17a, SBS17b, SBS29$^*$, SBS40 \\ 
Kidney-RCC.clearcell & Kidney-RCC.clearcell & SBS1$^*$, SBS2, SBS5$^*$, SBS13, SBS22, SBS29$^*$, SBS40, SBS41 \\ 
Kidney-RCC.papillary & Kidney-RCC.papillary & SBS1$^*$, SBS2, SBS5$^*$, SBS13, SBS22, SBS29$^*$, SBS40, SBS41 \\ 
Liver-HCC & Liver-HCC & SBS1$^*$, SBS4$^*$, SBS5$^*$, SBS6, SBS9, SBS12, SBS14, SBS16, SBS17a, SBS17b, SBS18, SBS19, SBS22, SBS24, SBS26, SBS28, SBS29$^*$, SBS30, SBS31$^*$, SBS35$^*$, SBS40, SBS53, SBS54$^*$, SBS56$^*$ \\ 
Lung-SCC & Lung-SCC & SBS1$^*$, SBS2, SBS4$^*$, SBS5$^*$, SBS8, SBS13 \\ 
Lymph-BNHL & Lymph-BNHL & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS6, SBS9, SBS13, SBS17a, SBS17b, SBS34, SBS36, SBS37, SBS40, SBS56$^*$ \\ 
Lymph-CLL & Lymph-CLL & SBS1$^*$, SBS5$^*$, SBS9, SBS40 \\ 
Ovary-AdenoCA & Ovary-AdenoCA & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS8, SBS13, SBS18, SBS26, SBS35$^*$, SBS39, SBS40, SBS41 \\ 
Panc-AdenoCA & Panc-AdenoCA & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS6, SBS8, SBS13, SBS17a, SBS17b, SBS18, SBS20, SBS26, SBS28, SBS30, SBS40, SBS51$^*$ \\ 
Panc-Endocrine & Panc-Endocrine & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS6, SBS8, SBS9, SBS11$^*$, SBS13, SBS26, SBS30, SBS36, SBS39 \\ 
Prost-AdenoCA & Prost-AdenoCA & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS8, SBS13, SBS18, SBS33, SBS37, SBS40, SBS41, SBS45$^*$, SBS52$^*$, SBS58$^*$ \\ 
Skin-Melanoma.cutaneous & Skin-Melanoma.cutaneous & SBS1$^*$, SBS2, SBS5$^*$, SBS7a$^*$, SBS7b$^*$, SBS7c$^*$, SBS7d$^*$, SBS13, SBS17a, SBS17b, SBS38, SBS40, SBS58$^*$ \\ 
Stomach-AdenoCA & Stomach-AdenoCA & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS9, SBS13, SBS15, SBS17a, SBS17b, SBS18, SBS20, SBS21, SBS26, SBS28, SBS40, SBS41, SBS43$^*$, SBS44, SBS51$^*$, SBS58$^*$ \\ 
Thy-AdenoCA & Thy-AdenoCA & SBS1$^*$, SBS2, SBS5$^*$, SBS13, SBS40, SBS58$^*$ \\ 
Uterus-AdenoCA & Uterus-AdenoCA & SBS1$^*$, SBS2, SBS3, SBS5$^*$, SBS6, SBS10a, SBS10b, SBS13, SBS14, SBS15, SBS26, SBS28, SBS40, SBS44 \\ 
\hline
\end{tabular}
\end{table}