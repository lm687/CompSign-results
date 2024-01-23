B4_samples_dirs <- list.files("../../data/assessing_models_simulation/inference_results/TMB/nlminb/", full.names = T)
B4_samples_dirs <- grep('nPCAWG6_*', B4_samples_dirs, value = T)
B4_samples <- sapply(B4_samples_dirs, list.files, full.names=T)
B4_samples

B4_samples <- B4_samples[!grepl('low4nlambdaPCAWG6|low5nlambdaPCAWG6', B4_samples)]
B4_samples <- sapply(unlist(B4_samples), readRDS)


pvals <- sapply(B4_samples, wald_TMB_wrapper)

names(pvals) <- gsub("REDM.*", "REDM", gsub("_betainterceptPCAWG6_betaslopePCAWG6_covmatFULLPCAWG6", "",
                     gsub("multiple_GenerationJnorm_nPCAWG6_", "", basename(names(pvals)))))
names(pvals)

pvals <- pvals[grep('diagREDM', names(pvals))]
table(names(pvals))

## different number of mutations
pvals1 <- pvals[grep('nlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_diagREDM|lownlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_diagREDM|low2nlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_diagREDM', names(pvals))]

## different number of samples
pvals2 <- pvals[grep('nlambdaPCAWG6_lambdaPCAWG6_dPCAWG6_0_diagREDM|nlambdaPCAWG6_lowlambdaPCAWG6_dPCAWG6_0_diagREDM|nlambdaPCAWG6_highlambdaPCAWG6_dPCAWG6_0_diagREDM', names(pvals))]

sapply(split(pvals1, f = names(pvals1)), function(i) mean(i <= 0.05, na.rm=T)) ## all DA
