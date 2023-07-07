rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
require(TMB)
require(reshape2)
require(ggplot2)
source("../../../2_inference_TMB/helper_TMB.R")

# idx_dataset_betaintercept = 1 ## idx for beta intercept true parameters
# idx_dataset_betaslope = 2 ## idx for beta slope true parameters
# idx_dataset_cov = 1 ## idx for beta slope true parameters
idx_dataset_betaintercept = "1d4" ## idx for beta intercept true parameters
idx_dataset_betaslope = "1d4" ## idx for beta slope true parameters
idx_dataset_cov = "1d4" ## idx for beta slope true parameters

plot_only_converged <- TRUE
if(plot_only_converged){
  add_convergence <- '_onlyconverged_'
}else{
  add_convergence <- 'withnonconverged_'
}

model <- 'fullREDMsinglelambda_'
name_dataset <- "multiple_GenerationHnorm_"

# model <- "fullREM_"
# name_dataset <- "multiple_generationMGnorm_"


# model <- 'diagREDMsinglelambda'
# model <- 'fullREDMsinglelambda_'
# model <- 'diagREDM'
# model <- "fullREM"
# model <- "fullREDM_"
# model <- "diagREM"
# model <- "FEDMsinglelambda"
optimiser <- 'nlminb'

# name_dataset <- "multiple_GenerationHnorm_"
# name_dataset <- "multiple_GenerationCnormdiagRE_"
# name_dataset <- "multiple_GenerationCnorm_"
# name_dataset <- "multiple_GenerationGnorm_"
# name_dataset <- "multiple_GenerationMGnorm_"
# name_dataset <- "multiple_generationMGnorm_"
# name_dataset <- "multiple_generationMGnorm_80_180_100_6_0_"
# name_dataset <- "multiple_GenerationDMFE1_80_180_100_3_0_"
name_dataset0 <- paste0(strsplit(name_dataset, '_')[[1]][1:2], sep = '_', collapse = '')


fles = list.files("../../../../data/assessing_models_simulation/datasets/", full.names = T)
fles <- fles[grep(paste0(name_dataset, "*"), fles)]
fles <- fles[grep(paste0("betaintercept", idx_dataset_betaintercept, "_betaslope", idx_dataset_betaslope, "_covmat", idx_dataset_cov, "*"), fles)]
x <- lapply(fles, readRDS)
lst <- list.files(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/", optimiser, "/"), full.names = T)
lst <- lst[grepl(name_dataset, lst)]
lst <- lst[grepl(model, lst)]
lst <- lst[grep(paste0("betaintercept", idx_dataset_betaintercept, "_betaslope", idx_dataset_betaslope, "_covmat", idx_dataset_cov, "*"), lst)]

all_pd <- lapply(lst, function(i){x <- readRDS(i); try(x$pdHess)})
all_pd[sapply(all_pd, typeof) == 'character'] = FALSE
all_pd_list <- as.vector(unlist(all_pd))

table(all_pd_list)
lst0 <- lst
if(plot_only_converged) lst = lst[all_pd_list]
runs <- lapply(lst, readRDS)

summaries = lapply(runs, function(i){
  summary <- summary.sdreport(i)
  summary <- summary[!grepl("u_large", rownames(summary)),]
  # summary[,rownames(summary) == "logs_sd_RE"] = exp(summary[,rownames(summary) == "logs_sd_RE"])
  # summary[,rownames(summary) == "log_lambda"] = exp(summary[,rownames(summary) == "log_lambda"])
  # rownames(summary) = gsub("logs_", "", rownames(summary))
  # rownames(summary) = gsub("log_", "", rownames(summary))
  summary
})

## by default, the columns are sorted. Therefore, the estimates also need to be resorted
# for(it_runs in 1:length(fles)){
#   .y = load_PCAWG(fles[[it_runs]], typedata = "signatures", simulation = T, path_to_data = NA)$Y
#   .order_y = order(colSums(.y), decreasing = F)
#   ## the last category is dropped, so we only care about the relative order of the first d-1 categories
#   .order_y = order(.order_y[-length(.order_y)])
#   rep(1:(sum(rownames(summaries[[it_runs]]) == 'beta')/2), each=2)
#   python_like_select_rownames(summaries[[it_runs]], 'beta')
#   #   summaries[[it_runs]][grepl('beta', rownames(summaries[[it_runs]]))] = 
#   
#   summaries[[it_runs]]
# }
  
summaries2 = sapply(summaries, function(j) j[,1])
# if(typeof(summaries2) == "list"){
#   summaries2 <- do.call('cbind', summaries2)
# }
summaries_melt = data.frame(melt(summaries2), stringsAsFactors = F)
summaries_melt$Var1 = as.character(summaries_melt$Var1)
summaries_melt$idx_param = rep(1:sum(summaries_melt$Var2 == 1), length(lst))
summaries_melt[summaries_melt$Var1 == "log_lambda", "value"] = exp(summaries_melt[summaries_melt$Var1 == "log_lambda", "value"])
if(model %in% c("fullREM_")){
  # summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"] = exp(summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"])
  summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"] = (exp(summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"]))
}else if(model %in% c("fullREDMsinglelambda_")){
  summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"] = exp(summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"])**2
}else{
  stop()
}

summaries_melt[summaries_melt$Var1 == "log_lambda", "Var1"] = "lambda"
summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "Var1"] = "sd_RE"


if(name_dataset0 %in% c("multiple_GenerationCnormdiagRE_", "multiple_GenerationCnorm_")){
  sds = rep(x[[1]]$sd_RE, x[[1]]$d-1)
}else if(name_dataset0 %in% c("multiple_GenerationMGnorm_", "multiple_generationMGnorm_")){
  sds = x[[1]]$sd_RE
}else if(name_dataset0 %in% c("multiple_GenerationHnorm_")){
  cov_mat_true = readRDS(paste0("../../../../data/assessing_models_simulation/additional_files/multiple_fixed_covmat", idx_dataset_cov, ".RDS"))
  sds <- diag(cov_mat_true)
  covs <- cov_mat_true[upper.tri(cov_mat_true)]
  covs_true <- TRUE
}

if(model %in% c("fullREM_")){
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                rep(0,((x[[1]]$d-1)**2-(x[[1]]$d-1))/2), ## covariances RE
                sds) ##sd RE ## this is particular to GenerationCnorm
}else if(model %in% c("diagREM_")){
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                sds) ##sd RE ## this is particular to GenerationCnorm
}else if(model %in% c("fullREDM_", "fullREDMsinglelambda_")){
  if(model == "fullREDM_"){
    overdisp <- x[[1]]$lambda
    if(name_dataset0 %in% c("multiple_GenerationMGnorm_", "multiple_generationMGnorm_")){
      overdisp <- c(NA,NA)
    }
  }else if(model=="fullREDMsinglelambda_"){
    if(!is.na(x[[1]]$lambda)){
      stopifnot(x[[1]]$lambda[1] == x[[1]]$lambda[2])
    }
    overdisp <- x[[1]]$lambda[1]
  }
  
  if(covs_true){
    ## we have covariances
    covs <- covs
  }else{
    covs <- rep(0,((x[[1]]$d-1)**2-(x[[1]]$d-1))/2)
  }
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                sds, ##sd RE ## this is particular to GenerationCnorm
                covs, ## covariances RE
                overdisp)
}else if(model %in% c("diagREDM_", "diagREDMsinglelambda_" )){
  if(model == "diagREDM_"){
    if(name_dataset0 %in% c("multiple_GenerationMGnorm_", "multiple_generationMGnorm_")){
      overdisp = c(NA, NA)
      }else{
      overdisp <- x[[1]]$lambda
    }
  }else if(model=="diagREDMsinglelambda_"){
    stopifnot(x[[1]]$lambda[1] == x[[1]]$lambda[2])
    overdisp <- x[[1]]$lambda[1]
  }
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                sds, ##sd RE ## this is particular to GenerationCnorm
                overdisp)
}else if(model %in% c("FEDMsinglelambda_" )){
  overdisp <- x[[1]]$lambda
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                overdisp)
}else{
  stop('<Specify correct model>')
}
summaries_melt$subtract = summaries_melt$value - rep(true_vals, length(lst))

ggplot(summaries_melt, aes(x=idx_param, y=subtract, group=idx_param))+
  geom_abline(slope = 0, intercept = 0, lty='dashed', col='blue')+
  geom_boxplot()+facet_wrap(~Var1, scales = "free", nrow=1)+labs(x="Parameter", y="Bias")
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/bias_and_coverage/setsim_", 
            name_dataset, optimiser, '_', model, idx_dataset_betaintercept, '_',
            idx_dataset_betaslope, '_', idx_dataset_cov, add_convergence, "bias.pdf"), width = 10, height = 3.5)
# summaries_mat = matrix(summaries_melt$value, nrow=sum(summaries_melt$Var2 == 1))


(x[[1]]$beta)[1,1]
sapply(summaries, function(j) j[1,1])

(x[[1]]$beta)[2,1]
sapply(summaries, function(j) j[2,1])

boxplot(sapply(summaries, function(j) j[2,1]))
abline(h = (x[[1]]$beta)[2,1])

boxplot(summaries_melt[c(T,rep(F,10)),'value'])
abline(h = true_vals[1])

x[[1]]$lambda
exp(sapply(summaries, function(j) j[rownames(j) == "log_lambda",1]))

## compute coverage
## for each estimate, check if the true value is in the 95% confidence interval computed from the standard errors
summaries[[1]]
true_vals
confints <- sapply(summaries, function(it_run){
sapply(1:length(true_vals), function(i){
  confintint = c(it_run[i,1]-1.96*it_run[i,2],it_run[i,1]+1.96*it_run[i,2])
  (true_vals[i] > confintint[1]) & (true_vals[i] < confintint[2])
  } )
})
rownames(confints) = rownames(summaries[[1]])

df_coverage <- cbind.data.frame(parameter=make.names(rownames(confints), unique = T), CI=apply(confints, 1, mean, na.rm=T))
df_coverage$type_param = gsub("\\..*", "", df_coverage$parameter)
df_coverage$type_param[df_coverage$type_param == "beta"] = c('beta_intercept', 'beta_slope')
ggplot(df_coverage,
       aes(x=parameter, y=CI, group=1))+geom_abline(slope = 0, intercept = 0.95, col='blue', lty='dashed')+
  geom_line()+geom_point()+facet_wrap(.~type_param, scales = "free_x", nrow=1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/bias_and_coverage/setsim_",name_dataset,
              optimiser,'_',  model, idx_dataset_betaintercept, '_', idx_dataset_betaslope, '_', idx_dataset_cov, add_convergence,
              "coverage.pdf"), width = 10, height = 3.5)

# summaries[[1]]
# true_vals

# sapply(x, `[`, 'beta')

as.vector(x[[1]]$beta)
python_like_select_rownames(summaries[[1]], 'beta')[,1]

as.vector(x[[1]]$lambda)
exp(python_like_select_rownames(summaries[[1]], 'log_lambda')[1])

## What converged, what didn't?
## (what are the difference between the runs?)

runs_nonconverged <- lapply(lst0[!all_pd_list], readRDS)
runs

runs_nonconverged_list <- lapply(runs_nonconverged, function(j) try(j$par.fixed))
runs_nonconverged_list <- runs_nonconverged_list[sapply(runs_nonconverged_list, typeof) == "double"]

comparison_conv_notconv <- rbind(cbind.data.frame(melt(runs_nonconverged_list), conv='not_conv'),
      cbind.data.frame(melt(lapply(runs, function(j) try(j$par.fixed))), conv='converged'))
comparison_conv_notconv$param_name = make.names(names(runs[[1]]$par.fixed), unique = T)

comparison_conv_notconv$param_type <- gsub("\\..*", "", comparison_conv_notconv$param_name)
p <- ggplot(comparison_conv_notconv, aes(x=interaction(conv, param_name), y=value, col=conv))+
  geom_boxplot()+facet_wrap(.~param_type, scales = "free", nrow=1)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

gp <- ggplotGrob(p)
# get gtable columns corresponding to the facets (5 & 9, in this case)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]

# get the number of unique x-axis values per facet (1 & 3, in this case)
x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                function(l) length(l$range$range))

# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
gp$widths[facet.columns] <- gp$widths[facet.columns] * c(1, 2, 0.3, 1)
# gp <- ggplotGrob(a)

pdf(paste0("../../../../results/results_TMB/simulated_datasets/bias_and_coverage/setsim_",name_dataset,
           optimiser,'_',  model, idx_dataset_betaintercept, '_', idx_dataset_betaslope, '_', idx_dataset_cov, add_convergence,
           "comparison_nonconv.pdf"), width = 10, height = 3.5)
grid::grid.draw(gp)
dev.off()
ggplot(comparison_conv_notconv, aes(x=interaction(conv, param_name), y=value, col=conv))+
  geom_violin()+geom_point(alpha=0.1)+facet_wrap(.~param_name, scales = "free", nrow=3)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/bias_and_coverage/setsim_",name_dataset,
              optimiser,'_',  model, idx_dataset_betaintercept, '_', idx_dataset_betaslope, '_', idx_dataset_cov, add_convergence,
              "comparison_nonconv2.pdf"), width = 15, height = 5.5)
