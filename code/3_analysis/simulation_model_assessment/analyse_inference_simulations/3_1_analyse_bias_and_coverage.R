
library(optparse)

local <- F
if(local){
  rm(list = ls())
  # setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("/Users/morril01/Documents/PhD/GlobalDA/code")
  
  optimiser <- 'nlminb'
  
  # multiple_GenerationJnorm_200_180_20_3_0_fullREDM_betaintercept1d2_betaslope1d2_covmat1d2_dataset
  
  ## setsim_multiple_GenerationJnorm_nlminb_200_180_2_5_0_diagREDM_betaintercept1d4_betaslope1d4_covmat1d4_onlyconverged_betacorrelation.pdf
  # name_dataset <- 'multiple_GenerationMGnorm_'
  name_dataset <- 'multiple_GenerationJnorm_'
  
  # model <- "fullREDM_"
  # model <- "singleREDM_"
  model <- "diagREDM_"
  
  # idx_dataset_betaintercept <- "betaintercept3"
  # idx_dataset_betaslope <- "betaslope3"
  # idx_dataset_cov <- "sdRE1"
  idx_dataset_betaintercept <- "betaintercept1d4"
  idx_dataset_betaslope <- "betaslope1d4"
  idx_dataset_cov <- "covmat1d4"
  # idx_dataset_betaintercept <- "betaintercept1d2"
  # idx_dataset_betaslope <- "betaslope1d2"
  # idx_dataset_cov <- "covmat1d2"
  
  arg_n <- 200
  arg_nlambda <- 180
  # arg_lambda <- 2
  # arg_lambda <- 20
  arg_lambda <- 100
  
  arg_d <- 5
  arg_beta_gamma_shape <- 0
  
  lst <- list.files(paste0("../data/assessing_models_simulation/inference_results/TMB/", optimiser, "/"), full.names = T)
  lst <- lst[grepl(name_dataset, lst)]
  lst <- lst[!grepl('_NC.RDS', lst)]
  lst <- lst[grepl(model, lst)]
  lst <- lst[grep(paste0(idx_dataset_betaintercept, "_", idx_dataset_betaslope, "_", idx_dataset_cov, "*"), lst)]
  # lst <- lst[grepl('200_180_2_5_0_', lst)]
  # lst <- lst[grepl('200_180_20_3_0', lst)]
  # lst <- lst[grepl('200_180_100_5_0', lst)]
  lst <- lst[grepl(paste0(arg_n, '_', arg_nlambda, '_', arg_lambda, '_', arg_d, '_', arg_beta_gamma_shape), lst)]
  # lst <- lst[grep(paste0(dataset_name_split[1:7], collapse = '_'), lst)]
  length(lst)
  
  opt <- list()
  opt$model <- gsub("_$", "", model)
  opt$input <- lst
  opt$dataset_generation <- gsub("multiple_", "", gsub("_$", "", name_dataset))
  opt$multiple_runs <- T
  opt$run_nonconverged <- F
  opt$optimiser <- optimiser
  opt$n = arg_n
  opt$nlambda = arg_nlambda
  opt$lambda = arg_lambda
  opt$d = arg_d
  opt$beta_gamma_shape <- arg_beta_gamma_shape
  
  # dataset <- "multiple_GenerationHnorm_80_180_9_6_0_fullREDMsinglelambda_betaintercept2_betaslope2_covmat2"
  # dataset <- "multiple_GenerationMGnorm_80_180_100_6_0_fullREDMonefixedlambda_betaintercept3_betaslope3_sdRE1"
  # dataset <- "multiple_GenerationMGnorm_80_180_100_6_0_fullREDMonefixedlambda3_betaintercept3_betaslope3_sdRE1"
  # dataset <- "multiple_GenerationMGnorm_80_180_100_6_0_diagREDM_betaintercept3_betaslope3_sdRE1"
  # dataset <- "multiple_GenerationMGnorm_200_180_100_6_0_fullREDM_betaintercept3_betaslope3_sdRE1"
  # dataset <- "multiple_GenerationCnormdiagRE_80_180_9_6_0_betaintercept2_betaslope2"
  
}else{
  option_list = list(
    make_option(c("--input_list"), type="character", default=NA,
                help="Text with small description of the type of simulation being carried out", metavar="character"),
    # make_option(c("--output_folder_name"), type="character", default=NA,
    #             help="Name of the folder for the output files", metavar="character"),
    # make_option(c("--output_string"), type="character", default=NA,
    #             help="String to be attached to the output files", metavar="character"),
    make_option(c("--model"), type="character", default=NA,
                help="Name of model", metavar="character"),
    make_option(c("--dataset_generation"), type="character", default=NA,
                help="Generation name", metavar="character"),
    make_option(c("--multiple_runs"), type="logical", default=F,
                help="Boolean: are we analysing multiple runs?", metavar="logical"),
    make_option(c("--optimiser"), type="character", default='nlminb',
                help="Name of optimiser used by TMB", metavar="character"),
    make_option(c("--run_nonconverged"), type="logical", default=F,
                help="Boolean: are we analysing runs that did not converge?", metavar="logical"),
    make_option(c("--beta_intercept_input"), type="character", default=NA,
                help="Fixed intercept for the betas", metavar="character"),
    make_option(c("--beta_slope_input"), type="character", default=NA,
                help="Fixed slope for the betas", metavar="character"),
    make_option(c("--sdRE_input"), type="character", default=NA,              
                help="Fixed standard deviations and covariances for RE", metavar="character"),
    make_option(c("--d"), type="numeric", default=NA,
                help="Number of features", metavar="numeric"),
    make_option(c("--n"), type="numeric", default=NA,
                help="Number of samples", metavar="numeric"),
    make_option(c("--nlambda"), type="numeric", default=NA,
                help="Parameter lambda for Poisson draws of number of mutations in sample", metavar="numeric"),
    make_option(c("--beta_gamma_shape"), type="numeric", default=NA,
                help="Shape parameter for gamma distribution for beta (i.e. slope coefficient for changes in exposure between groups)", metavar="numeric"),
    make_option(c("--lambda"), type="numeric", default=0,
                help="Overdispersion parameter", metavar="numeric")
  );
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  
  opt$input <- strsplit(opt$input, ' ')[[1]]
  
  idx_dataset_betaintercept <- opt$beta_intercept_input
  idx_dataset_betaslope <- opt$beta_slope_input
  idx_dataset_cov <- opt$sdRE_input
  
}

require(TMB)
require(reshape2)
require(ggplot2)
require(GGally)
source("2_inference_TMB/helper_TMB.R")

input_run <- sub("_dataset[^_dataset]+$", "", basename(opt$input))
if(length(table(input_run)) > 1){
  warning('There are runs of mutiple parameters included in this analysis\n')
}

plot_only_converged <- !opt$run_nonconverged

if(plot_only_converged){
  add_convergence <- '_onlyconverged_'
}else{
  add_convergence <- 'withnonconverged_'
}

model <- opt$model
if(opt$multiple_runs){
  name_dataset <- paste0('multiple_', opt$dataset_generation, '_') 
}else{
  name_dataset <- paste0(opt$dataset_generation, '_') 
}

first_part_output <- paste0("../results/results_TMB/simulated_datasets/bias_and_coverage/setsim_", 
                            name_dataset, opt$optimiser, '_', opt$n, '_', opt$nlambda,  '_', opt$lambda,  '_', opt$d,
                            '_', opt$beta_gamma_shape,  '_', model,  '_', idx_dataset_betaintercept, '_',
                            idx_dataset_betaslope, '_', idx_dataset_cov, add_convergence)
# "../results/results_TMB/simulated_datasets/bias_and_coverage/setsim_multiple_{datasetgeneration}_nlminb_{n}_{nlambda}_{lmbda}_{d}_{beta_intensity}_{model}_{fixed_beta_intercept}_{fixed_beta_slope}_{sdRE_input}_onlyconverged_coverage_beta.pdf"

# name_dataset0 <- paste0(strsplit(name_dataset, '_')[[1]][1:2], sep = '_', collapse = '')

## Load runs
print(opt$input)
cat("Loading runs\n")

lst <- opt$input

all_pd <- lapply(lst, function(i){x <- readRDS(i); try(x$pdHess)})
all_pd[sapply(all_pd, typeof) == 'character'] = FALSE
all_pd_list <- as.vector(unlist(all_pd))

table(all_pd_list)

lst0 <- lst
if(plot_only_converged) lst = lst[all_pd_list]
runs <- lapply(lst, readRDS)
tryerror <- sapply(runs, typeof) == "list"
runs <- runs[tryerror]
lst <- lst[tryerror]

print(paste0("../data/assessing_models_simulation/datasets/",
             gsub(paste0(opt$model, "_"), "", basename(lst[1]))))

### load datasets
cat("Loading datasets...\n")
x <- lapply(lst, function(i) (readRDS(paste0("../data/assessing_models_simulation/datasets/",
                                             gsub(paste0(opt$model, "_"), "", basename(i))))))
cat("... datasets loaded.\n")

cat("Creating summaries of runs\n")
summaries = lapply(runs, function(i){
  summary <- summary.sdreport(i)
  summary <- summary[!grepl("u_large", rownames(summary)),]
  summary <- summary[!grepl("u_random_effects", rownames(summary)),] ## for  "singleREDM"
  summary
})

summaries2 = sapply(summaries, function(j) j[,1])

cat("Melting summaries of runs\n")

summaries_melt = data.frame(melt(summaries2), stringsAsFactors = F)
summaries_melt$Var1 = as.character(summaries_melt$Var1)
summaries_melt$idx_param = rep(1:sum(summaries_melt$Var2 == 1), length(lst))

cat("Transforming estimates of runs\n")

## log_lambda simply need to be exponentiated (e^x) to compare it to the simulated value
summaries_melt[summaries_melt$Var1 == "log_lambda", "value"] = exp(summaries_melt[summaries_melt$Var1 == "log_lambda", "value"])

## logs_sd_RE simply need to be exponentiated (e^x) to compare it to the simulated value
## if we want the variances, we need it to the power of two
if("logs_sd_RE" %in% unique(summaries_melt$Var1)){
  summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"] = (exp(summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"]))
  summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"] = exp(summaries_melt[summaries_melt$Var1 == "logs_sd_RE",
                                                                                    "value"])**2
}

d <- length(python_like_select_name(runs[[1]]$par.fixed, 'logs_sd_RE'))
if("cov_par_RE" %in% unique(summaries_melt$Var1)){
  cat('Getting estimated covariances\n')
  for(idx in unique(summaries_melt$Var2)){
    .x <- L_to_cov(summaries_melt[(summaries_melt$Var1 == "cov_par_RE") & (summaries_melt$Var2 == idx), "value"], d=d)
    summaries_melt[(summaries_melt$Var1 == "cov_par_RE") & (summaries_melt$Var2 == idx), "value"] <- .x[upper.tri(.x)]
  }
  ## checking positive semi-definiteness
  give_first_col <- function(i){
    if(is.null(dim(i))){i[1]}else{i[,1]}
  }
  mvtnorm::rmvnorm(n=10, mean = rep(0, d), sigma = L_to_cov(give_first_col(python_like_select_rownames(summaries[[1]], 'cov_par_RE')), d=d))
  
}

summaries_melt[summaries_melt$Var1 == "log_lambda", "Var1"] = "lambda"
summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "Var1"] = "var_RE"


if(opt$dataset_generation %in% c("GenerationCnormdiagRE", "GenerationCnorm")){
  ## shared standard deviations for random effect
  sds = rep(x[[1]]$sd_RE, x[[1]]$d-1)
  covs_true <- F
}else if(opt$dataset_generation %in% c("GenerationMGnorm", "generationMGnorm")){
  ## independent random effects; several random effects
  sds = x[[1]]$sd_RE
  covs_true <- F
}else if(opt$dataset_generation %in% c("GenerationHnorm", "GenerationJnorm")){
  ## full, unconstrained, random effects matrix
  cat('Reading covariance matrix\n')
  cov_mat_true = readRDS(paste0("../data/assessing_models_simulation/additional_files/multiple_fixed_", idx_dataset_cov, ".RDS"))
  sds <- diag(cov_mat_true)
  covs <- cov_mat_true[upper.tri(cov_mat_true)]
  covs_true <- TRUE
}else{
  stop('Specify correct <dataset_generation>\n')
}

if(model %in% c("fullREM")){
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                rep(0,((x[[1]]$d-1)**2-(x[[1]]$d-1))/2), ## covariances RE
                sds) ##sd RE ## this is particular to GenerationCnorm
  # }else if(model %in% c("diagREM")){
  #   true_vals = c(as.vector(x[[1]]$beta), ## betas
  #                 sds) ##sd RE ## this is particular to GenerationCnorm
}else if(model %in% c("fullREDM", "fullREDMsinglelambda")){
  if(model == "fullREDM"){
    overdisp <- x[[1]]$lambda
    if(opt$dataset_generation %in% c("GenerationMGnorm", "generationMGnorm")){
      overdisp <- c(NA,NA)
    }
  }else if(model=="fullREDMsinglelambda"){
    if(!is.na(x[[1]]$lambda)){
      stopifnot(x[[1]]$lambda[1] == x[[1]]$lambda[2])
    }
    overdisp <- x[[1]]$lambda[1]
  }
  
  if(covs_true){
    ## we have covariances
    covs <- covs
  }else{
    ## if we don't have covariances, we put zeros in the off-diagonal elements of the covariance matrix
    covs <- rep(0,((x[[1]]$d-1)**2-(x[[1]]$d-1))/2)
  }
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                sds, ##sd RE ## this is particular to GenerationCnorm
                covs, ## covariances RE
                overdisp)
}else if(model %in% c("diagREDM", "diagREDMsinglelambda" )){
  if(model == "diagREDM"){
    if(opt$dataset_generation %in% c("GenerationMGnorm", "generationMGnorm")){
      overdisp = c(NA, NA)
    }else{
      overdisp <- x[[1]]$lambda
    }
  }else if(model=="diagREDMsinglelambda"){
    stopifnot(x[[1]]$lambda[1] == x[[1]]$lambda[2])
    overdisp <- x[[1]]$lambda[1]
  }
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                sds, ##sd RE ## this is particular to GenerationCnorm
                overdisp)
}else if(model=="singleREDM"){
  if(opt$dataset_generation %in% c("GenerationMGnorm", "generationMGnorm")){
    overdisp = c(NA, NA)
  }else{
    overdisp <- x[[1]]$lambda
  }
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                NA, ## single logsd for the single random effects. They are not simulated like this in most datasets, so I add NA. If I am using a generation where random effects are simulated with a single value, I should add it here
                overdisp)
}else if(model %in% c("FEDMsinglelambda" )){
  overdisp <- x[[1]]$lambda
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                overdisp)
}else if(model == "fullREDMonefixedlambda"){
  true_vals <- c(as.vector((x[[1]]$beta)),
                 rep(NA, length(python_like_select_name(runs[[1]]$par.fixed, 'cov_par_RE'))),
                 sds, c((x[[1]]$lambda)))
}else if(model == "fullREDMonefixedlambda3"){
  true_vals <- c(as.vector((x[[1]]$beta)),
                 rep(NA, length(python_like_select_name(runs[[1]]$par.fixed, 'cov_par_RE'))),
                 sds, c((x[[1]]$lambda)))
}else{
  stop('<Specify correct model>')
}

summaries_melt$true = rep(true_vals, length(lst))
summaries_melt$subtract = summaries_melt$value - summaries_melt$true

plot(summaries_melt$true[summaries_melt$Var1 == 'beta'],
     summaries_melt$value[summaries_melt$Var1 == 'beta']) ### correlation between true and estimated beta
abline(coef = c(0,1), lty='dashed', col='blue')

ggplot(cbind.data.frame(summaries_melt[summaries_melt$Var1 == 'beta',], intslope=c('Intercept', 'Slope') ),
       aes(x=true, y=value, shape=intslope))+geom_point()+
  geom_abline(coef = c(0,1), lty='dashed', col='black')+theme_bw()+labs(shape="", x='True beta coefficient', y='Estimated beta coefficient')+
  theme(legend.position = "bottom")#+scale_color_manual(values = c(''))
ggsave(paste0(first_part_output, "betacorrelation.pdf"), width = 2.5, height = 2.8)

ggplot(summaries_melt, aes(x=idx_param, y=subtract, group=idx_param))+
  geom_abline(slope = 0, intercept = 0, lty='dashed', col='blue')+
  geom_boxplot()+facet_wrap(~Var1, scales = "free", nrow=1)+labs(x="Parameter", y="Bias")+
  theme_bw()
ggsave(paste0(first_part_output, "bias.pdf"), width = 10, height = 3.5)
# summaries_mat = matrix(summaries_melt$value, nrow=sum(summaries_melt$Var2 == 1))

ggplot(cbind.data.frame(summaries_melt[summaries_melt$Var1 == 'beta',],  intslope=c('Intercept', 'Slope')),
       aes(x=idx_param, y=subtract, group=idx_param, lty=intslope))+
  geom_abline(slope = 0, intercept = 0, lty='dashed', col='blue')+
  geom_boxplot()+#+facet_wrap(~Var1, scales = "free", nrow=1)+
  labs(x="Beta", y="Bias", lty="")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+theme(legend.position = "bottom")#+scale_color_manual(values = c('#c2b2b2', '#72081d'))
ggsave(paste0(first_part_output, "bias_betas.pdf"), width = 2.5, height = 2.8)

ggplot(cbind.data.frame(summaries_melt[summaries_melt$Var1 == 'beta',],  intslope=c('Intercept', 'Slope')),
       aes(x=(1+idx_param) %/% 2, y=subtract, group=idx_param))+
  geom_abline(slope = 0, intercept = 0, lty='dashed', col='blue')+
  geom_boxplot()+#+facet_wrap(~Var1, scales = "free", nrow=1)+
  labs(x="Beta", y="Bias", lty="")+
  theme_bw()+
  facet_wrap(.~intslope, drop = T)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+theme(legend.position = "bottom")#+scale_color_manual(values = c('#c2b2b2', '#72081d'))
ggsave(paste0(first_part_output, "bias_betas_v2.pdf"), width = 2.2, height = 2)

ggplot(cbind.data.frame(summaries_melt[summaries_melt$Var1 == 'beta',],  intslope=c('Intercept', 'Slope')),
       aes(x=(1+idx_param) %/% 2, y=subtract, group=idx_param))+
  geom_abline(slope = 0, intercept = 0, lty='dashed', col='blue')+
  geom_boxplot()+#+facet_wrap(~Var1, scales = "free", nrow=1)+
  labs(x="Beta", y="Bias", lty="")+
  theme_bw()+
  facet_wrap(.~intslope, ncol = 1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+theme(legend.position = "bottom")#+scale_color_manual(values = c('#c2b2b2', '#72081d'))
ggsave(paste0(first_part_output, "bias_betas_v2_vert.pdf"), width = 2, height = 3)


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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_bw()
ggsave(paste0(first_part_output, "coverage.pdf"), width = 10, height = 3.5)

df_coverage$type_param[df_coverage$type_param == "beta_intercept"] <- "Intercept"
df_coverage$type_param[df_coverage$type_param == "beta_slope"] <- "Slope"
ggplot(df_coverage[df_coverage$type_param %in% c("Intercept",  "Slope"),],
       aes(x=parameter, y=CI, group=1, shape=type_param))+geom_abline(slope = 0, intercept = 0.95, col='blue', lty='dashed')+
  geom_line()+geom_point()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_bw()+lims(y=c(0,1))+
  theme(legend.position = "bottom")+
  labs(shape="", x='Beta', y='Coverage')+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(paste0(first_part_output, "coverage_beta.pdf"), width = 2.5, height = 2.8)

draws_binom <- rbinom(n = 1000, size = sum(all_pd_list), prob = 0.95)
hist(draws_binom)
var(draws_binom)
pbinomprobs <- data.frame(t(quantile(draws_binom, probs = c(0.025, 0.975))/ sum(all_pd_list)))
colnames(pbinomprobs) <- c('x1', 'x2')


df_coverage[df_coverage$parameter == 'beta','parameter'] <- 'beta.0'

pbinomprobs <- data.frame(t(apply(df_coverage[df_coverage$type_param %in% c("Intercept",  "Slope"),c('parameter', 'type_param')],
                                  1, function(i) unlist(c(i, pbinomprobs)))))
colnames(pbinomprobs) <- c('parameter', 'type_param', 'x1', 'x2')
pbinomprobs$x1 <- as.numeric(pbinomprobs$x1)
pbinomprobs$x2 <- as.numeric(pbinomprobs$x2)

transformparam <- function(parameter) 1+as.numeric(gsub('.*[.]', '', parameter)) %/% 2

if(min(df_coverage[grepl('beta', df_coverage$parameter),'CI']) < 0.3){
  min_y = min(df_coverage[grepl('beta', df_coverage$parameter),'CI'])
  min_y_2 = min(df_coverage[grepl('beta', df_coverage$parameter),'CI'])
}else{
  min_y = 0.3
  min_y_2 = 0.7
}

library(cowplot)
A <- ggplot(df_coverage[df_coverage$type_param %in% c("Intercept"),])+
  geom_ribbon(data = pbinomprobs[pbinomprobs$type_param %in% c("Intercept"),], aes(x=transformparam(parameter), ymin=x1, ymax=x2, group=type_param),
              fill='yellow', alpha=0.9)+
  geom_abline(slope = 0, intercept = 0.95, col='blue', lty='dashed')+
  geom_line(aes(x=(1+as.numeric(gsub('.*[.]', '', parameter)) %/% 2), y=CI,
                group=1))+geom_point(aes(x=transformparam(parameter), y=CI,
                                         group=1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_bw()+#lims(y=c(0,1))+
  theme(legend.position = "bottom")+
  facet_wrap(.~type_param)+
  labs(shape="", x='Beta', y='Coverage')+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylim(min=c(min_y, 1))
if(min_y == 0.3){
  A <- A+annotate(geom="label",x=-Inf,y=0.315,label="   //", fill="white", col='black', label.size = NA)+
    annotate(geom="label",x=Inf,y=0.315,label="//   ", fill="white", col='black', label.size = NA)+
    geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=0.31, ymax=0.32), col='black', fill="white", alpha=0.2, lty=1)
}
B <- ggplot(df_coverage[df_coverage$type_param %in% c("Slope"),])+
  geom_ribbon(data = pbinomprobs[pbinomprobs$type_param %in% c("Slope"),], aes(x=transformparam(parameter), ymin=x1, ymax=x2, group=type_param),
              fill='yellow', alpha=0.9)+
  geom_abline(slope = 0, intercept = 0.95, col='blue', lty='dashed')+
  geom_line(aes(x=(1+as.numeric(gsub('.*[.]', '', parameter)) %/% 2), y=CI,
                group=1))+geom_point(aes(x=transformparam(parameter), y=CI,
                                         group=1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_bw()+#lims(y=c(0,1))+
  theme(legend.position = "bottom")+
  scale_y_continuous(position = "right", limits = c(min_y_2, 1))+
  facet_wrap(.~type_param)+
  labs(shape="", x='Beta', y='Coverage')+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=element_blank())
if(min_y == 0.3){
  B <- B+  annotate(geom="label",x=-Inf,y=0.715,label="   //", fill="white", col='black', label.size = NA)+
    annotate(geom="label",x=Inf,y=0.715,label="//   ", fill="white", col='black', label.size = NA)+
    geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=0.71, ymax=0.715), col='black', fill="white", alpha=0.2, lty=1)
}


pdf(paste0(first_part_output, "coverage_beta_v3.pdf"), width = 2.2, height = 2)
plot_grid(A, B, nrow=1, rel_widths = c(1.1, 0.95))
dev.off()

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

cat('Comparing converged and non-converged')

if((length(runs_nonconverged_list) > 0) & (length(runs) > 0) ){
  comparison_conv_notconv <- rbind(cbind.data.frame(melt(runs_nonconverged_list), conv='not_conv'),
                                   cbind.data.frame(melt(lapply(runs, function(j) try(j$par.fixed))), conv='converged'))
  comparison_conv_notconv$param_name = make.names(names(runs[[1]]$par.fixed), unique = T)
  
  comparison_conv_notconv$param_type <- gsub("\\..*", "", comparison_conv_notconv$param_name)
  p <- ggplot(comparison_conv_notconv, aes(x=interaction(conv, param_name), y=value, col=conv))+
    geom_boxplot()+facet_wrap(.~param_type, scales = "free", nrow=1)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme_bw()
  
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
  
  pdf(paste0(first_part_output, "comparison_nonconv.pdf"), width = 10, height = 3.5)
  grid::grid.draw(gp)
  dev.off()
  ggplot(comparison_conv_notconv, aes(x=interaction(conv, param_name), y=value, col=conv))+
    geom_violin()+geom_point(alpha=0.1)+facet_wrap(.~param_name, scales = "free", nrow=3)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+theme_bw()
  ggsave(paste0(first_part_output, "comparison_nonconv2.pdf"), width = 15, height = 5.5)
  
  cat('Reading beta intercept\n')
  true_beta_intercept <- readRDS(paste0("../data/assessing_models_simulation/additional_files/multiple_fixed_", idx_dataset_betaintercept, ".RDS"))
  print(true_beta_intercept)
  cat('Reading beta slope\n')
  true_beta_slope <- readRDS(paste0("../data/assessing_models_simulation/additional_files/multiple_fixed_", idx_dataset_betaslope, ".RDS"))
  print(true_beta_slope)
  boxplot(python_like_select_rownames(sapply(runs, function(i) i$par.fixed), 'beta')) ## THIS IS NOT THE BIAS!
  boxplot(t(python_like_select_rownames(sapply(runs, function(i) i$par.fixed), 'beta')))
  
  betas <- python_like_select_rownames(sapply(runs, function(i) i$par.fixed), 'beta')
  
  print(dim(betas))
  
  print('True betas')
  print(rbind(true_beta_intercept, true_beta_slope))
  
  print(betas)
  
  betas_subtract <- sweep(betas, 1, as.vector(rbind(true_beta_intercept, true_beta_slope)), '-')
  boxplot(betas_subtract)
  
  true_beta_intercept
  true_beta_slope
  
  x[[1]]$beta
  matrix(python_like_select_name(runs[[1]]$par.fixed, 'beta'), nrow=2)
  
  summaries_melt[summaries_melt$idx_param == 1,]
  
  ## check any potential correlations between estimates, for identifiability
  
  all_pars <- do.call('cbind', sapply(runs, `[`, 'par.fixed'))
  if(dim(all_pars)[1] < 10){
    pairs(t(all_pars))
    ggpairs(data.frame(t(all_pars)), aes(alpha = 0.4))
    ggsave(paste0(first_part_output, "correlationestimates.pdf"), width = 15, height = 15)
    ##
    model
  }
}
