#-------------------------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
library(parallel) ## mclapply
library(latex2exp)
library(dplyr)
library(ggplot2)
library(reshape2)
library(xtable)
library(scales) ## for alpha transparency in plots
source("helper_functions.R")
source("../../2_inference_TMB/helper_TMB.R")
set.seed(1245)
#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#

set.seed(234)

softmax = function(x){
  if(is.null(dim(x))){
    ## vector
    sum_x = sum(exp(x))
    exp(x)/sum_x
  }else{
    ## matrix
    sum_x = rowSums(exp(x))
    sweep(exp(x), 1, sum_x, '/')
  }
}

reflect_matrix = function(m){
  m[nrow(m):1,]
}

compile("FE_multinomial.cpp", "-std=gnu++17")
# compile("sea.cpp", "-std=gnu++17")
dyn.load(dynlib("FE_multinomial"))
# set.seed(123)

give_data_and_parameters = function(n, d){
  beta = rbind(rgamma(d-1, 1, 1),rgamma(d-1, 6, 2))
  # beta = rbind(rgamma(d-1, 1, 1), 0)
  # beta[2,d-1] = 0 ## this is set to zero as it is redundant. It is not redundant!
  x = cbind(1, sample(c(0,1), n, replace = TRUE))
  # x = cbind(1, sample(c(0), n, replace = TRUE))
  # x = x[order(x[,2]),]
  true_theta = softmax(cbind(x %*% beta + matrix(rnorm(n*(d-1), mean = 0, sd=0.1), ncol=(d-1)), 0))
  # image(t(true_theta)) 
  # pheatmap::pheatmap(true_theta)
  data <- list(Y = t(apply(true_theta, 1, function(p) rmultinom(1, 300, p))), n=n, d=d, x=x)
  # dev.off(); image(t(data$Y))
  
  return(list(data=data, true_theta=true_theta, beta=beta, x=x))
  
}

give_data_and_parameters_complex = function(n, d, betashape_intersect, betashape_slope){
  ## slightly more complicated version of give_data_and_parameters
  beta = rbind(rgamma(n = d-1, shape = betashape_intersect, rate = betashape_intersect),
               rgamma(n = d-1, shape = betashape_slope, rate = betashape_slope))
  x = cbind(1, sample(c(0,1), n, replace = TRUE))
  true_theta = softmax(cbind(x %*% beta + matrix(rnorm(n*(d-1), mean = 0, sd=0.1), ncol=(d-1)), 0))
  data <- list(Y = t(apply(true_theta, 1, function(p) rmultinom(1, 300, p))), n=n, d=d, x=x)

  return(list(data=data, true_theta=true_theta, beta=beta, x=x))
  
}

give_data_and_parameters_give_coefs = function(matrix_beta, matrix_x){
  d = ncol(matrix_beta) + 1
  n = nrow(matrix_x)
  true_theta = softmax(cbind(matrix_x %*% matrix_beta + matrix(rnorm(n*(d-1), mean = 0, sd=0.1), ncol=(d-1)), 0))
  data <- list(Y = t(apply(true_theta, 1, function(p) rmultinom(1, 300, p))), n=n, d=d, x=matrix_x)
  return(list(data=data, true_theta=true_theta, beta=matrix_beta, x=matrix_x))
}

inference_random_start = function(bool_identical_initial_beta, d, data, x, true_beta){
  
  if(bool_identical_initial_beta){
    parameters <- list( beta = array(rep(runif(1, min = -4, max = 4), 2*(d-1)), dim = c(2,d-1)) )
  }else{
    parameters <- list(
      # beta = array(rbind(rep(1,d-1), rep(NA,d-1)), dim = c(2,d-1))
      # beta = array(rep(runif(1, min = -4, max = 4), 2*(d-1)), dim = c(2,d-1))
      beta = array(c(rep(runif(1, min = -4, max = 4), (d-1)),
                     rep(runif(1, min = -4, max = 4), (d-1))),
                   dim = c(2,d-1))
    )
  }
  
  obj <- MakeADFun(data, parameters, DLL="FE_multinomial")
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  obj$he()    ## <-- Analytical hessian
  sdreport(obj)
  
  # betas_estimate = rbind(opt$par[grep('beta_intersect', names(opt$par))],
  #                        c(opt$par[grep('beta_slope', names(opt$par))]))
  betas_estimate = matrix(opt$par, nrow=2)
  # betas_estimate
  # parameters$beta
  # beta
  
  return(list(x=x, betas_estimate = matrix(opt$par, nrow=2),
              true_beta = true_beta,
              initial_beta = parameters$beta,
              stderr = give_stderr(sdreport(obj)),
              tmb_object=sdreport(obj)))
}

replicate_inference = function(n=100, d=2, nreplicas = 30){
  simulated_data = give_data_and_parameters(n = n, d = d)
  inference = c(replicate(nreplicas, inference_random_start(F, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE),
                replicate(nreplicas, inference_random_start(T, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE))
  return(inference)
}

replicate_inference_complex = function(n=100, d=2, nreplicas = 30, betashape_intersect, betashape_slope){
  # wrapper for give_data_and_parameters_complex and its inference
  simulated_data = give_data_and_parameters_complex(n = n, d = d, betashape_intersect, betashape_slope)
  inference = c(replicate(nreplicas, inference_random_start(F, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE),
                replicate(nreplicas, inference_random_start(T, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE))
  return(inference)
}

replicate_inference_varysinglecoef = function(n=100, d=2, nreplicas = 30, betashape_intersect, betashape_slope){
  # wrapper for creating datasets and inferring parameters in which only one parameter varies
  x = cbind(1, sample(c(0,1), n, replace = TRUE))
  print(x)
  simulated_data = give_data_and_parameters_give_coefs(matrix_beta = rbind(betashape_intersect, betashape_slope),
                                                       matrix_x = x)
  inference = c(replicate(nreplicas, inference_random_start(F, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE))
  return(inference)
}

create_plot = function(n=100, d=2, simulated_data=NULL){
  if(is.null(simulated_data)){
    ## if no simulated data is provided, I generate random data here
    simulated_data = give_data_and_parameters(n = n, d = d)
  }
  
  nreplicas = 30
  inference = c(replicate(nreplicas, inference_random_start(F, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE),
                replicate(nreplicas, inference_random_start(T, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE))
  
  lims = lapply(1:2, function(beta_idx) sapply(list(c(sapply(inference, function(i) i$betas_estimate)[beta_idx,],
                       sapply(inference, function(i) i$initial_beta)[beta_idx,],
                       sapply(inference, function(i) i$true_beta)[beta_idx,])), function(j) c(min(j), max(j))))
  plot(t(sapply(inference, function(i) i$betas_estimate)), col=alpha(factor(rep(c(1,2), each=nreplicas)), 0.2),
       pch='E', xlab='Beta intercept', ylab='Beta slope', xlim=lims[[1]], ylim=lims[[2]])
  points(t(sapply(inference, function(i) i$initial_beta)), col=alpha(factor(rep(c(1,2), each=nreplicas)), 0.2), pch=1, cex=0.5)
  points(t(sapply(inference, function(i) i$true_beta)), col='blue', pch='T')
  
  return(inference)
}

#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#

pdf("../../../../CDA_in_Cancer/results/ProjectDA/TMB/bias_beta_FE_multinomial.pdf",
    width = 10)
par(mfrow=c(2,3))
results_inference = replicate(6, create_plot())
dev.off()
## two examples of sets of coefficients whihc should give equally good fits

## I only select the first one from each dataset run since the estimates are the same for any initial value
plot(t(sapply(1:ncol(results_inference), function(dataset_idx) c(results_inference[2,dataset_idx][[1]]$betas_estimate[1,],
                                                               results_inference[2,dataset_idx][[1]]$true_beta[1,]))),
     xlab = "True intercept", ylab='Estimated intercept', main='There is a bias in underestimating\n the intercept')
abline(coef = c(0,1), col='blue', lty='dashed')


#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#
## Now width d>2
results_inference_multivariate = replicate(50, replicate_inference(n = 100, d = 4))
results_inference_multivariate

par(mfrow=c(1,2))
plot(t(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) c(results_inference_multivariate[2,dataset_idx][[1]]$betas_estimate[1,1],
                                                                              results_inference_multivariate[2,dataset_idx][[1]]$true_beta[1,1]))),
     xlab = "True intercept", ylab='Estimated intercept', main='Intercept is not biased\n in the multivariate case\n (one logR)')
abline(coef = c(0,1), col='blue', lty='dashed')

plot(t(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) c(results_inference_multivariate[2,dataset_idx][[1]]$betas_estimate[2,1],
                                                                              results_inference_multivariate[2,dataset_idx][[1]]$true_beta[2,1]))),
     xlab = "True intercept", ylab='Estimated intercept', main='Slope is well estimated\n in the multivariate case\n (one logR)')
abline(coef = c(0,1), col='blue', lty='dashed')

## pairs plots
# for intercept
intercepts_multivariate = t(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) c(results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[1,],
                                                                         results_inference_multivariate[1,dataset_idx][[1]]$true_beta[1,])))
colnames(intercepts_multivariate) = apply(expand.grid(paste0('beta intercept ', 1:(ncol(intercepts_multivariate)/2)), c(' estimate', ' true')), 1, paste0, collapse="")
pairs(intercepts_multivariate)
plot_pairs_with_identity(intercepts_multivariate, pch=19)

## Compute the bias (if any - there isn't)
#dev.off()
par(mfrow=c(2,3))
# intercept
apply(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) results_inference_multivariate[1,dataset_idx][[1]]$true_beta[1,] - results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[1,]),
      1, hist)
# slope
apply(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) results_inference_multivariate[1,dataset_idx][[1]]$true_beta[2,] - results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[2,]),
      1, hist)

bias_intercept = apply(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) results_inference_multivariate[1,dataset_idx][[1]]$true_beta[1,] - results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[1,]),
                   1, mean)
bias_slope = apply(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) results_inference_multivariate[1,dataset_idx][[1]]$true_beta[2,] - results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[2,]),
                       1, mean)
bias_intercept ## essentially zero: no bias
bias_slope ## essentially zero: no bias
#-------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------#

## What happens if there are zero values, either in the intersect or the slope?
## (in the case of more complicated models, these zero betas are over or under-estimated, consistently)

#-------------------------------------------------------------------------------------------------------------------#
## In the intersect
give_data_and_parameters_complex(n = 100, d = 5, betashape_intersect = 0, betashape_slope = 1)

results_inference_multivariate_zeros = replicate(50, replicate_inference_complex(n = 100, d = 4,
                                                                                 betashape_intersect = 0,
                                                                                 betashape_slope = 1))
par(mfrow=c(1,2), mar=c(2,2,2,2))
plot(t(sapply(1:ncol(results_inference_multivariate_zeros), function(dataset_idx) c(results_inference_multivariate_zeros[2,dataset_idx][[1]]$betas_estimate[1,1],
                                                                                    results_inference_multivariate_zeros[2,dataset_idx][[1]]$true_beta[1,1]))),
     xlab = "True intercept", ylab='Estimated intercept', main='Intercept is not biased\n in the multivariate case\n (one logR)')
abline(coef = c(0,1), col='blue', lty='dashed')

plot(t(sapply(1:ncol(results_inference_multivariate_zeros), function(dataset_idx) c(results_inference_multivariate_zeros[2,dataset_idx][[1]]$betas_estimate[2,1],
                                                                              results_inference_multivariate_zeros[2,dataset_idx][[1]]$true_beta[2,1]))),
     xlab = "True intercept", ylab='Estimated intercept', main='Slope is well estimated\n in the multivariate case\n (one logR)')
abline(coef = c(0,1), col='blue', lty='dashed')

intercepts_multivariate_zeros = t(sapply(1:ncol(results_inference_multivariate_zeros), function(dataset_idx) c(results_inference_multivariate_zeros[1,dataset_idx][[1]]$betas_estimate[1,],
                                                                                                   results_inference_multivariate_zeros[1,dataset_idx][[1]]$true_beta[1,])))
colnames(intercepts_multivariate) = apply(expand.grid(paste0('beta intercept ', 1:(ncol(intercepts_multivariate)/2)), c(' estimate', ' true')), 1, paste0, collapse="")
pairs(intercepts_multivariate)
plot_pairs_with_identity(intercepts_multivariate, pch=19)

## all betas; slope and not
all_betas_zero_intercept = rbind(t(sapply(1:ncol(results_inference_multivariate_zeros), function(dataset_idx) c(results_inference_multivariate_zeros[2,dataset_idx][[1]]$betas_estimate[1,1],
                                                                                     results_inference_multivariate_zeros[2,dataset_idx][[1]]$true_beta[1,1]))),
      t(sapply(1:ncol(results_inference_multivariate_zeros), function(dataset_idx) c(results_inference_multivariate_zeros[2,dataset_idx][[1]]$betas_estimate[2,1],
                                                                                     results_inference_multivariate_zeros[2,dataset_idx][[1]]$true_beta[2,1]))))
pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaintercept.pdf", height = 4)
par(mfrow=c(1,1), mar=c(4,4,4,4))
plot(all_betas_zero_intercept,
     xlab='Estimated beta values (intercept and slope pooled)',
     ylab='True beta values (intercept and slope pooled)',
     col=factor(rep(c(0,1), rep(nrow(all_betas_zero_intercept)/2, 2))), pch=18)
abline(coef = c(0,1), lty='dashed')
dev.off()

## clearly the zero intercepts are something over- and under-estimated. Does this correspond to an inversely
## (under- or over-estimated) beta for the corresponding slopes?

pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaintercept_correlation.pdf", width=3, height = 4)
plot((all_betas_zero_intercept[,1] - all_betas_zero_intercept[,2])[1:(nrow(all_betas_zero_intercept)/2)],
(all_betas_zero_intercept[,1] - all_betas_zero_intercept[,2])[(1+nrow(all_betas_zero_intercept)/2):nrow(all_betas_zero_intercept)],
xlab='Signed error in intercept recovery', ylab='Signed error in slope recovery', pch=18)
abline(coef=c(0,-1), lty='dashed')
dev.off()

## there is a slight negative correlation, which corresponds to, indeed, higher overestimates to a beta parameter
## when its intercept/slope counterpart is underestimated

#-------------------------------------------------------------------------------------------------------------------#
# Zeros in the slope
results_inference_multivariate_zerosslope = replicate(50, replicate_inference_complex(n = 100, d = 4,
                                                                                 betashape_intersect = 1,
                                                                                 betashape_slope = 0))

## here I am only selecting the first log-ratio slopes (betas_estimate[x,1])
all_betas_zero_slope = rbind(t(sapply(1:ncol(results_inference_multivariate_zerosslope), function(dataset_idx) c(results_inference_multivariate_zerosslope[2,dataset_idx][[1]]$betas_estimate[1,1],
                                                                                                                 results_inference_multivariate_zerosslope[2,dataset_idx][[1]]$true_beta[1,1]))),
                                 t(sapply(1:ncol(results_inference_multivariate_zerosslope), function(dataset_idx) c(results_inference_multivariate_zerosslope[2,dataset_idx][[1]]$betas_estimate[2,1],
                                                                                                                     results_inference_multivariate_zerosslope[2,dataset_idx][[1]]$true_beta[2,1]))))
pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope.pdf", height = 4)
par(mfrow=c(1,1), mar=c(4,4,4,4))
plot(all_betas_zero_slope,
     xlab='Estimated beta values (intercept and slope pooled)',
     ylab='True beta values (intercept and slope pooled)',
     col=factor(rep(c(0,1), rep(nrow(all_betas_zero_slope)/2, 2))), pch=18)
abline(coef = c(0,1), lty='dashed')
dev.off()

pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation.pdf", width=3, height = 4)
plot((all_betas_zero_slope[,1] - all_betas_zero_slope[,2])[1:(nrow(all_betas_zero_slope)/2)],
     (all_betas_zero_slope[,1] - all_betas_zero_slope[,2])[(1+nrow(all_betas_zero_slope)/2):nrow(all_betas_zero_slope)],
     xlab='Signed error in intercept recovery', ylab='Signed error in slope recovery', pch=18)
abline(coef=c(0,-1), lty='dashed')
dev.off()

## Is the error in the recovery of zeros worse than in the recovery of non-zeros?
L1_results_inference_multivariate_zerosslope = do.call('cbind', lapply(results_inference_multivariate_zerosslope[1,],
                        function(i){
                          i$betas_estimate - results_inference_multivariate_zerosslope[1,1][[1]]$true_beta
                          }))
plot(density(abs(L1_results_inference_multivariate_zerosslope[1,])))
lines(density(abs(L1_results_inference_multivariate_zerosslope[2,])), col='blue')

#-------------------------------------------------------------------------------------------------------------------#
# Recovering the same simulated dataset, multiple times
# by default, in the simulated datasets we have thirty runs of inference (with different initial values)
# for each simulated dataset, but for the analysis above we have only used the first run (i.e. in $betas_estimate[1,1]
# and its equivalents, we always use the first column) FALSE! IS IT results_inference_multivariate_zerosslope[[XX]] instead?

# Here I look a at a single dataset, and all the runs of inference
# the index 2 refers to the second run, in the case below
dataset_idx_single = 1
results_inference_multivariate_zerosslope[2,dataset_idx_single][[1]]$betas_estimate[1,]
results_inference_multivariate_zerosslope[2,dataset_idx_single][[1]]$true_beta[1,]

## we have 30 replicates (30 runs) for random initial values and 30 runs for equal initial values
## all these 60 runs are good for assessing how the estimates change from run to run
dim(results_inference_multivariate_zerosslope)
single_dataset_several_runs = lapply(1:dim(results_inference_multivariate_zerosslope)[2], function(dataset_idx_single){
  cbind.data.frame(beta_intercept_est = as.vector(sapply(1:dim(results_inference_multivariate_zerosslope)[1], function(idx_replicate){
  results_inference_multivariate_zerosslope[idx_replicate, dataset_idx_single][[1]]$betas_estimate[1,]
})),
beta_intercept_true = as.vector(sapply(1:dim(results_inference_multivariate_zerosslope)[1], function(idx_replicate){
  results_inference_multivariate_zerosslope[idx_replicate, dataset_idx_single][[1]]$true_beta[1,]
})),
beta_slope_est = as.vector(sapply(1:dim(results_inference_multivariate_zerosslope)[1], function(idx_replicate){
  results_inference_multivariate_zerosslope[idx_replicate, dataset_idx_single][[1]]$betas_estimate[2,]
})),
beta_slope_true = as.vector(sapply(1:dim(results_inference_multivariate_zerosslope)[1], function(idx_replicate){
  results_inference_multivariate_zerosslope[idx_replicate, dataset_idx_single][[1]]$betas_estimate[1,]
})))})

# with [c(T,F,F)] I only select one of the betas at a time
plot(single_dataset_several_runs[[1]]$beta_intercept_est[c(T,F,F)], single_dataset_several_runs[[1]]$beta_intercept_true[c(T,F,F)])
abline(v = single_dataset_several_runs[[1]]$beta_intercept_true[c(T,F,F)][1])
hist(single_dataset_several_runs[[1]]$beta_intercept_est[c(T,F,F)])

single_dataset_several_runs_concat = do.call('rbind', single_dataset_several_runs)

pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation_all_runs.pdf", width=3, height = 4)
plot(single_dataset_several_runs_concat$beta_intercept_est - single_dataset_several_runs_concat$beta_intercept_true,
     single_dataset_several_runs_concat$beta_slope_est - single_dataset_several_runs_concat$beta_slope_true)
dev.off()

#-------------------------------------------------------------------------------------------------------------------#
# Showing the standard error too
# since here the zero coefficients are the slope, I only select the slopes
# and I only select the first run, since they all give the same results
single_dataset_several_runs_concat_stderr = do.call('rbind', lapply(1:dim(results_inference_multivariate_zerosslope)[2], function(dataset_idx){
  cbind.data.frame(stderr_beta_slope=select_slope_2(results_inference_multivariate_zerosslope[2,dataset_idx][[1]]$stderr, v=F),
        est_beta_slope=results_inference_multivariate_zerosslope[2,dataset_idx][[1]]$betas_estimate[2,],
        true_beta_slope=results_inference_multivariate_zerosslope[2,dataset_idx][[1]]$true_beta[2,],
        idx_dataset=dataset_idx)
}))
select_slope_2(results_inference_multivariate_zerosslope[3,dataset_idx_single][[1]]$stderr, v=F)
results_inference_multivariate_zerosslope[2,dataset_idx_single][[1]]$betas_estimate[2,]

## why is there a bias? it should be centered around zero
pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_est_and_stderr.pdf", width=4, height = 4)
plot(0,0, xlim = c(0,dim(single_dataset_several_runs_concat_stderr)[1]),
     ylim=c(min(single_dataset_several_runs_concat_stderr$est_beta_slope-single_dataset_several_runs_concat_stderr$stderr_beta_slope),
            max(single_dataset_several_runs_concat_stderr$est_beta_slope+single_dataset_several_runs_concat_stderr$stderr_beta_slope)), col='white',
     xlab='Index of beta and dataset', ylab='Beta estimate CI')
segments(x0=1, x1=dataset_idx, y0=0, y1=0, lty='dashed')
for(dataset_idx in 1:dim(single_dataset_several_runs_concat_stderr)[1]){
  segments(x0 = dataset_idx, y0 = single_dataset_several_runs_concat_stderr[dataset_idx,]$est_beta_slope - single_dataset_several_runs_concat_stderr[dataset_idx,]$stderr_beta_slope,
         x1 = dataset_idx, y1 = single_dataset_several_runs_concat_stderr[dataset_idx,]$est_beta_slope + single_dataset_several_runs_concat_stderr[dataset_idx,]$stderr_beta_slope,
  col=as.factor(dataset_idx))
}
dev.off()

#-------------------------------------------------------------------------------------------------------------------#
# So, when does the problem become grave?
# Running with several beta_shapes

betashape_slopes_vec = rgamma(n = 60, rate = 0.3, shape = 0.1)
hist(betashape_slopes_vec, breaks = 40)
results_inference_multivariate_zerosslope_severalbeta = mclapply(betashape_slopes_vec, function(betashape_slope_idx){replicate_inference_complex(n = 100, d = 4,
                                                                                      betashape_intersect = 1,
                                                                                      betashape_slope = betashape_slope_idx, nreplicas = 2)})
all_betas_zero_slope_severalbeta = cbind(do.call('rbind', lapply(1:length(results_inference_multivariate_zerosslope_severalbeta), function(dataset_idx) cbind(results_inference_multivariate_zerosslope_severalbeta[[dataset_idx]][[1]]$betas_estimate[1,],
                                                                                                                                         results_inference_multivariate_zerosslope_severalbeta[[dataset_idx]][[1]]$true_beta[1,],
                                                                                                                                         results_inference_multivariate_zerosslope_severalbeta[[dataset_idx]][[1]]$stderr))),
                                         do.call('rbind', lapply(1:length(results_inference_multivariate_zerosslope_severalbeta), function(dataset_idx) cbind(results_inference_multivariate_zerosslope_severalbeta[[dataset_idx]][[1]]$betas_estimate[2,],
                                                                                                                             results_inference_multivariate_zerosslope_severalbeta[[dataset_idx]][[1]]$true_beta[2,],
                                                                                                                             results_inference_multivariate_zerosslope_severalbeta[[dataset_idx]][[1]]$stderr))))
colnames(all_betas_zero_slope_severalbeta) = c('Beta_intercept_est', 'Beta_intercept_true', 'Beta_intercept_stderr', 'Beta_slope_est', 'Beta_slope_true', 'Beta_slope_stderr')
all_betas_zero_slope_severalbeta = as.data.frame(all_betas_zero_slope_severalbeta)
# this can be done since the number of features is the same in all datasets
all_betas_zero_slope_severalbeta$idx_dataset = rep(1:length(results_inference_multivariate_zerosslope_severalbeta), each=dim(all_betas_zero_slope_severalbeta)[1]/length(results_inference_multivariate_zerosslope_severalbeta))

par(mfrow=c(1,2))
plot(all_betas_zero_slope_severalbeta[,'Beta_slope_est'],
     all_betas_zero_slope_severalbeta[,'Beta_slope_true'])
plot(all_betas_zero_slope_severalbeta[,'Beta_slope_est'],
     all_betas_zero_slope_severalbeta[,'Beta_slope_true'],
     ylim=c(-.5,.4), xlim=c(-.5,.4))
     # ylim=c(-.5,.4), xlim=c(-.5,.4))
abline(coef=c(0,1), lty='dashed'); abline(v=0, lty='dashed')
# zoomed in image


# I try to sort the dataframe all_betas_zero_slope_severalbeta in some way that can explain in which cases is the
# recovery successful and in which it's not.
# I need to sort it by the intersect, in some way or another

# group by mean intercept coefficient within a dataset
order_datasets = all_betas_zero_slope_severalbeta %>% group_by(all_betas_zero_slope_severalbeta$idx_dataset)  %>% mutate(mean=mean(Beta_intercept_true)) %>% dplyr::pull(mean) %>% order
# group by each intercept coefficient
order_datasets = order(all_betas_zero_slope_severalbeta$Beta_intercept_est)

all_betas_zero_slope_severalbeta = all_betas_zero_slope_severalbeta[order_datasets,]
cols = as.factor( (sapply(1:dim(all_betas_zero_slope_severalbeta)[1], function(dataset_idx) all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est - 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr) < 0 ) & (sapply(1:dim(all_betas_zero_slope_severalbeta)[1], function(dataset_idx) all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est + 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr) > 0))

## and now plot the standard errors
pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_est_and_stderr_severalbeta.pdf", width=6, height = 4)
plot(0,0, xlim = c(0,dim(all_betas_zero_slope_severalbeta)[1]),
     ylim=c(min(all_betas_zero_slope_severalbeta$Beta_slope_est-1.96*all_betas_zero_slope_severalbeta$Beta_slope_stderr),
            max(all_betas_zero_slope_severalbeta$Beta_slope_est+1.96*all_betas_zero_slope_severalbeta$Beta_slope_stderr)), col='white',
     xlab='Index of beta and dataset', ylab='Beta estimate CI')
for(dataset_idx in 1:dim(all_betas_zero_slope_severalbeta)[1]){
  segments(x0 = dataset_idx, y0 = all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est - 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr,
           x1 = dataset_idx, y1 = all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est + 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr,
           col=cols[dataset_idx])
}
dev.off()

# group by mean intercept coefficient for each beta independently
order_datasets = order(all_betas_zero_slope_severalbeta$Beta_intercept_true)
all_betas_zero_slope_severalbeta = all_betas_zero_slope_severalbeta[order_datasets,]
cols = as.factor( (sapply(1:dim(all_betas_zero_slope_severalbeta)[1], function(dataset_idx) all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est - 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr) < 0 ) & (sapply(1:dim(all_betas_zero_slope_severalbeta)[1], function(dataset_idx) all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est + 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr) > 0))
## and now plot the standard errors
pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_est_and_stderr_severalbeta_intercept_beta_independently.pdf", width=6, height = 4)
plot(0,0, xlim = c(0,dim(all_betas_zero_slope_severalbeta)[1]),
     ylim=c(min(all_betas_zero_slope_severalbeta$Beta_slope_est-1.96*all_betas_zero_slope_severalbeta$Beta_slope_stderr),
            max(all_betas_zero_slope_severalbeta$Beta_slope_est+1.96*all_betas_zero_slope_severalbeta$Beta_slope_stderr)), col='white',
     xlab='Index of beta and dataset', ylab='Beta estimate CI')
for(dataset_idx in 1:dim(all_betas_zero_slope_severalbeta)[1]){
  segments(x0 = dataset_idx, y0 = all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est - 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr,
           x1 = dataset_idx, y1 = all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est + 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr,
           col=cols[dataset_idx])
}
dev.off()

# group by mean slopes coefficient for each beta independently
order_datasets = order(all_betas_zero_slope_severalbeta$Beta_slope_true)
all_betas_zero_slope_severalbeta = all_betas_zero_slope_severalbeta[order_datasets,]
cols = as.factor( (sapply(1:dim(all_betas_zero_slope_severalbeta)[1], function(dataset_idx) all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est - 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr) < 0 ) & (sapply(1:dim(all_betas_zero_slope_severalbeta)[1], function(dataset_idx) all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est + 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr) > 0))
## and now plot the standard errors
pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_est_and_stderr_severalbeta_slope_beta_independently.pdf", width=6, height = 4)
plot(0,0, xlim = c(0,dim(all_betas_zero_slope_severalbeta)[1]),
     ylim=c(min(all_betas_zero_slope_severalbeta$Beta_slope_est-1.96*all_betas_zero_slope_severalbeta$Beta_slope_stderr),
            max(all_betas_zero_slope_severalbeta$Beta_slope_est+1.96*all_betas_zero_slope_severalbeta$Beta_slope_stderr)), col='white',
     xlab='Index of beta and dataset', ylab='Beta estimate CI')
for(dataset_idx in 1:dim(all_betas_zero_slope_severalbeta)[1]){
  segments(x0 = dataset_idx, y0 = all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est - 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr,
           x1 = dataset_idx, y1 = all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_est + 1.96*all_betas_zero_slope_severalbeta[dataset_idx,]$Beta_slope_stderr,
           col=cols[dataset_idx])
}
dev.off()


#-------------------------------------------------------------------------------------------------------------------#
# How does this problem depend on the number of categories?
set.seed(1345)
d_its = 2:5
results_inference_multivariate_zerosslope_d = sapply(d_its, function(d_it) replicate(500, replicate_inference_complex(n = 100, d = d_it,
                                                                                      betashape_intersect = 1,
                                                                                      betashape_slope = 0, nreplicas = 2)))
## save or read in
saveRDS(results_inference_multivariate_zerosslope_d, "../../../data/robjects_cache/tmb_results_simulations/FEmultinomial_results_inference_multivariate_zerosslope_d.RDS")
results_inference_multivariate_zerosslope_d = readRDS('../../../data/robjects_cache/tmb_results_simulations/FEmultinomial_results_inference_multivariate_zerosslope_d.RDS')

## each column of results_inference_multivariate_zerosslope_d contains a different d
dim(results_inference_multivariate_zerosslope_d)
sapply(1:dim(results_inference_multivariate_zerosslope_d)[2], function(col){
  ncol(results_inference_multivariate_zerosslope_d[1,col][[1]]$initial_beta)
})

all_betas_zero_slope_severald = lapply(1:(dim(results_inference_multivariate_zerosslope_d)[2]), function(d_it){
  lapply(1:(dim(results_inference_multivariate_zerosslope_d)[1]), function(it){
    results_inference_multivariate_zerosslope_d[it,d_it][[1]]$betas_estimate[2,] - results_inference_multivariate_zerosslope_d[it,d_it][[1]]$true_beta[2,]
})
})
all_betas_zero_slope_severald_melt = melt(all_betas_zero_slope_severald)
all_betas_zero_slope_severald_melt$L1 = d_its[all_betas_zero_slope_severald_melt$L1]

## remove clear outliers
all_betas_zero_slope_severald_melt_subsample = all_betas_zero_slope_severald_melt %>% group_by(L1) %>% mutate(top=quantile(abs(value), probs = c(0.99))) %>% ungroup() %>% mutate(outliers=(abs(value) > top)) %>% filter(!outliers)

## subsample to get the same number of observations for each d
all_betas_zero_slope_severald_melt_subsample = all_betas_zero_slope_severald_melt_subsample[as.vector(sapply(d_its, function(i) sample(x = which(all_betas_zero_slope_severald_melt$L1 == i), size = min(table(all_betas_zero_slope_severald_melt$L1)), replace = FALSE))),]

ggplot(all_betas_zero_slope_severald_melt_subsample, aes(x=L1, col=L2, group=L1, y=abs(value)))+geom_boxplot()
ggplot(all_betas_zero_slope_severald_melt_subsample, aes(x=L1, group=L1, y=abs(value)))+geom_violin()+geom_jitter(size=0.1)+theme_bw()+
  geom_line(data=all_betas_zero_slope_severald_melt_subsample %>% group_by(L1) %>% summarize(mean=mean(abs(value))),
            aes(x=L1, y=mean))+
  geom_line(data=all_betas_zero_slope_severald_melt_subsample %>% group_by(L1) %>% summarize(mean=mean(abs(value))),
            aes(x=L1, y=mean, group=1), col='blue')+
  geom_label(data=all_betas_zero_slope_severald_melt_subsample %>% group_by(L1) %>% summarize(mean=mean(abs(value))),
             aes(x=L1, y=mean, label=round(mean, 4)))+
  labs(y=latex2exp::TeX('L1 norm of error of zero beta slope'), x='Number of categories')
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation_depending_on_d3.pdf", width=4, height = 4)

ggplot(all_betas_zero_slope_severald_melt_subsample, aes(x=L1, group=L1, y=log(abs(value))))+geom_violin()+geom_jitter(size=0.1)+theme_bw()+
  geom_line(data=all_betas_zero_slope_severald_melt_subsample %>% group_by(L1) %>% summarize(mean=mean(abs(value))),
            aes(x=L1, y=log(mean)))+
  geom_line(data=all_betas_zero_slope_severald_melt_subsample %>% group_by(L1) %>% summarize(mean=mean(abs(value))),
            aes(x=L1, y=log(mean), group=1), col='blue')+
  geom_label(data=all_betas_zero_slope_severald_melt_subsample %>% group_by(L1) %>% summarize(mean=mean(abs(value))),
             aes(x=L1, y=log(mean), label=round(mean, 4)))+
  labs(y=latex2exp::TeX('L1 norm of error of zero beta slope'), x='Number of categories')
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation_depending_on_d3_log.pdf", width=4, height = 4)

## with the same data, put the fraction of betas for which zero is in the +/- 1.96std error of the estimate (CI)
## Every for rows correspond to the same dataset, as there are (two replicas) x (two different starting points) = 4 runs for each simulated dataset

zero_in_stderr_CI = function(vec){
  stopifnot(length(vec) == 2)
 (0 > min(vec)) & (0 < max(vec))
}

give_fraction_in_CI = function(run_result_list){
  if(typeof(run_result_list) == "logical"){
    return(NA)
  }else{
    ## for each beta slope that is zero, see if zero is included in +/- 1.96std error (CI)
    ## note that give_stderr by default only gives the stderr of beta slopes
    ## (in this case this is good enough because our zeros are in the beta slopes)
    return(sum(apply(rbind( run_result_list[[1]]$betas_estimate[2,] - 1.96*run_result_list[[1]]$stderr,
          run_result_list[[1]]$betas_estimate[2,] + 1.96*run_result_list[[1]]$stderr), 2, zero_in_stderr_CI)))
  }
}
results_inference_multivariate_zerosslope_d_outer_num_zero_in_cred = outer(1:ncol(results_inference_multivariate_zerosslope_d), 1:nrow(results_inference_multivariate_zerosslope_d),
                                                                           Vectorize(function(i,j) give_fraction_in_CI(results_inference_multivariate_zerosslope_d[j,i])))

all_betas_zero_slope_severald_melt_fractionCI = melt(t(results_inference_multivariate_zerosslope_d_outer_num_zero_in_cred)/(sapply(d_its, function(i) rep(i, nrow(results_inference_multivariate_zerosslope_d)))-1))
all_betas_zero_slope_severald_melt_fractionCI$Var2 = d_its[all_betas_zero_slope_severald_melt_fractionCI$Var2]

ggplot(all_betas_zero_slope_severald_melt_fractionCI,
       aes(x=Var2, y=value, group=Var2))+
  geom_violin()+geom_jitter(alpha=0.6, size=0.2)+
  geom_line(data=all_betas_zero_slope_severald_melt_fractionCI %>% group_by(Var2) %>% summarize(mean=mean(abs(value), na.rm = T)),
            aes(x=Var2, y=(mean), group=1), col='blue')+
  labs(y='Fraction of beta slopes with zero in \ntheir est CI', x='Number of categories')+theme_bw()
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation_depending_on_d3_fracinCI.pdf", width=4, height = 4)

## Plot the fraction of simulated datasets which have at least one beta that doesn't contain zero in their confidence interval,
## as we vary d
## if the values below are zero, it means that zero is included in all the betas
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset = melt((sapply(d_its, function(i) rep(i, nrow(results_inference_multivariate_zerosslope_d)))-1) - t(results_inference_multivariate_zerosslope_d_outer_num_zero_in_cred))
## if 1, there are some betas which don't contain 0
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset$value[all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset$value > 0] = 1
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset$Var2 = d_its[all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset$Var2]
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset$value = plyr::revalue(as.factor(all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset$value), c("0"="Zero in all betas' CI", "1"="Zero not in all betas's CI"))
ggplot()+
  geom_jitter(data=all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset,
                     aes(x=Var2,y=value))+
  geom_label(data=all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset %>% group_by(Var2) %>% summarize(percentage=table(value)[1]/n()) %>% mutate(value=NA),
             aes(x=Var2,label=percentage, y="Zero in all betas' CI"))+
  labs(x='d')
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation_depending_on_d3_fracinCI_per_dataset.pdf", width=4, height = 4)

## the fraction of betas the CI of which doesn't include zero is
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset
## the fraction of datasets with some beta the CI of which doesn't include zero is


ggplot(all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset)+geom_bar(aes(group=value,
                                                                                       y=factor(Var2, levels=rev(d_its)), fill=value), position='stack')+
  labs(y='d')+theme_bw()+theme(legend.position = "bottom", legend.title = element_blank())
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation_depending_on_d3_fracinCI_per_dataset_b.pdf", width=4, height = 4)


## Do the same but using other packages available in R
replicate_inference_complex_rpackages = function(n=100, d=2, nreplicas = 30, betashape_intersect, betashape_slope, give_stderr_only_slope=T){
  # wrapper for give_data_and_parameters_complex and its inference
  simulated_data = give_data_and_parameters_complex(n = n, d = d, betashape_intersect, betashape_slope)
  if(d == 2){
    ## binomial
    model0 <- glm(simulated_data$data$Y ~ simulated_data$data$x[,2], family="binomial")
    if(give_stderr_only_slope){
      stders = summary(model0)$coefficients[2, 'Std. Error']
    }else{
      stders = summary(model0)$coefficients[, 'Std. Error']
    }
    return(list(list(x=simulated_data$data$x, betas_estimate = matrix(matrix(model0$coefficients, nrow=2), nrow=2),
                true_beta = simulated_data$beta,
                initial_beta = NA,
                stderr = stders)))
  }else if(d>2){
    require(nnet)
    ## multinomial
    model0 <- nnet::multinom(simulated_data$data$Y ~ simulated_data$data$x[,2])
    if(give_stderr_only_slope){
      stders = summary(model0)$standard.errors[,2]
    }else{
      stders = t(summary(model0)$standard.errors)
    }
    return(list(list(x=simulated_data$data$x, betas_estimate = t(coefficients(model0)),
                true_beta = simulated_data$beta,
                initial_beta = NA,
                stderr = stders)))
  }else{
    stop('d has to be > 1')
  }
}

results_inference_multivariate_zerosslope_d_rpackages = sapply(d_its, function(d_it) replicate(100, replicate_inference_complex_rpackages(n = 100, d = d_it,
                                                                                                                      betashape_intersect = 1,
                                                                                                                      betashape_slope = 0, nreplicas = 2)))
saveRDS(results_inference_multivariate_zerosslope_d_rpackages, "../../../data/robjects_cache/tmb_results_simulations/FEmultinomial_results_inference_multivariate_zerosslope_d_rpackages.RDS")


results_inference_multivariate_zerosslope_d_outer_num_zero_in_cred_rpackages = outer(1:ncol(results_inference_multivariate_zerosslope_d_rpackages), 1:nrow(results_inference_multivariate_zerosslope_d_rpackages),
                                                                           Vectorize(function(i,j) give_fraction_in_CI(results_inference_multivariate_zerosslope_d_rpackages[j,i])))

all_betas_zero_slope_severald_rpackages = lapply(1:(dim(results_inference_multivariate_zerosslope_d_rpackages)[2]), function(d_it){
  lapply(1:(dim(results_inference_multivariate_zerosslope_d_rpackages)[1]), function(it){
    results_inference_multivariate_zerosslope_d_rpackages[it,d_it][[1]]$betas_estimate[2,] - results_inference_multivariate_zerosslope_d_rpackages[it,d_it][[1]]$true_beta[2,]
  })
})
all_betas_zero_slope_severald_melt_rpackages = melt(all_betas_zero_slope_severald_rpackages)
all_betas_zero_slope_severald_melt_rpackages$L1 = d_its[all_betas_zero_slope_severald_melt_rpackages$L1]

## remove clear outliers
all_betas_zero_slope_severald_melt_subsample_rpackages = all_betas_zero_slope_severald_melt_rpackages %>% group_by(L1) %>% mutate(top=quantile(abs(value), probs = c(0.99))) %>% ungroup() %>% mutate(outliers=(abs(value) > top)) %>% filter(!outliers)
## subsample to get the same number of observations for each d
all_betas_zero_slope_severald_melt_subsample_rpackages = all_betas_zero_slope_severald_melt_subsample_rpackages[as.vector(sapply(d_its, function(i) sample(x = which(all_betas_zero_slope_severald_melt_rpackages$L1 == i), size = min(table(all_betas_zero_slope_severald_melt_rpackages$L1)), replace = FALSE))),]

ggplot(all_betas_zero_slope_severald_melt_subsample_rpackages, aes(x=L1, group=L1, y=abs(value)))+geom_violin()+geom_jitter(size=0.1)+theme_bw()+
  geom_line(data=all_betas_zero_slope_severald_melt_subsample_rpackages %>% group_by(L1) %>% summarize(mean=mean(abs(value))),
            aes(x=L1, y=mean))+
  geom_line(data=all_betas_zero_slope_severald_melt_subsample_rpackages %>% group_by(L1) %>% summarize(mean=mean(abs(value))),
            aes(x=L1, y=mean, group=1), col='blue')+
  geom_label(data=all_betas_zero_slope_severald_melt_subsample_rpackages %>% group_by(L1) %>% summarize(mean=mean(abs(value))),
             aes(x=L1, y=mean, label=round(mean, 4)))+
  labs(y=latex2exp::TeX('L1 norm of error of zero beta slope'), x='Number of categories')
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation_depending_on_d_rpackages.pdf", width=4, height = 4)


all_betas_zero_slope_severald_melt_fractionCI_rpackages = melt(t(results_inference_multivariate_zerosslope_d_outer_num_zero_in_cred_rpackages)/(sapply(d_its, function(i) rep(i, ncol(results_inference_multivariate_zerosslope_d_outer_num_zero_in_cred_rpackages)))-1))
all_betas_zero_slope_severald_melt_fractionCI_rpackages$Var2 = d_its[all_betas_zero_slope_severald_melt_fractionCI_rpackages$Var2]

ggplot(all_betas_zero_slope_severald_melt_fractionCI_rpackages,
       aes(x=Var2, y=value, group=Var2))+
  geom_violin()+geom_jitter(alpha=0.6, size=0.2)+
  geom_line(data=all_betas_zero_slope_severald_melt_fractionCI_rpackages %>% group_by(Var2) %>% summarize(mean=mean(abs(value), na.rm = T)),
            aes(x=Var2, y=(mean), group=1), col='blue')+
  labs(y='Fraction of beta slopes with zero in \ntheir est CI', x='Number of categories')+theme_bw()
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation_depending_on_d3_fracinCI_rpckages.pdf", width=4, height = 4)

all_betas_zero_slope_severald_melt_fractionCI_both = melt(list(all_betas_zero_slope_severald_melt_fractionCI, all_betas_zero_slope_severald_melt_fractionCI_rpackages),id.vars=c('Var1', 'Var2', 'value'))
all_betas_zero_slope_severald_melt_fractionCI_both$L1 = c('TMB model', 'R packages')[all_betas_zero_slope_severald_melt_fractionCI_both$L1]

ggplot(all_betas_zero_slope_severald_melt_fractionCI_both,
       aes(x=interaction(L1,Var2), y=value, group=interaction(Var2, L1), col=L1))+
  scale_colour_manual(values = c("R packages" = "#ff9999", "TMB model" = "#4166f5"))+
  geom_violin()+geom_jitter(alpha=0.6, size=0.2)+
  geom_line(data=cbind.data.frame(all_betas_zero_slope_severald_melt_fractionCI %>% group_by(Var2) %>% summarize(mean=mean(abs(value), na.rm = T)), L1='TMB model'),
            aes(x=interaction(L1,Var2), y=(mean), group=1), col='blue')+
  geom_line(data=cbind.data.frame(all_betas_zero_slope_severald_melt_fractionCI_rpackages %>% group_by(Var2) %>% summarize(mean=mean(abs(value), na.rm = T)), L1='R packages'),
            aes(x=interaction(L1,Var2), y=(mean), group=1), col='red')+
  geom_point(data=cbind.data.frame(all_betas_zero_slope_severald_melt_fractionCI %>% group_by(Var2) %>% summarize(mean=mean(abs(value), na.rm = T)), L1='TMB model'),
            aes(x=interaction(L1,Var2), y=(mean), group=1), col='blue')+
  geom_point(data=cbind.data.frame(all_betas_zero_slope_severald_melt_fractionCI_rpackages %>% group_by(Var2) %>% summarize(mean=mean(abs(value), na.rm = T)), L1='R packages'),
            aes(x=interaction(L1,Var2), y=(mean), group=1), col='red')+
  labs(y='Fraction of beta slopes with zero in \ntheir CI', x='Number of categories')+theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_correlation_depending_on_d3_fracinCI_comparison.pdf", width=6, height = 5)


## Plot the fraction of simulated datasets which have at least one beta that doesn't contain zero in their confidence interval,
## as we vary d, for the R packages
## if the values below are zero, it means that zero is included in all the betas
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_rpackages = melt((sapply(d_its, function(i) rep(i, nrow(results_inference_multivariate_zerosslope_d_rpackages)))-1) - t(results_inference_multivariate_zerosslope_d_outer_num_zero_in_cred_rpackages))
## if 1, there are some betas which don't contain 0
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_rpackages$value[all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_rpackages$value > 0] = 1
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_rpackages$Var2 = d_its[all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_rpackages$Var2]
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_rpackages$value = plyr::revalue(as.factor(all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_rpackages$value), c("0"="Zero in all betas' CI", "1"="Zero not in all betas's CI"))


all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_both = melt(list(all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset,
          all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_rpackages), id.vars=c('value', 'Var2'))
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_both$L1 = c('TMB model', 'R packages')[all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_both$L1]
ggplot(all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset_both)+
  geom_bar(aes(group=value,
               y=factor(Var2, levels=rev(d_its)), fill=value), position='stack')+
  facet_wrap(.~L1, scale='free_x')+
  labs(y='d')+theme_bw()+theme(legend.position = "bottom", legend.title = element_blank())
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_zerobetaslope_per_dataset_both.pdf", width=6, height = 3)


## Wald test
pvals_TMB = apply(results_inference_multivariate_zerosslope_d, 2, function(j){
  sapply(j, function(i) wald_TMB_wrapper(i$tmb_object, verbatim = F))
})
pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_distrib_pvals.pdf", height=5)
hist(as.vector(pvals_TMB), breaks=60,
     xlab='P-values under null', main='Distribution of p-values in non-differentially abundant datasets')
dev.off()

pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_distrib_pvals_per_d.pdf", height=3)
par(mfrow=c(1,ncol(results_inference_multivariate_zerosslope_d)))
sapply(1:ncol(pvals_TMB), function(j) {
  hist(pvals_TMB[,j], breaks=10,
       xlab='P-values under null', main=paste0('d=', d_its[j]))
})
dev.off()

sum(pvals_TMB*prod(dim(pvals_TMB)) <= 0.05)
sum(pvals_TMB*prod(dim(pvals_TMB)) > 0.05)

pvals_TMB_adj = apply(pvals_TMB, 2, p.adjust)

ggplot(melt(pvals_TMB_adj <= 0.05), aes(x=Var2, y=value, group=Var2))+geom_jitter()+
  geom_label(data=melt(pvals_TMB_adj <= 0.05) %>% group_by(Var2) %>% summarize(mean=mean(value)), aes(x=Var2, y = TRUE, label=mean))

ggplot()+
  geom_bar(data=melt(pvals_TMB_adj <= 0.05), aes(y=d_its[Var2], fill=value),col='black', position = 'stack')+
  scale_fill_manual(values = c("#e9ffdb", "#e4717a"))+
  geom_text(data=melt(pvals_TMB_adj <= 0.05) %>% group_by(Var2) %>% summarize(mean=round(1-(mean(value, na.rm = T))), 5) %>% mutate(value=300),
             aes(y=d_its[Var2], x = value, label=mean))+
  theme_bw()+ labs(y='d')+theme(legend.position = "bottom")+guides(fill=FALSE)+
  ggtitle('Fraction of false positives (red) and true negatives (green)')
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_false_negatives_per_d_adjusted.pdf", height=3, width = 6)

ggplot()+
  geom_bar(data=melt(pvals_TMB <= 0.05), aes(y=d_its[Var2], fill=value),col='black', position = 'stack')+
  scale_fill_manual(values = c("#e9ffdb", "#e4717a"))+
  geom_text(data=melt(pvals_TMB <= 0.05) %>% group_by(Var2) %>% summarize(mean=round(1-(mean(value, na.rm = T)), 5)) %>% mutate(value=300),
            aes(y=d_its[Var2], x = value, label=mean))+
  theme_bw()+ labs(y='d')+theme(legend.position = "bottom")+guides(fill=FALSE)+
  ggtitle('Fraction of false positives (red) and true negatives (green)')
ggsave("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_false_negatives_per_d.pdf", height=3, width = 6)

#----------------------------------------------------------------------------------#
## the fraction of betas the CI of which doesn't include zero is
all_betas_zero_slope_severald_melt_fractionCI %>% group_by(Var2) %>% summarize(mean=mean(value, na.rm = T)) %>% select(mean) %>% t
## the fraction of datasets with some beta the CI of which doesn't include zero is
all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset %>% group_by(Var2) %>% summarize(mean=table(value)[1]/length(value)) %>% select(mean) %>% t
## the fraction of differentially abundant datasets
melt(pvals_TMB <= 0.05) %>% group_by(Var2) %>% summarize(mean=mean(value, na.rm = T))%>% select(mean) %>% t
## put everything in a table
gsub(x = xtable::xtable(rbind.data.frame(fraction_CI_with_zero = all_betas_zero_slope_severald_melt_fractionCI %>% group_by(Var2) %>% summarize(mean=mean(value, na.rm = T)) %>% select(mean) %>% t,
                 fraction_datasets_all_CI_with_zero = all_betas_zero_slope_severald_melt_fractionCI_all_zero_in_dataset %>% group_by(Var2) %>% summarize(mean=table(value)[1]/length(value)) %>% select(mean) %>% t,
                 fraction_differentially_abundant = melt(pvals_TMB <= 0.05) %>% group_by(Var2) %>% summarize(mean=mean(value, na.rm = T))%>% select(mean) %>% t)),
     pattern = "\\_", replacement = " ")
#----------------------------------------------------------------------------------#

## Vary a single coefficient at a time

single_coef_vec = seq(0.1, 4, length.out = 200)
scaling_factors = c(1, 5, 7)
single_coef_res = lapply(scaling_factors, function(scalefact) mclapply(single_coef_vec, function(j) replicate_inference_varysinglecoef(n = 100, d = 4, nreplicas = 1,
                                   betashape_intersect = runif(3)*scalefact,
                                   betashape_slope = c(0.2,j,2))))
pdf("../../../results/TMB/recovery_zero_betas_TMB/FE_multinomial_single_coef_variance_interceptscaling.pdf",
    height=3, width = 8)
par(mfrow=c(1,length(scaling_factors)))
sapply(1:length(scaling_factors), function(i){
  plot(single_coef_vec, sapply(single_coef_res[[i]], function(i) i[[1]]$betas_estimate[2,2]),
     main=paste0('Scaling factor for intercept: ', scaling_factors[i]), xlab='True varying slope coefficient',
     ylab='Estimated varying slope coefficient', pch=6, cex=0.8)
abline(coef = c(0,1), col='blue', lwd=2)
})
dev.off()


## Confirm that the confidence intervals are correct: for the 95% one, in 95% of the runs (simulating under the same parameter),
## the true parameter should be included in the confidence interval. I am going to simulate multiple datasets created with the
## same parameters


multiple_runs = replicate(20, replicate(n = 300, replicate_inference_varysinglecoef(n = 100, d = 4, nreplicas = 1,
                                  betashape_intersect = runif(3),
                                  betashape_slope = runif(3))))
saveRDS(multiple_runs, "../../../data/robjects_cache/tmb_results_simulations/multiple_runs.RDS")
## I take the confidence interval for each of the runs, and see for how many the true value lies inside
# multiple_runs[[1]]$betas_estimate
# give_stderr(multiple_runs[[1]]$tmb_object, only_slopes = T, only_betas = T)
# multiple_runs[[1]]$true_beta[2,]
# multiple_runs[[1]]$betas_estimate[2,]

ci_bool = apply(multiple_runs, 2, function(k) mclapply(1:length(k) , function(j) give_params_in_CI(k[[j]]$betas_estimate[2,],
                            give_stderr(k[[j]]$tmb_object, only_slopes = T, only_betas = T),
                            k[[j]]$true_beta[2,])))
## do the test
tests = apply(multiple_runs, 2, function(k) mclapply(1:length(k) , function(j) wald_TMB_wrapper(k[[j]]$tmb_object, v=F) ))

## ~80%???
par(mfrow=c(1,3))
plot(sort(sapply(ci_bool, function(k) sum(sapply(k, sum) == 3)/length(k))))
abline(h=mean(sapply(ci_bool, function(k) sum(sapply(k, sum) == 3)/length(k))))

plot(sort(sapply(tests, function(k) sum(k < 0.05)/length(k))))
abline(h=mean(sapply(tests, function(k) sum(k < 0.05)/length(k))))

plot(sapply(ci_bool, function(k) sum(sapply(k, sum) == 3)/length(k)),
     sapply(tests, function(k) sum(k < 0.05)/length(k)))