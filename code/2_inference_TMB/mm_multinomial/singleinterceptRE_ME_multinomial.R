#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
source("helper_functions.R")
# set.seed(1245)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
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

rsq = function (x, y) cor(x, y) ^ 2
#-------------------------------------------------------------------------------------------------#

stop('This is not full ME! It only has one vector of random intercepts')

#-------------------------------------------------------------------------------------------------#
TMB::compile("ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("ME_multinomial"))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
d = 4 # opt$Nk ## number of signatures
n = 10 #opt$Ns ## number of samples
beta_gamma_shape = 2  ##opt$hyperparam_shape ## shape parameter for the beta
sd_RE = 0.3 ## standard deviation for random effects
lambda = c(200, 200) ## overdispersion scalars. Lower value -> higher overdispersion
Nm_lambda = 300 ## opt$Nm_lambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)

## Group effects
## covariate matrix
X_sim = matrix(NA, nrow=2, ncol=2*n)
## the samples are split into two groups
X_sim[1,] = 1
X_sim[2,] = rep(c(0,1), each=n)
beta = matrix(0, nrow=2, ncol=d-1)
beta[2,] = rgamma(n = d-1, shape = beta_gamma_shape, rate = beta_gamma_shape) ## for the coefficients

## Random effects
Z_sim0 = matrix(0, nrow=n, ncol=n)
diag(Z_sim0) = 1
Z_sim = t(rbind(Z_sim0, Z_sim0))
u = matrix(NA, nrow=n, ncol=1)
u[,1] = rnorm(n = n, mean = 0, sd = sd_RE)

## lambdas: overdispersion
lambdas = c(rep(lambda[1], n), rep(lambda[2], n))

## create alpha
theta = softmax( cbind(t(X_sim)%*%beta + t(Z_sim)%*%replicate(d-1, u, simplify = TRUE), 0) )

W = t(apply(theta, 1, rmultinom, n = 1, size = Nm_lambda))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
data <- list(Y = W, n=n*2, d=d, x=t(X_sim), z = t(Z_sim), num_individuals=n)
# dev.off(); image(t(data$Y))

parameters <- list(
  beta = array(c(rep(runif(1, min = -4, max = 4), (d-1)),
                 rep(runif(1, min = -4, max = 4), (d-1))),
               dim = c(2,d-1)),
  u_random_effects = matrix(rep(0, n)),
  logSigma_RE=0.5
)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
obj <- MakeADFun(data, parameters, DLL="ME_multinomial", random = "u_random_effects")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
# obj$he()    ## <-- Analytical hessian
sdreport(obj)
#-------------------------------------------------------------------------------------------------#
beta

#-------------------------------------------------------------------------------------------------#
## Comparing true estimate values
par(mfrow=c(1,2))
plot(beta, opt$par[grep(names(opt$par), pattern = 'beta')])
abline(a = c(0,1))
plot(u, opt$par[grep(names(opt$par), pattern = 'u')])
rsq(x = as.vector(u),  y = opt$par[grep('u', names(opt$par))])
abline(a = c(0,1))
rsq(x = as.vector(beta),  y = opt$par[grep('beta', names(opt$par))])
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
create_dataset = function(d=4, n=10, beta_gamma_shape=2, sd_RE=0.3, lambda=c(200,200), Nm_lambda=300){
  # d = 4 # opt$Nk ## number of signatures
  # n = 10 #opt$Ns ## number of samples
  # beta_gamma_shape = 2  ##opt$hyperparam_shape ## shape parameter for the beta
  # sd_RE = 0.3 ## standard deviation for random effects
  # lambda = c(200, 200) ## overdispersion scalars. Lower value -> higher overdispersion
  # Nm_lambda = 300 ## opt$Nm_lambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)
  
  ## Group effects
  ## covariate matrix
  X_sim = matrix(NA, nrow=2, ncol=2*n)
  ## the samples are split into two groups
  X_sim[1,] = 1
  X_sim[2,] = rep(c(0,1), each=n)
  X_sim = t(as.matrix(X_sim[2,]))
  beta = matrix(0, nrow=2, ncol=d-1)
  beta[2,] = rgamma(n = d-1, shape = beta_gamma_shape, rate = beta_gamma_shape) ## for the coefficients
  beta = t(as.matrix(beta[2,]))
  
  ## Random effects
  Z_sim0 = matrix(0, nrow=n, ncol=n)
  diag(Z_sim0) = 1
  Z_sim = t(rbind(Z_sim0, Z_sim0))
  u = matrix(NA, nrow=n, ncol=1)
  u[,1] = rnorm(n = n, mean = 0, sd = sd_RE)
  # u[,1] = rep(0, n)
  
  ## lambdas: overdispersion
  lambdas = c(rep(lambda[1], n), rep(lambda[2], n))
  
  ## create alpha
  theta = softmax( cbind(t(X_sim)%*%beta + t(Z_sim)%*%replicate(d-1, u, simplify = TRUE), 0) )
  
  W = t(apply(theta, 1, rmultinom, n = 1, size = Nm_lambda))

  data <- list(Y = W, n=n*2, d=d, x=t(X_sim), z = t(Z_sim), num_individuals=n)
  return(list(data=data, true_beta=beta, true_u=u))
}

inference_random_start = function(d, n, data, true_u, true_beta){
  parameters <- list(
    beta = array(runif(d-1, min = -4, max = 4),
                 dim = c(1,d-1)),
    u_random_effects = matrix(runif(n, min = -1, max = 1)),
    logSigma_RE=1,
    log_lambda = 2
  )

  obj <- MakeADFun(data, parameters, DLL="ME_multinomial", random = "u_random_effects")
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  # obj$he()    ## <-- Analytical hessian
  sdreport(obj) ### Hessian for fixed effects was not positive definite! is it because I am including an intercept in the fixed effects?

  return(list(x=t(data$x),
              betas_estimate = opt$par[grep(names(opt$par), pattern = 'beta')],
              u_estimate = opt$par[grep(names(opt$par), pattern = 'u')],
              true_beta = true_beta,
              true_u=as.vector(true_u),
              initial_beta = parameters$beta,
              initial_u = as.vector(parameters$u)))
}

replicate_inference = function(d=4, n=10, beta_gamma_shape=2, sd_RE=0.3, lambda=c(200,200), Nm_lambda=300, nreplicas=30){
  
  data_list = create_dataset(d, n, beta_gamma_shape, sd_RE, lambda, Nm_lambda)
  
  inference = c(replicate(nreplicas, inference_random_start(d = d, n = n, data = data_list[['data']],
                                                            data_list[['true_u']], data_list[['true_beta']]), simplify = FALSE))
  return(list(data_list, inference))
}
#-------------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------------#
create_plot = function(d=4, n=10, beta_gamma_shape=2, sd_RE=0.3, lambda=c(200,200), Nm_lambda=300){
  
  data_list = create_dataset(d, n, beta_gamma_shape, sd_RE, lambda, Nm_lambda)
  
  # dev.off(); image(t(data$Y))
  
  nreplicas = 30
  inference = c(replicate(nreplicas, inference_random_start(d = d, n = n, data = data_list[['data']],
                                                            data_list[['true_u']], data_list[['true_beta']]), simplify = FALSE))
  
  ## we are only plotting the first two betas
  lims = lapply(1:(d-1), function(beta_idx) sapply(list(c(sapply(inference, function(i) i$betas_estimate)[beta_idx,],
                                                          sapply(inference, function(i) i$initial_beta)[beta_idx,],
                                                          sapply(inference, function(i) i$true_beta)[beta_idx,])), function(j) c(min(j), max(j))))
  plot(t(sapply(inference, function(i) i$betas_estimate)[1:2,]), col=scales::alpha(factor(rep(c(1,2), each=nreplicas)), 0.2),
       pch='E', xlab='Beta 1', ylab='Beta 2', xlim=lims[[1]], ylim=lims[[2]], main='Beta (1 and 2)')
  points(t(sapply(inference, function(i) i$initial_beta)[1:2,]), col=scales::alpha(factor(rep(c(1,2), each=nreplicas)), 0.2), pch=1, cex=0.5)
  points(t(sapply(inference, function(i) i$true_beta)[1:2,]), col='blue', pch='T')
  
  lims_u = lapply(1:(d-1), function(idx) sapply(list(c(sapply(inference, function(i) i$u_estimate)[idx,],
                                                       sapply(inference, function(i) i$initial_u)[idx,],
                                                       sapply(inference, function(i) i$true_u)[idx,])), function(j) c(min(j), max(j))))
  plot(t(sapply(inference, function(i) i$u_estimate)[1:2,]), col=scales::alpha(factor(rep(c(1,2), each=nreplicas)), 0.2),
       pch='E', xlab='u 1', ylab='u 2', xlim=lims_u[[1]], ylim=lims_u[[2]], main = 'u (1 and 2)')
  points(t(sapply(inference, function(i) i$initial_u)[1:2,]), col=scales::alpha(factor(rep(c(1,2), each=nreplicas)), 0.2), pch=1, cex=0.5)
  points(t(sapply(inference, function(i) i$true_u)[1:2,]), col='blue', pch='T')
  
  # pairs(cbind(t(sapply(inference, function(i) i$u_estimate)),
  #       t(sapply(inference, function(i) i$betas_estimate))), pch=9)
  return(inference)
  
}
# data_list = create_dataset(d, n, beta_gamma_shape, sd_RE, lambda, Nm_lambda)
# 
# # dev.off(); image(t(data$Y))
# 
# nreplicas = 30
# inference = c(replicate(nreplicas, inference_random_start(d = d, n = n, data = data_list[['data']],
#                                                           data_list[['true_u']], data_list[['true_beta']]), simplify = FALSE))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
nreplicas_per_dataset = 2

grid_params = expand.grid(list(d=5,#d=4:10,
                               beta_gamma_shape = c(1, 2),
                               sd_RE=0.1, #sd_RE = c(0.1, 0.3, 1),
                               n=30, #n = c(30, 80),
                               lambda = c(as.character(paste0(c(200,200), collapse = ",")),
                                          as.character(paste0(c(100,200), collapse = ",")),
                                          as.character(paste0(c(10,10), collapse = ","))),
                               Nm_lambda = 300, nreplicas = nreplicas_per_dataset), stringsAsFactors = FALSE)

results_inference_more = apply(grid_params, 1, function(i) try(replicate_inference(d = as.numeric(i["d"]),
                                                                                    n = as.numeric(i["n"]),
                                                                                    beta_gamma_shape = as.numeric(i["beta_gamma_shape"]),
                                                                                    sd_RE = as.numeric(i["sd_RE"]),
                                                          lambda = as.numeric(strsplit(as.character(i["lambda"]), ",")[[1]]),
                                                          Nm_lambda = as.numeric(i["Nm_lambda"]), nreplicas = as.numeric(i["nreplicas"]))))

# results_inference_more = replicate(40, try(replicate_inference(d = 4, n = 10, beta_gamma_shape = 2, sd_RE = 0.3,
#                                                                lambda = c(200, 200), Nm_lambda = 300, nreplicas = nreplicas_per_dataset)))

# results_inference_more_dataset = results_inference_more[1,]
# results_inference_more_inference = t(results_inference_more[2,])
results_inference_more_dataset = lapply(results_inference_more, function(i) i[[1]])
results_inference_more_inference = lapply(results_inference_more, function(i) i[[2]])
results_inference_more_inference = results_inference_more_inference[!sapply(results_inference_more_inference, typeof) == "character"] ## errors
results_inference_more_inference = t(do.call('rbind', results_inference_more_inference))


## we are looking at the beta for the slope only, not the intercept
## as many rows as nreplicas*num_datasets, and as many columns as (d-1)*nreplicates
## but what about the true beta and the betas estimate?
## this can only be done if we have the same number of features in all runs
all_beta_scatter0 = lapply(1:ncol(results_inference_more_inference),
       function(dataset_idx) c(results_inference_more_inference[,dataset_idx][[1]]$true_beta,
                               results_inference_more_inference[,dataset_idx][[1]]$betas_estimate))
all_beta_scatter = do.call('rbind', all_beta_scatter0)
## alternative
all_beta_scatter = do.call('rbind', lapply(all_beta_scatter0, function(i) cbind(i[1:(length(i)/2)], i[(length(i)/2 + 1):(length(i))]) ))
plot_pairs_with_identity(all_beta_scatter)


rsq_df = sapply(all_beta_scatter0, function(i) rsq(i[1:(length(i)/2)], i[(length(i)/2 + 1):(length(i))]))
rsq_df = cbind.data.frame(rsq=rsq_df,
               d=sapply(results_inference_more_dataset, function(i) i$data$d),
               n=sapply(results_inference_more_dataset, function(i) i$data$n))

## comparison of recovered values for each of the beta slopes. replicates are pooled
par(mfrow=c(1,ncol(all_beta_scatter)/nreplicas_per_dataset))
sapply(lapply(1:(ncol(all_beta_scatter)/nreplicas_per_dataset), function(idx_reps){
  all_beta_scatter[,idx_reps - 1 + seq(by = ncol(all_beta_scatter)/nreplicas_per_dataset, from = 1, to = ncol(all_beta_scatter))]
}), function(i) {plot(i); abline(coef = c(0,1))})

dim(results_inference_more_inference)

## replicates are in consecutive rows, e.g. the first (d-1) columns are shared, and the last (d-1) columns are equivalent estimates
all_beta_scatter[1:2,]

## for which runs is there best recovery?
dmin1 = ncol(all_beta_scatter)/2

dmin1_datasets = apply(results_inference_more_inference, 2, function(i) length(i[[1]][["betas_estimate"]]))

par(mfrow=c(1,1), mar=c(3,3,3,3))
plot(sapply(1:nrow(all_beta_scatter), function(idx_row_sim){
  rsq(all_beta_scatter[idx_row_sim,1:dmin1], all_beta_scatter[idx_row_sim,(1+dmin1):(2*dmin1)])
}), sapply(results_inference_more_dataset, function(i) i$data$n))

plot(sapply(1:nrow(all_beta_scatter), function(idx_row_sim){
  rsq(all_beta_scatter[idx_row_sim,1:dmin1], all_beta_scatter[idx_row_sim,(1+dmin1):(2*dmin1)])
}), sapply(results_inference_more_dataset, function(i) i$data$d))


recovery_df2 = melt(lapply(1:ncol(results_inference_more_inference),
       function(dataset_idx) c(results_inference_more_inference[,dataset_idx][[1]]$true_beta,
                               results_inference_more_inference[,dataset_idx][[1]]$betas_estimate)))
recovery_df2 = cbind(recovery_df2, rep( sapply(results_inference_more_dataset, function(i) i$data$d),  dmin1_datasets*2))


recovery_df = cbind.data.frame(rsq=sapply(1:nrow(all_beta_scatter), function(idx_row_sim){
  rsq(all_beta_scatter[idx_row_sim,1:dmin1], all_beta_scatter[idx_row_sim,(1+dmin1):(2*dmin1)])
}),
d=sapply(results_inference_more_dataset, function(i) i$data$d),
n = sapply(results_inference_more_dataset, function(i) i$data$n))

ggplot(recovery_df, aes(x=d, group=d, y=rsq, col=n))+
  geom_boxplot()+geom_jitter()
ggsave("../../../results/assessing_models/TMB_multinomial_RSS_d.pdf", width = 4, height = 4)

ggplot(rsq_df, aes(x=d, group=d, y=rsq, col=n))+
  geom_boxplot()+geom_jitter()
ggsave("../../../results/assessing_models/TMB_multinomial_RSS_d_2.pdf", width = 4, height = 4)

ggplot(rsq_df, aes(x=n, group=n, y=rsq, col=n))+
  geom_boxplot()+geom_jitter()
ggsave("../../../results/assessing_models/TMB_multinomial_RSS_n_2.pdf", width = 4, height = 4)


plot(sapply(1:nrow(all_beta_scatter), function(idx_row_sim){
  rsq(all_beta_scatter[idx_row_sim,1:dmin1], all_beta_scatter[idx_row_sim,(1+dmin1):(2*dmin1)])
}), sapply(results_inference_more_dataset, function(i) mean(abs(i$true_beta))))

par(mfrow=c(1,1))
keep = ((all_beta_scatter[,1] - all_beta_scatter[,2])**2) < 100
all_beta_scatter_no_outliers = all_beta_scatter[keep,]
plot(all_beta_scatter_no_outliers, col=rep( sapply(results_inference_more_dataset, function(i) i$data$d),  dmin1_datasets*2)[keep]
)

ggplot(cbind.data.frame(beta=all_beta_scatter_no_outliers, d=rep( sapply(results_inference_more_dataset, function(i) i$data$d),  dmin1_datasets*2)[keep]),
       aes(x=beta.1, y=beta.2, col=d))+geom_point(size=.6)+theme_bw()+facet_wrap(.~d)
ggsave("../../../results/assessing_models/TMB_multinomial_RSS_d_3.pdf", width = 4, height = 4)
