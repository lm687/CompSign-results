#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
library(ggplot2)
source("helper_functions.R")
# set.seed(1234)
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

#-------------------------------------------------------------------------------------------------#
TMB::compile("ME_LNM.cpp", "-std=gnu++17")
dyn.load(dynlib("ME_LNM"))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
d = 2 # opt$Nk ## number of signatures
n = 20 #opt$Ns ## number of samples
beta_gamma_shape = 0.05  ##opt$hyperparam_shape ## shape parameter for the beta
sd_RE = 1.3 ## standard deviation for random effects
lambda = c(200, 200) ## overdispersion scalars. Lower value -> higher overdispersion
Nm_lambda = 300 ## opt$Nm_lambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)

## Group effects
## covariate matrix
X_sim = matrix(NA, nrow=2, ncol=2*n)
## the samples are split into two groups
X_sim[1,] = 1
X_sim[2,] = rep(c(0,1), each=n)
X_sim = t(as.matrix(X_sim[2,])) ## remove intercept
beta = matrix(0, nrow=2, ncol=d-1)
beta[2,] = rgamma(n = d-1, shape = beta_gamma_shape, rate = beta_gamma_shape) ## for the coefficients
beta = t(as.matrix(beta[2,]))

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
data <- list(Y = W, x=t(X_sim), z = t(Z_sim))
# dev.off(); image(t(data$Y))

parameters <- list(
  beta = matrix(rep(runif(1, min = -4, max = 4), (d-1)),
                nrow = 1, byrow=TRUE),
  u_random_effects = matrix(rep(1, n)),
  logSigma_RE=1,
  theta_prime = matrix(0, nrow=nrow(W), ncol=d-1)
  # u = matrix(rep(1, n)),
  # theta_prime = matrix(0, nrow=nrow(W), ncol=d-1)
)

#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
obj <- MakeADFun(data, parameters, DLL="ME_LNM", random = "u_random_effects")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
sdreport(obj)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## Comparing true estimate values
par(mfrow=c(1,2))
plot(beta, opt$par[grep(names(opt$par), pattern = 'beta')])
abline(coef = c(0,1))
# plot(u, opt$par[grep(names(opt$par), pattern = 'u')])
# rsq(x = as.vector(u),  y = opt$par[grep('u', names(opt$par))])
# abline(coef = c(0,1))
rsq(x = as.vector(beta),  y = opt$par[grep('beta', names(opt$par))])
#-------------------------------------------------------------------------------------------------#

## Note change in scale
# plot(u, opt$par[names( opt$par) == 'u'])

t(X_sim) %*% (matrix(opt$par[names(opt$par) == "beta"]))

# theta_est = softmax(cbind(t(Z_sim)%*%replicate(d-1, opt$par[names(opt$par) == "u"], simplify = TRUE), 0))

# dev.off(); par(mfrow=c(1,2))
# image(t(theta), zlim = c(0,1))
# image(t(theta_est), zlim = c(0,1))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
create_dataset = function(d = 2, n = 20, beta_gamma_shape = 0.02, sd_RE = 1.3, lambda = c(200, 200), Nm_lambda = 300){
  ## Group effects
  ## covariate matrix
  X_sim = matrix(NA, nrow=2, ncol=2*n)
  ## the samples are split into two groups
  X_sim[1,] = 1
  X_sim[2,] = rep(c(0,1), each=n)
  X_sim = t(as.matrix(X_sim[2,])) ## remove intercept
  beta = matrix(0, nrow=2, ncol=d-1)
  beta[2,] = rgamma(n = d-1, shape = beta_gamma_shape, rate = beta_gamma_shape) ## for the coefficients
  beta = t(as.matrix(beta[2,]))
  
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
  data <- list(Y = W, n=n*2, d=d, x=t(X_sim), z = t(Z_sim), num_individuals=n)
  return(list(data=data, true_beta=beta, true_u=u, true_theta=theta))
}

inference_random_start = function(d, n, data, true_u, true_beta, true_theta){
  parameters <- list(
    
    beta = array(#c(rep(runif(1, min = -4, max = 4), (d-1)),
      runif(d-1, min = -4, max = 4),
      dim = c(1,d-1)),
    u_random_effects = matrix(runif(n, min = -1, max = 1)),
    logSigma_RE=0,
    theta_prime = matrix(0, nrow=nrow(data$Y), ncol=d-1)
  )
  #-------------------------------------------------------------------------------------------------#
  
  #-------------------------------------------------------------------------------------------------#
  obj <- MakeADFun(data, parameters, DLL="ME_LNM", random = "u_random_effects")
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  # obj$he()    ## <-- Analytical hessian
  sdreport(obj) ### Hessian for fixed effects was not positive definite! is it because I am including an intercept in the fixed effects?
  #-------------------------------------------------------------------------------------------------#
  
  return(list(x=t(data$x),
              z = t(data$z),
              betas_estimate = opt$par[grep(names(opt$par), pattern = 'beta')],
              u_estimate = opt$par[grep(names(opt$par), pattern = 'u')],
              true_beta = true_beta,
              true_u=as.vector(true_u),
              initial_beta = parameters$beta,
              initial_u = as.vector(parameters$u),
              true_theta=true_theta))
}


replicate_inference = function(d=4, n=10, beta_gamma_shape=2, sd_RE=0.3, lambda=c(200,200), Nm_lambda=300, nreplicas=30){
  data_list = create_dataset(d, n, beta_gamma_shape, sd_RE, lambda, Nm_lambda)
  inference = c(replicate(nreplicas, inference_random_start(d = d, n = n, data = data_list[['data']], true_theta=data_list[['true_theta']],
                                                            data_list[['true_u']], data_list[['true_beta']]), simplify = FALSE))
  return(inference)
}

results_inference_more = replicate(60, replicate_inference(d=4, n=20, beta_gamma_shape=0.2, sd_RE=1.3,
                                                           lambda = c(200, 200), Nm_lambda = 300,
                                                           nreplicas = 2))
## I have no idea what this is
plot_pairs_with_identity(t(sapply(1:ncol(results_inference_more),
                                  function(dataset_idx) c(results_inference_more[1,dataset_idx][[1]]$true_beta,
                                                              results_inference_more[1,dataset_idx][[1]]$betas_estimate))), pch=19)

plot_pairs_with_identity(t(sapply(1:ncol(results_inference_more),
                                  function(dataset_idx) cbind(results_inference_more[1,dataset_idx][[1]]$true_u,
                                                              results_inference_more[1,dataset_idx][[1]]$u_estimate))))
all_u_scatter = do.call('rbind', lapply(1:ncol(results_inference_more),
                                        function(dataset_idx) cbind(results_inference_more[1,dataset_idx][[1]]$true_u,
                                                                    results_inference_more[1,dataset_idx][[1]]$u_estimate)))
par(mfrow=c(1,1), mar=c(2,2,2,2)); plot(all_u_scatter); abline(coef=c(0,1), col='blue')

##' there clearly are two types of runs: those which have well recovered and those which are extremely close to zero
##' Now I try to find which parameters/data give this poorly recovered second group
hist(all_u_scatter[,2], breaks = 100) ## these are the datapoints in the middle that I need to select
par(mfrow=c(1,1), mar=c(2,2,2,2)); plot(all_u_scatter, col=); abline(coef=c(0,1), col='blue')



## conclusion: very large beta does not recover well the parameters and sets random effects to zero

lapply(1:ncol(results_inference_more), function(i) results_inference_more[1,i][[1]]$u_estimate)

## good accuracy requires high sd for random effects, and low beta. otherwise the random effects go to zero

results_inference_more[1,1][[1]]$true_theta
results_inference_more[1,1][[1]]$u_estimate
results_inference_more[1,1][[1]]$true_u
results_inference_more[1,1][[1]]$x

# When are the random effects not essentially zerOo? (i.e. well recovered)
sort(sapply(1:dim(results_inference_more)[2], function(i) sum(results_inference_more[1,i][[1]]$u_estimate)**2))


# RSS
rss = outer(1:nrow(results_inference_more), 1:ncol(results_inference_more), Vectorize(function(i,j){
  sum((results_inference_more[i,j][[1]]$true_u - results_inference_more[i,j][[1]]$u_estimate)**2)
}))

## within a dataset, it doesn't change

sort(apply(rss, 2, mean))
order(apply(rss, 2, mean))

give_rss = function(i, j){sum((i-j)**2)}
give_rss(results_inference_more[1,3][[1]]$u_estimate,
         results_inference_more[1,3][[1]]$true_u)

sapply(1:ncol(results_inference_more), function(i) results_inference_more[1,i][[1]]$true_u)
sapply(1:ncol(results_inference_more), function(i) results_inference_more[1,i][[1]]$u_estimate)

## if we're analysing several of the same dataset
# (results_inference_more[1,1])
# all_u_scatter = do.call('rbind', lapply(1:nrow(results_inference_more),
#                                         function(run_idx) cbind(results_inference_more[run_idx,1][[1]]$true_u,
#                                                                     results_inference_more[run_idx,1][[1]]$u_estimate)))
# par(mfrow=c(1,1), mar=c(2,2,2,2)); plot(all_u_scatter); abline(coef=c(0,1), col='blue')
# 
# plot_pairs_with_identity(t(sapply(1:ncol(results_inference_more),
#                                   function(dataset_idx) cbind(results_inference_more[1,dataset_idx][[1]]$true_u,
#                                                               results_inference_more[1,dataset_idx][[1]]$u_estimate))))
# 
# 
# sapply(1:nrow(results_inference_more), function(i) results_inference_more[i,1][[1]]$u_estimate)
# results_inference_more[1,1][[1]]$true_u
