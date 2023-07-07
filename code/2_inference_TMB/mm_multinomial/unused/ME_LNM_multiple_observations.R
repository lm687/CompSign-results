## Same as ME_LNM.R but with multiple observations for each individual and group
## Still
# logSigma_RE -1.044348e+01          NaN
# Warning:
#   Hessian of fixed effects was not positive definite.
# Maximum gradient component: 29.91629 
# Warning message:
#   In sqrt(diag(object$cov.fixed)) : NaNs produced
#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
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

#-------------------------------------------------------------------------------------------------#
compile("ME_LNM.cpp", "-std=gnu++17")
dyn.load(dynlib("ME_LNM"))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
d = 2 # opt$Nk ## number of signatures
nobs_per_person_per_group = 2
n = 20 #opt$Ns ## number of samples
beta_gamma_shape = 2  ##opt$hyperparam_shape ## shape parameter for the beta
sd_RE = 2 #0.3 ## standard deviation for random effects
lambda = c(200, 200) ## overdispersion scalars. Lower value -> higher overdispersion
Nm_lambda = 300 ## opt$Nm_lambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)

## Group effects
## covariate matrix
X_sim = matrix(NA, nrow=2, ncol=n*2*nobs_per_person_per_group)
## the samples are split into two groups
X_sim[1,] = 1
X_sim[2,] = rep(c(0,1), each=n*nobs_per_person_per_group)
X_sim = t(as.matrix(X_sim[2,])) ## remove intercept
beta = matrix(0, nrow=2, ncol=d-1)
beta[2,] = rgamma(n = d-1, shape = beta_gamma_shape, rate = beta_gamma_shape) ## for the coefficients
beta = t(as.matrix(beta[2,]))

## Random effects
Z_sim0 = matrix(0, nrow=n, ncol=n)
diag(Z_sim0) = 1
Z_sim = t(do.call('rbind', replicate(2*nobs_per_person_per_group, Z_sim0, simplify = FALSE)))
Z_sim = Z_sim/4 ## trying to scale the covarariates
u = matrix(NA, nrow=n, ncol=1)
u[,1] = rnorm(n = n, mean = 0, sd = sd_RE)

## lambdas: overdispersion
lambdas = c(rep(lambda[1], n), rep(lambda[2], n))

## create alpha
theta = softmax( cbind(t(X_sim)%*%beta + t(Z_sim)%*%replicate(d-1, u, simplify = TRUE), 0) )

W = t(apply(theta, 1, rmultinom, n = 1, size = Nm_lambda))
#-------------------------------------------------------------------------------------------------#
image(t(W))

#-------------------------------------------------------------------------------------------------#
data <- list(Y = W, x=t(X_sim), z = t(Z_sim))
# dev.off(); image(t(data$Y))

parameters <- list(
  # beta = matrix(c(rep(runif(1, min = -4, max = 4), (d-1)),
  #   rep(runif(1, min = -4, max = 4), (d-1))),
  # nrow = 2, byrow=TRUE),
  beta = matrix(rep(runif(1, min = -4, max = 4), (d-1)),
                nrow = 1, byrow=TRUE),
  u = matrix(rep(1, n)),
  logSigma_RE=1
  # mu_par = rep(0, (d-1)*(d-1)-(d-1))/2
  # cov_par = rep(0, (d-1)*(d-1)-(d-1))/2
)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
obj <- MakeADFun(data, parameters, DLL="ME_LNM")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)
#-------------------------------------------------------------------------------------------------#

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

