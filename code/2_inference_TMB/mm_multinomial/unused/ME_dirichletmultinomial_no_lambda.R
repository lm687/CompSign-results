#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
library(scales)
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
compile("ME_dirichletmultinomial_no_lambda.cpp", "-std=gnu++17")
dyn.load(dynlib("ME_dirichletmultinomial_no_lambda"))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
# create_plot = function(){
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
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
data <- list(Y = W, n=n*2, d=d, x=t(X_sim), z = t(Z_sim), num_individuals=n)
# dev.off(); image(t(data$Y))

# inference_random_start = function(){
parameters <- list(
  beta = array(#c(rep(runif(1, min = -4, max = 4), (d-1)),
    runif(d-1, min = -4, max = 4),
    # rep(runif(1, min = -4, max = 4), (d-1)),
    dim = c(1,d-1)),
  # u = matrix(rep(0, n)),
  u = matrix(runif(n, min = -1, max = 1)),
  logSigma_RE=0
)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
obj <- MakeADFun(data, parameters, DLL="ME_dirichletmultinomial_no_lambda")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj) ### Hessian for fixed effects was not positive definite and there is no estimate for variance of variance of RE
