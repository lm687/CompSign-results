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
compile("ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("ME_dirichletmultinomial"))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
create_plot = function(){
  d = 3 # opt$Nk ## number of signatures
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
  
  inference_random_start = function(){
    parameters <- list(
      beta = array(#c(rep(runif(1, min = -4, max = 4), (d-1)),
        runif(d-1, min = -4, max = 4),
        # rep(runif(1, min = -4, max = 4), (d-1)),
        dim = c(1,d-1)),
      # u = matrix(rep(0, n)),
      # u = matrix(runif(n, min = -1, max = 1)),
      u = matrix(rnorm(n, mean = 0, sd = rgamma(n = 1, shape = 1.2, rate = 1))),
      logSigma_RE=0,
      log_lambda = 2
    )
    #-------------------------------------------------------------------------------------------------#
    
    #-------------------------------------------------------------------------------------------------#
    obj <- MakeADFun(data, parameters, DLL="ME_dirichletmultinomial")
    obj$hessian <- TRUE
    opt <- do.call("optim", obj)
    opt
    opt$hessian ## <-- FD hessian from optim
    obj$he()    ## <-- Analytical hessian
    report = sdreport(obj) ### Hessian for fixed effects was not positive definite! is it because I am including an intercept in the fixed effects?
    #-------------------------------------------------------------------------------------------------#
    # report whether hessian is positive definite
    return(list(parameters=parameters, pdHess=report$pdHess))
  }
  
  reps_inference = replicate(n = 300, expr = inference_random_start(), simplify = FALSE)
  # sapply(reps_inference, function(i) i$pdHess) ## whether the Hessian of fixed effects was positive definite.
  
  # plot(t(sapply(reps_inference, function(i) i$parameters$beta)), pch=19,
  #      col=factor(sapply(reps_inference, function(i) i$pdHess)))
  # 
  # ## first two individuals
  # plot(t(sapply(reps_inference, function(i) i$parameters$u)[1:2,]), pch=19,
  #      col=factor(sapply(reps_inference, function(i) i$pdHess)))
  # 
  pca = prcomp(do.call('rbind', lapply(reps_inference, function(i) unlist(i$parameters))))
  plot(pca$x[,1:2], pch=19,
       col=factor(sapply(reps_inference, function(i) i$pdHess)), main='PCA')
  return(reps_inference)
}

par(mfrow=c(2,3))
## 6 columns (for 6 independent cases) and 300 (for each case, 300 sets of initial values)
results_replicates = replicate(6, create_plot())

logistic_regression = function(idx_dataset){
  df = cbind(response=sapply(results_replicates[,idx_dataset], function(i) i$pdHess),
        do.call('rbind', lapply(results_replicates[,idx_dataset], function(i) unlist(i$parameters))))
  return(glm(formula = as.formula(paste0('response~',paste0(colnames(df)[-1], collapse = '+'))), family = "binomial", 
           data = data.frame(df)))
}
    
logistic_regression_res = lapply(1:dim(results_replicates)[2], logistic_regression)
logistic_regression_res_coef = sapply(logistic_regression_res, function(i) i$coefficients)
par(mfrow=c(1,1))
plot(as.vector(apply(logistic_regression_res_coef, 1, as.vector)), type = "h", col=rep(1:dim(logistic_regression_res_coef)[1], each=dim(logistic_regression_res_coef)[2]))


