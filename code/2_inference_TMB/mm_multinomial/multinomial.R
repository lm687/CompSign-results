rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
source("../helper_TMB.R")
compile("multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("multinomial"))
# set.seed(123)
n = 50 #100
d = 4
true_theta = MCMCpack::rdirichlet(1, rep(1, d))
data <- list(Y = t(rmultinom(n, 300, true_theta)), n=n, d=d)
parameters <- list(theta_prime=rep(1, d-1))
obj <- MakeADFun(data, parameters, DLL="multinomial")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
# obj$he()    ## <-- Analytical hessian
sdreport(obj)

closure = function(v) v/sum(v)

theta_est = closure(exp(c( opt$par, 0)))

plot(true_theta, theta_est)

true_theta


data2 <- list(Y = rbind(data$Y, c(20, 40, 10, 27)), n=n+1, d=d)
data3 <- list(Y = rbind(data$Y, c(20, 40, 10, 27)*4), n=n+1, d=d)
give_inf_res = function(dta){
  obj <- MakeADFun(dta, parameters, DLL="multinomial")
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  # obj$he()    ## <-- Analytical hessian
  sdreport(obj)
}
inf2 <- give_inf_res(data2)
inf3 <- give_inf_res(data3)

sdreport(obj)$par.fixed
inf2$par.fixed ## the addition of a sample with smaller mutational burden changes more the estimates
inf3$par.fixed

## data
# Estimate Std. Error
# theta_prime -0.6814748  0.1871377
# theta_prime -0.3155504  0.1670145
# theta_prime  0.2578137  0.1444141
softmax(c(-0.6814748, -0.3155504, 0.2578137, 0))

# data2
# Estimate Std. Error
# theta_prime 20.335478   8237.675
# theta_prime -2.963194  37169.382
# theta_prime -2.963194  37169.382
# Maximum gradient component: 1.625963e-08 
softmax(c(20.335478, -2.963194, -2.963194, 0))

# data3
# Estimate Std. Error
# theta_prime  29.12309   334054.8
# theta_prime -10.87691 76711528.4
# theta_prime -10.87691 76711528.4
softmax(c(29.12309, -10.87691, -10.87691, 0))

