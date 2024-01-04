b <- readRDS("/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_betaintercept2.RDS")
length(b)

a <- cbind(runif(100), runif(100), runif(100), runif(100))
a <- cbind(a, a[,1]+runif(100, min = -0.3, max = -0.3))

cv <- cov(a)
saveRDS(cv, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_covmat1.RDS")

c <- cbind(runif(100, min = -7, max = 7), runif(100), runif(100, min = -3, max = 3),  runif(100, min = -7, max = 7))
c <- cbind(c[,1]+runif(100, min = -0.3, max = -0.3), c)
cv2 <- cov(c)
cv2

saveRDS(cv2, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_covmat2.RDS")

beta_6 <- runif(6, min = -4, max = 6)
saveRDS(beta_6, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_betaintercept1d7.RDS")
beta_6 <- runif(6, min = -3, max = 2)
saveRDS(beta_6, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_betaslope1d7.RDS")

readRDS("/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_betaintercept1d7.RDS")
readRDS("/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_betaslope1d7.RDS")

c <- cbind(runif(100, min = -3, max = 1), runif(100), runif(100, min = -3, max = 3),  runif(100, min = -7, max = 7))
c <- cbind(c[,1]+runif(100, min = -1.3, max = 2.3), c)
c <- cbind(-c[,2]+runif(100, min = -2.3, max = -0.3), c)
cv2 <- cov(c)
cv2
saveRDS(cv2, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_covmat1d7.RDS")

cv3 <- diag(exp(runif(n = 6)))
saveRDS(cv3, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_covmat2d7.RDS")

##----

beta_4 <- runif(4, min = -4, max = 6)
saveRDS(beta_4, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_betaintercept1d4.RDS")
beta_4 <- runif(4, min = -3, max = 2)
saveRDS(beta_4, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_betaslope1d4.RDS")

readRDS("/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_betaintercept1d4.RDS")
readRDS("/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_betaslope1d4.RDS")

c <- cbind(runif(100, min = -3, max = 1), runif(100))
c <- cbind(c[,1]+runif(100, min = -1.3, max = 2.3), c)
c <- cbind(-c[,2]+runif(100, min = -2.3, max = -0.3), c)
cv2 <- cov(c)
cv2
saveRDS(cv2, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_covmat1d4.RDS")

cv3 <- diag(exp(runif(n = 6)))
saveRDS(cv3, "/Users/morril01/Documents/PhD/GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_covmat2d7.RDS")
