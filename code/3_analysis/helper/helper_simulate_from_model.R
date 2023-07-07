simulate_from_DM = function(beta_coefs, RE_coefs, lambda, coefficient_overdispersion){
  alphabar = softmax(cbind(sapply(1:(length(beta_coefs)/2),
                                  function(some_dummy_idx){
                                    give_z_matrix(length(RE_coefs) * 2) %*% RE_coefs}) +
                             give_x_matrix(length(RE_coefs) * 2) %*% matrix(beta_coefs, nrow=2) , 0))
  ## here there is a single lambda
  alpha_mat = alphabar*rep(exp(lambda), each=n)*coefficient_overdispersion
  return(t(apply(alpha_mat, 1, MCMCpack::rdirichlet, n=1)))
}

simulate_from_M_RE = function(beta_coefs, RE_coefs){
  RE_coefs = matrix(RE_coefs, ncol=length(beta_coefs)/2)
  n = length(RE_coefs)/(length(beta_coefs)/ 2)
  theta = softmax(cbind( (give_z_matrix(2*n) %*% RE_coefs) +
                           (give_x_matrix(2*n)) %*% matrix(beta_coefs, nrow=2) , 0))
  ## here there is a single lambda
  return(theta)
}

simulate_from_DM_RE = function(beta_coefs, RE_coefs, lambda, coefficient_overdispersion){
  RE_coefs = matrix(RE_coefs, ncol=length(beta_coefs)/2)
  n = length(RE_coefs)/(length(beta_coefs)/ 2)
  alphabar = softmax(cbind( (give_z_matrix(2*n) %*% RE_coefs) +
                              (give_x_matrix(2*n)) %*% matrix(beta_coefs, nrow=2) , 0))
  ## here there is a single lambda
  alpha_mat = alphabar*rep(exp(lambda), each=n)*coefficient_overdispersion
  return(t(apply(alpha_mat, 1, MCMCpack::rdirichlet, n=1)))
}

simulate_from_DM_RE_altpar = function(beta_coefs, RE_coefs, lambda){
  RE_coefs = matrix(RE_coefs, ncol=length(beta_coefs)/2)
  n = length(RE_coefs)/(length(beta_coefs)/ 2)
  alphabar = softmax(cbind( (give_z_matrix(2*n) %*% RE_coefs) +
                              (give_x_matrix(2*n)) %*% matrix(beta_coefs, nrow=2) , 0))
  ## here there is a single lambda
  alpha_mat = alphabar*rep(1/exp(lambda), each=n)
  return(t(apply(alpha_mat, 1, MCMCpack::rdirichlet, n=1)))
}

## Simulate data with confidence intervals, to choose between M ME and DM ME

## if a CT has not run well, return NA


## DM applies to any RE structure: the important thing is that the RE coefficients have the same
## dimensionality regardless, as they are realisations
give_true_in_conf_int = function(idx_dataset, model, dataset_TMB){
  # sim_confint = sapply(1:20, function(idx_dataset){
  # sapply(6, function(idx_dataset){
  ## idx_dataset=9 is an example of a successful DM run
  # print(c(is.null(dataset_TMB[[idx_dataset]]),
  #         typeof(dataset_TMB[[idx_dataset]]) == "character",
  #         !(dataset_TMB[[idx_dataset]]$pdHess)))
  if(is.null(dataset_TMB[[idx_dataset]])){
    return(NA)
  }else{
    if(typeof(dataset_TMB[[idx_dataset]]) %in% c("character", "logical")){
      return(NA)
    }else{
      if(!(dataset_TMB[[idx_dataset]]$pdHess)){
        # no good convergence
        return(NA)
      }else{
        beta_coefs = python_like_select_name(dataset_TMB[[idx_dataset]]$par.fixed, 'beta')
        RE_coefs = dataset_TMB[[idx_dataset]]$par.random
        lambda = python_like_select_name(dataset_TMB[[idx_dataset]]$par.fixed, 'log_lambda')
        
        if(model == 'DM'){
          sim_thetas = replicate(1e3, simulate_from_DM_RE(beta_coefs, RE_coefs, lambda, coefficient_overdispersion))
        }else if(model == 'DM_altpar'){
          sim_thetas = replicate(1e3, simulate_from_DM_RE_altpar(beta_coefs, RE_coefs, lambda))
        }else if(model == 'M'){
          thetas_M = simulate_from_M_RE(beta_coefs, RE_coefs)
          sim_thetas = replicate(1e3, thetas_M)
        }else{
          stop('Indicate a correct <model>')
        }
        
        ## from observed
        matrices = slot(count_objects[[idx_dataset]], 'count_matrices_active')
        if(sum(sapply(matrices, length)) == 0){
          matrices = slot(count_objects[[idx_dataset]], 'count_matrices_all')
        }
        
        matrices = do.call('rbind', lapply(matrices, round))
        vec_observed = as.vector(matrices)
        
        mut_toll = rowSums(matrices)
        ## simulate with multinomial from these thetas
        multinom_draws0 = lapply(1:dim(sim_thetas)[3], function(simulation_idx){
          sapply(1:nrow(sim_thetas[,,simulation_idx]), function(i){
            rmultinom( n = 1, size = mut_toll[i], 
                       prob = sim_thetas[i,,simulation_idx])
          })
        })
        multinom_draws = sapply(multinom_draws0, function(i) as.vector(t(i)))
        
        ## check if real values fall in the confidence interval, for each of the thetas
        
        multinom_draws[[1]]
        multinom_draws[[2]]
        in_conf_int = sapply(1:dim(multinom_draws)[1], function(i){
          confint_bounds = quantile(multinom_draws[i,], probs = c(0.025, 0.975))
          (vec_observed[i] >= confint_bounds[1]) & (vec_observed[i] <= confint_bounds[2])
        })
        
        ## plotting with subsample
        subset_pars = 1:length(vec_observed) #sample(1:length(vec_observed), size = 10)
        ccccc=melt(multinom_draws[subset_pars,sample(1:1000, size = 100)])
        # plot(ccccc$Var1, ccccc$value)
        # points(1:length(subset_pars), vec_observed[subset_pars], col='red')
        
        ## the ones that worked
        ccccc=melt(multinom_draws[subset_pars,sample(1:1000, size = 100)][in_conf_int,])
        # plot(ccccc$Var1, ccccc$value)
        # points(1:sum(in_conf_int), vec_observed[subset_pars][in_conf_int], col='red')
        
        ## the ones that didn't work
        ccccc=melt(multinom_draws[subset_pars,sample(1:1000, size = 100)][!in_conf_int,])
        # plot(ccccc$Var1, ccccc$value)
        # points(1:sum(!in_conf_int), vec_observed[subset_pars][!in_conf_int], col='red')
        
        return(in_conf_int)
        
        # ml_thetas = normalise_rw(do.call('rbind', matrices))
        # ml_thetas = replicate(20, ml_thetas)
        # dim(melt(sim_thetas))
        # dim(melt(ml_thetas))
        # plot(cbind(sim=melt(sim_thetas), ml=melt(ml_thetas))[,c('sim.value', 'ml.value')],
        #      main=names(count_objects[idx_dataset]), cex.main=.7)
        # abline(coef = c(0,1), col='blue', lty='dashed')
      }
    }
  }
}
