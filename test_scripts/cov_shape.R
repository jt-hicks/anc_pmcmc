cov_var_time <- function(mcmc_pars) {
  sapply(2:nrow(mcmc_pars), function(i){
    matrix <- cov(mcmc_pars[1:i,])
    shape <- sqrt(det(matrix))
  })
}