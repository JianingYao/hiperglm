#' calculate log-likelihood under a linear model
lm_log_likelihood <- function(design, outcome, beta, noise_var = 1) {
  est <- design %*% beta
  loglik <- -0.5 / noise_var * sum((outcome - est)^2)
  return(loglik)
}

#' calculate gradient under a linear model
lm_loglik_gradient <- function(design, outcome, beta, noise_var = 1) {
  est <- design %*% beta
  grad <- t(design) %*% (outcome - est) / noise_var
  return(grad)
}

#' implement BFGS to find MLE under a linear model
lm_mle_BFGS <- function(design, outcome, noise_var = 1) {
  init_coef <- rep(0, ncol(design))
  result <- stats::optim(
    par = init_coef, fn = lm_log_likelihood, gr = lm_loglik_gradient,
    design = design, outcome = outcome, noise_var = noise_var,
    method = "BFGS", control = list(fnscale = -1)
  )
  optim_converged <- (result$convergence == 0L)
  if (!optim_converged) {
    warning("Optimization did not converge. The estimates may be meaningless.")
  }
  beta_est <- result$par
  return(beta_est)
}
