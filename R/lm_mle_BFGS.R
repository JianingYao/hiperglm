#' calculate log-likelihood under a linear model
log_likelihood <- function(design, outcome, beta, noise_var = 1) {
  est <- design %*% beta
  loglik <- -0.5 / noise_var * sum((outcome - est)^2)
  return(loglik)
}

#' calculate gradient under a linear model
loglik_gradient <- function(design, outcome, beta, noise_var = 1) {
  est <- design %*% beta
  grad <- t(design) %*% (outcome - est) / noise_var
  return(grad)
}

#' implement BFGS to find MLE
lm_mle_BFGS <- function(design, outcome, noise_var = 1) {
  np <- ncol(design)
  result <- stats::optim(
    par = rep(0, np), fn = log_likelihood, gr = loglik_gradient,
    design = design, outcome = outcome, noise_var = noise_var,
    method = "BFGS", control = list(fnscale = -1)
  )
  beta_est <- result$par
  return(beta_est)
}
