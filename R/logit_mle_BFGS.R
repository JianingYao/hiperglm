#' implement BFGS to find MLE under a logistic regression model
logit_mle_BFGS <- function(design, outcome) {
  init_coef <- rep(0, ncol(design))
  obj_fn <- function(coef) {
    logit_log_likelihood(design, outcome, coef)
  }
  obj_grad <- function(coef) {
    logit_loglik_gradient(design, outcome, coef)
  }
  result <- stats::optim(
    par = init_coef, fn = obj_fn, gr = obj_grad,
    method = "BFGS", control = list(fnscale = 1)
  )
  optim_converged <- (result$convergence == 0L)
  if (!optim_converged) {
    warning("Optimization did not converge. The estimates may be meaningless.")
  }
  beta_est <- result$par
  return(beta_est)
}
