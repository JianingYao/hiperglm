#' check if two inputs are close enough
are_all_close <- function(v, w, abs_tol = 1e-6, rel_tol = 1e-6) {
  abs_diff <- abs(v - w)
  are_all_within_atol <- all(abs_diff < abs_tol)
  are_all_within_rtol <- all(abs_diff < rel_tol * pmax(abs(v), abs(w)))
  return(are_all_within_atol && are_all_within_rtol)
}

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

#' calculate sigmoid function for a logistic regression model
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' calculate log-likelihood under a logistic regression model
logit_log_likelihood <- function(design, outcome, beta) {
  est <- design %*% beta
  loglik <- sum(outcome * est - log(1 + exp(est)))
  return(loglik)
}

#' calculate gradient under a logistic regression model
logit_loglik_gradient <- function(design, outcome, beta) {
  est <- design %*% beta
  mu <- sigmoid(est)
  grad <- as.vector(t(design) %*% (outcome - mu))
  return(grad)
}

#' calculate weights based on design matrix and outcome
logit_weights <- function(design, beta) {
  est <- design %*% beta
  mu <- as.vector(sigmoid(est))
  weight <- mu * (1 - mu)
  return(weight)
}

#' calculate Hessian of the log-likelihood under a logistic regression model
logit_hessian <- function(design, beta) {
  weight <- diag(logit_weights(design, beta))
  hess <- -t(design) %*% weight %*% design
  return(hess)
}
