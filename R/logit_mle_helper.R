sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' calculate log-likelihood under a logistic regression model
logit_log_likelihood <- function(design, outcome, beta) {
  est <- design %*% beta
  mu <- sigmoid(est)
  loglik <- sum(outcome * log(mu) + (1 - outcome) * log(1 - mu))
  return(-loglik)
}

#' calculate gradient under a logistic regression model
logit_loglik_gradient <- function(design, outcome, beta) {
  est <- design %*% beta
  mu <- sigmoid(est)
  grad <- t(design) %*% (mu - outcome)
  return(grad)
}

#' calculate Hessian of the log-likelihood under a logistic regression model
logit_hessian <- function(design, beta) {
  est <- design %*% beta
  mu <- as.vector(sigmoid(est))
  W <- diag(mu * (1 - mu))
  Hess <- t(design) %*% W %*% design
  return(Hess)
}
