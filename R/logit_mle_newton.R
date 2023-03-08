#' implement Newton's method to find MLE under a logistic regression model
logit_mle_newton <- function(design, outcome, tol = NULL, max_iter = 1000) {
  n_pred <- ncol(design)
  if (is.null(tol)) {
    threshold <- 1e-5
    tol <- stats::qchisq(threshold, df = n_pred)
  }

  beta <- rep(0, n_pred)
  log_lik_old <- logit_log_likelihood(design, outcome, beta)
  diff <- Inf
  rel_diff <- Inf
  iter <- 1
  while ((iter <= max_iter) && ((diff > tol) || (rel_diff > tol))) {
    grad <- logit_loglik_gradient(design, outcome, beta)
    Hess <- logit_hessian(design, beta)
    beta <- beta - solve(Hess, grad)
    log_lik_new <- logit_log_likelihood(design, outcome, beta)
    diff <- abs(log_lik_new - log_lik_old)
    rel_diff <- diff / abs(log_lik_old)

    if ((diff <= tol) && (rel_diff <= tol)) {
      cat("MLE is found after ", iter, " iterations with convergence.")
      return(beta)
    } else if (iter == max_iter) {
      cat("MLE finder failed to converge after ", iter, " iterations. The estimates may be meaningless.")
      return(beta)
    } else {
      iter <- iter + 1
      log_lik_old <- log_lik_new
    }
  }
}
