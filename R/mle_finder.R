#' use Cholesky decomposition to perform pseudo-inverse for a linear model
mle_pinv <- function(design, outcome) {
  a <- crossprod(design)
  b <- crossprod(design, outcome)
  upper <- chol(a)
  z <- backsolve(upper, b, transpose = TRUE)
  beta <- as.vector(backsolve(upper, z))
  return(list(coef = beta))
}

#' implement Newton's method to find MLE under a logistic regression model
mle_newton <- function(design, outcome, tol = NULL, max_iter = 50) {
  n_pred <- ncol(design)
  if (is.null(tol)) {
    tol <- 1e-5
  }
  beta <- rep(0, n_pred)
  log_lik_old <- logit_log_likelihood(design, outcome, beta)
  iter <- 0
  converged <- FALSE
  while ((iter < max_iter) && (!converged)) {
    grad <- logit_loglik_gradient(design, outcome, beta)
    hess <- logit_hessian(design, beta)
    beta <- beta - solve(hess, grad)
    log_lik_new <- logit_log_likelihood(design, outcome, beta)
    converged <- are_all_close(log_lik_old, log_lik_new,
      abs_tol = tol, rel_tol = tol
    )
    iter <- iter + 1
    log_lik_old <- log_lik_new
  }
  if (converged) {
    cat("MLE is found after ", iter, " iterations with convergence.")
  } else {
    warning("MLE finder failed to converge after ", iter, " iterations.
            The estimates may be meaningless.")
  }
  return(list(
    coef = beta, info_mat = -hess, converged = converged, iter = iter
  ))
}

#' implement BFGS to find MLE under a linear or logistic regression model
mle_bfgs <- function(design, outcome, model, ...) {
  init_coef <- rep(0, ncol(design))
  if (model == "linear") {
    obj_fn <- function(coefs) {
      lm_log_likelihood(design, outcome, coefs, ...)
    }
    obj_grad <- function(coefs) {
      lm_loglik_gradient(design, outcome, coefs, ...)
    }
  }
  if (model == "logit") {
    obj_fn <- function(coefs) {
      logit_log_likelihood(design, outcome, coefs)
    }
    obj_grad <- function(coefs) {
      logit_loglik_gradient(design, outcome, coefs)
    }
  }
  result <- stats::optim(
    par = init_coef, fn = obj_fn, gr = obj_grad,
    method = "BFGS", control = list(fnscale = -1)
  )
  optim_converged <- (result$convergence == 0L)
  if (!optim_converged) {
    warning("Optimization did not converge. The estimates may be meaningless.")
  }
  return(list(coef = result$par, converged = optim_converged))
}
