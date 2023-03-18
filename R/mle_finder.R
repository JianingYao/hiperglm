#' use QR decomposition to solve a linear model
mle_qr <- function(design, outcome) {
  beta <- qr.solve(design, outcome)
  return(list(coef = beta))
}

#' perform each iteration of Newton's method to solve a logit model
take_one_newton_step <- function(design, outcome, beta, qr_solver = TRUE) {
  if (qr_solver) {
    est <- design %*% beta
    mu <- sigmoid(est)
    weight <- logit_weights(design, beta)
    x_tilde <- diag(sqrt(weight)) %*% design
    y_tilde <- sqrt(weight) * est + 1 / sqrt(weight) * (outcome - mu)
    beta <- qr.solve(x_tilde, y_tilde)
  } else {
    grad <- logit_loglik_gradient(design, outcome, beta)
    hess <- logit_hessian(design, beta)
    beta <- beta - solve(hess, grad)
  }
  return(beta)
}

#' implement Newton's method to find MLE under a logistic regression model
mle_newton <- function(design, outcome, tol = NULL, qr_solver = TRUE, max_iter = 50) {
  n_pred <- ncol(design)
  if (is.null(tol)) {
    tol <- 1e-5
  }
  beta <- rep(0, n_pred)
  log_lik_old <- logit_log_likelihood(design, outcome, beta)
  iter <- 0
  converged <- FALSE
  while ((iter < max_iter) && (!converged)) {
    beta <- take_one_newton_step(design, outcome, beta, qr_solver = qr_solver)
    log_lik_new <- logit_log_likelihood(design, outcome, beta)
    converged <- are_all_close(log_lik_old, log_lik_new,
      abs_tol = tol, rel_tol = tol
    )
    iter <- iter + 1
    log_lik_old <- log_lik_new
  }
  if (converged) {
    cat("MLE is found after ", iter, " iterations with convergence. \n")
  } else {
    warning("MLE finder failed to converge after ", iter, " iterations.
            The estimates may be meaningless. \n")
  }
  return(list(
    coef = beta, converged = converged, iter = iter
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
