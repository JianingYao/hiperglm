test_that("return TRUE if analytical and numerical gradient match in logistic regression", {
  n_obs <- 32
  n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 1918)
  design <- data$design
  outcome <- data$outcome
  set.seed(615)
  n_test <- 10
  grads_are_close <- TRUE
  for (i in 1:n_test) {
    if (!grads_are_close) break
    beta <- rnorm(n_pred)
    analytical_grad <- logit_loglik_gradient(design, outcome, beta)
    numerical_grad <- approx_grad(function(x) logit_log_likelihood(design, outcome, x), beta)
    grads_are_close <- are_all_close(analytical_grad, numerical_grad,
      abs_tol = 1e-3, rel_tol = 1e-3
    )
  }
  expect_true(grads_are_close)
})

test_that("return TRUE if bfgs coincides with glm outputs on logit model", {
  data <- simulate_data(32, 4, model = "logit", intercept = 1, seed = 1918)
  glm_out <- stats::glm(data$outcome ~ data$design + 0, family = binomial())
  hglm_out <- hiper_glm(data$design, data$outcome,
    model = "logit", method = "BFGS"
  )
  expect_true(are_all_close(coef(hglm_out), coef(glm_out),
    abs_tol = 1e-3, rel_tol = 1e-3
  ))
})

test_that("return TRUE if analytical and numerical Hessian matrix match in logistic regression", {
  n_obs <- 32
  n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 1918)
  design <- data$design
  outcome <- data$outcome
  set.seed(615)
  n_test <- 10
  grads_are_close <- TRUE
  for (i in 1:n_test) {
    if (!grads_are_close) break
    v <- rnorm(n_pred)
    beta <- rnorm(n_pred)
    eps <- 1e-8
    hess_v_approx <- (logit_loglik_gradient(design, outcome, beta + eps * v) - logit_loglik_gradient(design, outcome, beta - eps * v)) / (2 * eps)
    hess_v <- logit_hessian(design, beta) %*% v
    grads_are_close <- are_all_close(as.vector(hess_v_approx), as.vector(hess_v),
      abs_tol = 1e-3, rel_tol = 1e-3
    )
  }
  expect_true(grads_are_close)
})

test_that("return TRUE if newton coincides with bfgs and glm outputs on logit model", {
  n_obs <- 32
  n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 100)
  design <- data$design
  outcome <- data$outcome
  via_newton_lu_out <- hiper_glm(design, outcome,
    model = "logit", method = "newton", newton_qr = FALSE
  )
  via_bfgs_out <- hiper_glm(design, outcome,
    model = "logit", method = "BFGS"
  )
  via_newton_qr_out <- hiper_glm(design, outcome,
    model = "logit", method = "newton", newton_qr = TRUE
  )
  glm_out <- stats::glm(data$outcome ~ data$design + 0, family = binomial())
  expect_true(are_all_close(
    coef(via_newton_lu_out), coef(via_bfgs_out),
    abs_tol = 1e-3, rel_tol = 1e-3
  ))
  expect_true(are_all_close(
    coef(glm_out), coef(via_newton_lu_out),
    abs_tol = 1e-3, rel_tol = 1e-3
  ))
  expect_true(are_all_close(
    coef(via_newton_qr_out), coef(via_newton_lu_out),
    abs_tol = 1e-3, rel_tol = 1e-3
  ))
})
