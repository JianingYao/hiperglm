test_that("return TRUE if MLE via pseudo-inverse is close to truth", {
  set.seed(1234)
  design <- as.matrix(cbind(rep(1, 5), rnorm(5, 0, 1), rnorm(5, 10, 2)))
  truth <- as.matrix(c(5, -4, 2.5))
  outcome <- as.matrix(design %*% truth)
  expect_true(are_all_close(lm_mle_pinv(design, outcome), truth,
    abs_tol = 1e-6, rel_tol = 1e-6
  ))
})

test_that("return TRUE if MLE via BFSG is close to truth", {
  set.seed(1234)
  design <- as.matrix(cbind(rep(1, 5), rnorm(5, 0, 1), rnorm(5, 10, 2)))
  truth <- as.matrix(c(5, -4, 2.5))
  outcome <- as.matrix(design %*% truth)
  expect_true(are_all_close(lm_mle_BFGS(design, outcome, noise_var = 0.1),
    truth,
    abs_tol = 1e-6, rel_tol = 1e-6
  ))
})

test_that("return TRUE if analytical and numerical gradient match", {
  set.seed(1234)
  design <- as.matrix(cbind(rep(1, 5), rnorm(5, 0, 1), rnorm(5, 10, 2)))
  truth <- as.matrix(c(5, -4, 2.5))
  outcome <- as.matrix(design %*% truth)
  n_test <- 10
  grads_are_close <- TRUE
  for (i in 1:n_test) {
    if (!grads_are_close) break
    beta <- rnorm(3)
    analytical_grad <- loglik_gradient(design, outcome, beta)
    numerical_grad <- approx_grad(function(beta) log_likelihood(design, outcome, beta), beta)
    grads_are_close <- are_all_close(analytical_grad, numerical_grad,
      abs_tol = 1e-3, rel_tol = 1e-3
    )
  }
  expect_true(grads_are_close)
})

test_that("return TRUE if MLE via pseudo-inverse and BFGS match", {
  set.seed(1234)
  design <- as.matrix(cbind(rep(1, 5), rnorm(5, 0, 1), rnorm(5, 10, 2)))
  truth <- as.matrix(c(5, -4, 2.5))
  outcome <- as.matrix(design %*% truth)
  pinv_result <- hiper_glm(design, outcome, model = "linear", method = "pinv")
  BFGS_result <- hiper_glm(design, outcome, model = "linear", method = "BFGS")
  expect_true(are_all_close(coef(pinv_result), coef(BFGS_result),
    abs_tol = 1e-3, rel_tol = 1e-3
  ))
})
