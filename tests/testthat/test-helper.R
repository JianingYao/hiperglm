# source("./tests/testthat/helper.R")
test_that("are_all_close returns TRUE when inputs are close", {
  v <- rep(100, 100)
  w <- v + 1e-8
  expect_true(are_all_close(v, w, abs_tol = 1e-6, rel_tol = 1e-6))
}
)


test_that("are_all_close returns FALSE when the relative error is above rel_tol", {
  v <- rep(100, 100)
  w <- v + 1e-2
  expect_false(are_all_close(v, w, abs_tol = 1e-6, rel_tol = 1e-3))
}
)


test_that("are_all_close returns FALSE when the absolute error is above abs_tol", {
  v <- rep(100, 100)
  w <- v + 1e-3
  expect_false(are_all_close(v, w, abs_tol = 1e-3, rel_tol = 1e-6))
}
)
