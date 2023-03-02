#' use Cholesky decomposition to perform pseudo-inverse
lm_mle_pinv <- function(design, outcome){
  a <- crossprod(design)
  b <- crossprod(design, outcome)
  upper <- chol(a)
  z <- backsolve(upper, b, transpose = TRUE)
  beta <- as.vector(backsolve(upper, z))
  return(beta)
}

