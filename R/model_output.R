#' @export
coef.hglm <- function(hglm_out) {
  coef <- hglm_out$coef
  return(coef)
}

#' @export
vcov.hglm <- function(hglm_out) {
  warning("The function is yet to be implemented.")
  cat("The variance-covariance matrix is ...")
}

#' @export
print.hglm <- function(hglm_out) {
  warning("The function is yet to be implemented.")
  cat("The model is ...")
}
