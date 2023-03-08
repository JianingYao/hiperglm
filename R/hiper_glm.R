#' @export
hiper_glm <- function(design, outcome, model = "linear", method = "BFGS") {
  supported_model <- c("linear", "logit")
  if (!(model %in% supported_model)) {
    stop(sprintf("The model %s is not supported.", model))
  }
  if (model == "linear") {
    if (method == "pinv") {
      MLE <- lm_mle_pinv(design, outcome)
    } else if (method == "BFGS") {
      MLE <- lm_mle_BFGS(design, outcome)
    } else {
      stop("The function is yet to be implemented for other methods other than
           pseudo-inverse or BFGS for linear regression.")
    }
  }
  if (model == "logit") {
    hglm_out <- list()
    class(hglm_out) <- "hglm"
    if (method == "BFGS") {
      MLE <- logit_mle_BFGS(design, outcome)
    } else if (method == "newton") {
      MLE <- logit_mle_newton(design, outcome, tol = 1e-8)
    } else {
      stop("The function is yet to be implemented for other methods other than
           Newton's method or BFGS for logistic regression.")
    }
  }
  hglm_out <- list(coef = MLE)
  class(hglm_out) <- "hglm"
  return(hglm_out)
}
