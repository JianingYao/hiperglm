#' @export
hiper_glm <- function(design, outcome, model = "linear", method = "BFGS") {
  supported_model <- c("linear", "logit")
  if (!(model %in% supported_model)) {
    stop(sprintf("The model %s is not supported.", model))
  }
  if (model == "linear") {
    if (method == "BFGS") {
      hglm_out <- mle_bfgs(design, outcome, model)
    } else if (method == "pinv") {
      hglm_out <- mle_pinv(design, outcome)
    } else {
      stop(sprintf(
        "The method %s is not supported,
        please use \"BFGS\" or \"pinv\" for linear regression.", method
      ))
    }
  }
  if (model == "logit") {
    if (method == "BFGS") {
      hglm_out <- mle_bfgs(design, outcome, model)
    } else if (method == "newton") {
      hglm_out <- mle_newton(design, outcome, tol = 1e-8)
    } else {
      stop(sprintf(
        "The method %s is not supported,
        please use \"BFGS\" or \"newton\" for logistic regression.", method
      ))
    }
  }
  class(hglm_out) <- "hglm"
  return(hglm_out)
}
