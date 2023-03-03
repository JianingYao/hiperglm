#' @export
hiper_glm <- function(design, outcome, model = "linear", method = "BFGS") {
  supported_model <- c("linear", "logit")
  if (!(model %in% supported_model)) {
    stop(sprintf("The model %s is not supported.", model))
  }
  if (model == "linear") {
    hglm_out <- list()
    class(hglm_out) <- "hglm"
    if (method == "pinv") {
      hglm_out$coef <- lm_mle_pinv(design, outcome)
    } else if (method == "BFGS") {
      hglm_out$coef <- lm_mle_BFGS(design, outcome)
    } else {
      stop("The function is yet to be implemented for other methods other than pseudo-inverse or BFGS.")
    }
  }
  return(hglm_out)
}
