#' @include predictMVUE.R
#' @rdname predictMVUE
#' @export
predictFerguson <- function(object, newdata, Log10=FALSE) {
	# Coding history:
	#    2013Aug13 DLLorenz Original coding
	#    2014Dec29 DLLorenz Conversion to roxygen headers
	#
  Fact <- if(Log10) 2.30258509299405 else 1.
  if(missing(newdata))
    newdata <- eval(as.list(object$call)$data)
  firstguess <- predict(object, newdata, type = "response")
  ## Compute the bias correction factor
  BCF <- exp(0.5*(rmse(object)*Fact)^2)
  retval <- exp(firstguess * Fact) * BCF
  return(retval)
}
