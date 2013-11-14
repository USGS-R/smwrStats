# correct log-transformed response variable using Ferguson (MLE)
#
# Coding history:
#    2013Aug13 DLLorenz Original coding

predictFerguson <- function(object, newdata, Log10=FALSE) {
  Fact <- if(Log10) 2.30258509299405 else 1.
  if(missing(newdata))
    newdata <- eval(as.list(object$call)$data)
  firstguess <- predict(object, newdata, type = "response")
  ## Compute the bias correction factor
  BCF <- exp(0.5*(rmse(object)*Fact)^2)
  retval <- exp(firstguess * Fact) * BCF
  return(retval)
}
