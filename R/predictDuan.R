#' @include predictMVUE.R
#' @rdname predictMVUE
#' @export
predictDuan <- function(object, newdata, back.trans=exp) {
	# Coding history:
	#    2013Aug13 DLLorenz Original coding
	#    2014Dec29 DLLorenz Conversion to roxygen headers
	#
  if(missing(newdata))
    newdata <- eval(as.list(object$call)$data)
  firstguess <- predict(object, newdata, type = "response")
  resids <- residuals(object, type = "response")
  ## Clean up just in case there are NAs in the residuals
  resids <- resids[!is.na(resids)]
  ## Compute the bias correction factor for each observation
  ## Note that only for the log transform can the back transorm be 
  ## computed for each observation
  retval <- sapply(firstguess, function(x) 
    mean(back.trans(x + resids)))
  return(retval)
}
