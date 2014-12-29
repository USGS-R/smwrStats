#' Prediction Error Sum of Squares
#' 
#' Compute the prediction error sum of squares statistic (PRESS) (Helsel and
#' Hirsch, 2002) for a linear regression model.
#' 
#' 
#' @param model an object of class "lm" or the output from \code{lsfit}.
#' @return The prediction error sum of squares statistic.
#' @seealso \code{\link{multReg}}, \code{\link{lm}}
#' @references Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in
#' water resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' @keywords regression
#' @export press
press <- function(model) {
	# Coding history:
	#    2005Jul14 DLLorenz Initial dated verion
	#    2011Aug09 DLLorenz Conversion to R
	#    2012Nov06 DLLorenz Bug fix to account for influence returning 0 for NAs
	#    2014Dec29 DLLorenz convert to roxygen headers
  ##
  if(inherits(model, "lm")) {
    h <- influence(model)$hat
    r <- resid(model)
    retval <- sum((r/(1 - h))^2, na.rm=TRUE)
  }
  else if(!is.null(model$residuals) && !is.null(model$qr)) {
    ## Assumed lsfit
    h <- ls.diag(model)$hat
    ## lsfit returns the same length and corresponding values even if
    ## there are missing values, so we can compute the correct value!
    retval <- sum((model$residuals/(1 - h))^2, na.rm=TRUE)
  }
  else
    stop("Input must be a 'lm' or 'lsfit' object")
  return(retval)
}
