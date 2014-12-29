#' Utility Function for Safe Prediction
#' 
#' A utility to help \code{\link{model.frame.default}} create the right matrices 
#'when predicting from models with \code{trends} terms. Used only internally.
#' 
#' @rdname makepredictcall
#' @param var a variable.
#' @param call the term in the formula, as a call.
#' @return A replacement for \code{call} for the prediction variable.
#' @importFrom stats makepredictcall
#' @export
#' @method makepredictcall trends
makepredictcall.trends <- function(var, call) {
	# Coding history:
	#    2013Aug15 DLLorenz Original coding
	#    2014Dec22 DLLorenz Roxygen header
	#
  if (as.character(call)[1L] != "trends") 
    return(call)
  call$breaks <- attr(var, "breaks")
  call$boundary.breaks <- c(TRUE, TRUE) 
  return(call)
}