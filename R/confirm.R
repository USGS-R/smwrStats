#' Confirm an Analysis
#' 
#' Review and accept the results of an analysis
#' 
#' 
#' @aliases confirm confirm.default
#' @param x the object to be confirmed
#' @param \dots additional arguments required for specific methods
#' @return The object returned depends on the specific method.
#' @note The default method simply returns the object and issues a warning.
#' @keywords manip
#' @export confirm
confirm <- function(x, ...)
	# Coding history:
	#    2011Jul06 DLLorenz Original coding
	#    2014Dec22 DLLorenz Roxygen headers
  UseMethod("confirm")

#' @rdname confirm
#' @export
#' @method confirm default
confirm.default <- function(x, ...) {
  warning(paste("No known confirm method for object of class ", class(x), sep=''))
  invisible(x)
}
