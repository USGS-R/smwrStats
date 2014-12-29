#' Print Object
#' 
#' PRint the results of a receiver operator characteristics (ROC) for a logistic
#' regression model.
#' 
#' @param x an object of class "roc" from \code{roc}.
#' @param digits the number of significant digits to print numeric data.
#' @param \dots not used for method, required for other methods.
#' @return The object \code{x} is returned invisibly.
#' @note The printed output contains the area under to ROC curve.
#' @export
#' @method print roc
print.roc <- function(x, digits=3, ...) {
	## Arguments:
	##  x (roc object) the object to print
	##  digits (integer scalar) how many digits to use when prining
	##   ... (dots) not used, required for method function
	##
	cat("\nArea under the ROC curve: ", round(x$c.val, digits), "\n\n", sep='')
	invisible()
	return(x)
}