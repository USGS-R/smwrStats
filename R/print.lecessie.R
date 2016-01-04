#' Print Objects
#' 
#' Prints the results of a le Cessie-van Houwelingen test
#' (\code{leCessie.test}).
#' 
#' 
#' @param x an object of class "lecessie" from \code{leCessie.test}.
#' @param digits the number of significant digits to print numeric data.
#' @param \dots not used for method, required for other methods.
#' @return The object \code{x} is returned invisibly.
#' @note The printed ouput is very similar to the printed ouput for class
#' "htest," the output from a hypothesis test, but includes the additional
#' output summarizing the distance between observations, which can be useful
#' for comparing the reults with different bamdwidth settings.
#' @export
#' @method print lecessie
print.lecessie <- function(x, digits=4, ...) {
	##
	x.to.print <- x
	x.to.print$parameters <- round(x.to.print$parameters, digits) # fix this
	oldClass(x.to.print) <- "htest"
	print(x.to.print)
	cat("Distance between observations:\n")
	print(c(maximum=x$max.distance, bandwidth=x$bandwidth))
	cat("\n")
	invisible(x)
}
