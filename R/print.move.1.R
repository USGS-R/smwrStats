#' Print Objects
#' 
#' Print the results of a move.1 analysis (\code{move.1}).
#' 
#' 
#' @param x an object of class "move.1" from \code{move.1}.
#' @param digits the number of significant digits to print numeric data.
#' @param \dots additonal arguments for printing numeric values.
#' @return The object \code{x} is returned invisibly.
#' @note The printed output contains the call, the regression coefficients, and
#' selected statistics of the variables.
#' @export
#' @method print move.1
print.move.1 <- function(x, digits=4, ...) {
	##
	cat("Call:\n")
	dput(x$call)
	cat("\nCoefficients:\n")
	print(x$coefficients, digits=digits, ...)
	cat("\nStatistics of the variables:\nResponse (", x$var.names[1], "):\n",
			sep="")
	print(x$ystats, digits=digits, ...)
	cat("Predictor (", x$var.names[2], "):\n", sep="")
	print(x$xstats, digits=digits, ...)
	cat("Correlation coefficient: ", round(x$R, digits), 
			"\n                p-value: ", round(x$p.value, digits), "\n", sep="")
	if(!is.null(x$na.action)) {
		n.na <- length(x$na.action)
		if(n.na > 1)
			cat("  (", n.na, " observations deleted due to missing values)\n", sep="")
		else
			cat("  (", n.na, " observation deleted due to missing values)\n", sep="")
	}
	if(!is.null(x$cx) && !is.null(x$cy))
		cat(sum(x$cx | x$cy), " observations were left-consored\n", sep="")
	invisible(x)
}
