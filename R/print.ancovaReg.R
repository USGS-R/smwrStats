#' Print Objects
#' 
#' Print the results of a analysis of covariance (\code{ancovaReg}).
#' 
#' 
#' @param x an object of class "ancovaReg" from \code{ancovaReg}.
#' @param digits the number of significant digits to print numeric data.
#' @param \dots not used for method, required for other methods.
#' @return The object \code{x} is returned invisibly.
#' @note The printed output contains the ANOVA table for the orignal models,
#' the regression summary for the final model, variance inflation factors for
#' each explanatory variable in the final model, and selected test criteria
#' with observations that exceed one of more of the criteria.
#' @export
#' @method print ancovaReg
print.ancovaReg <- function(x, digits=3, ...) {
	##
	cat("\nOriginal model\n")
	print(x$aovtab, digits=digits, signif.stars=FALSE)
	cat("\n\n\nFinal model\n")
	print(x$parmests, digits=digits, signif.stars=FALSE)
	if(length(x$vif) > 1) {
		cat("\nVariance inflation factors\n")
		namvif <- format(names(x$vif), justify="left")
		valvif <- format(round(x$vif, 2), justify="right")
		for(i in seq(along=x$vif))
			cat(namvif[i], " ", valvif[i], "\n", sep="")	
	}
	cat("\nTest criteria\n")
	print(x$crit.val, digits=digits)
	if(any(x$flagobs)) {
		cat("\tObservations exceeding at least one test criterion\n")
		print(x$diagstats[x$flagobs,], digits=digits)
	}
	else
		cat("\tNo observations exceeded any test criteria\n")
	invisible(x)
}
