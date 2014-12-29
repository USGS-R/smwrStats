#' Print Objects
#' 
#' Print the results of a multiple regression diagnostic analysis
#' (\code{multReg}).
#' 
#' 
#' @param x an object of class "multReg" from \code{multReg}.
#' @param digits the number of significant digits to print numeric data.
#' @param \dots not used for method, required for other methods.
#' @return The object \code{x} is returned invisibly.
#' @note The printed output contains the regression model call; a summary of
#' the residuals; A table of the coefficients with their estiamtes, standard
#' errors, t vlaues, and attained probability levels; the residual standard
#' error; R-squared and F-statistical summaries; Model comparison statistics;
#' if more than one explanatory variable a type II sum-of-squares analysis of
#' variance table and variance inflation factors; and selected test criteria
#' with observations that exceed one of more of the criteria.
#' @export
#' @method print multReg
print.multReg <- function(x, digits=3, ...) {
	##
	print(x$parmests, digits=digits, signif.stars=FALSE)
	## Add model comparison stats:
	cat("press: ", signif(press(x$object), digits),
			"\n  AIC: ", signif(AIC(x$object), digits),
			"\n  BIC: ", signif(BIC(x$object), digits),
			"\n\n", sep="")
	if(length(x$vif) > 1L) {
		print(x$aovtab, digits=digits, signif.stars=FALSE) # really only needed for MLR
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
		dstats <- format(x$diagstats[x$flagobs,] , digits=digits)
		## Append * to each value that exceeds its criterion
		dstats$leverage <- paste(dstats$leverage, 
														 ifelse(x$diagstats[x$flagobs, "leverage"] > x$crit.val[1L],
														 			 "*", " "), sep="")
		dstats$cooksD <- paste(dstats$cooksD, 
													 ifelse(x$diagstats[x$flagobs, "cooksD"] > x$crit.val[2L],
													 			 "*", " "), sep="")
		dstats$dfits <- paste(dstats$dfits, 
													ifelse(abs(x$diagstats[x$flagobs, "dfits"]) > x$crit.val[3L],
																 "*", " "), sep="")
		print(dstats)
	}
	else
		cat("\tNo observations exceeded any test criteria\n")
	invisible(x)
}
