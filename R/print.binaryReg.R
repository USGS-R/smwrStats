#' Print Objects
#' 
#' Print the results of a logistic regression diagnotic analysis
#' (\code{binaryReg}).
#' 
#' The original regression model should be the output from \code{glm} with
#' \code{binomial} as the value for the \code{family} argument.
#' 
#' @param x an object of class "binaryreg" from \code{binaryReg}.
#' @param digits the number of significant digits to print numeric data.
#' @param \dots not used for method, required for other methods.
#' @return The object \code{x} is returned invisibly.comp2 Description of
#' 'comp2'
#' @note The printed output contains the original call; a summary of the
#' deviance residuals; the coefficients with their estimates, standard errors,
#' z scores, and attained p-values; a summary of the deviance comparison
#' between the model and null model (no explanatory varaibles), also known as
#' the G-squared or -2 log-likelihood; the response profile, which matches the
#' value of the response variable to 0 or 1 and the number of observations in
#' each group; the le Cessie-van Houwelingen and Hosmer-Lemeshow goodness of
#' fit tests; three assessments of the predictinve power: the MaFadden
#' R-squared and the adjusted R-squared, the classification table, the
#' concordance index, and the area under the reciever operating characteristics
#' curve; and selected test criteria with observations that exceed one of more
#' of the criteria.
#' @export
#' @method print binaryreg
print.binaryreg <- function(x, digits=4, ...) {
	##
	## Print the summary, any warning, the factor information, and G2
	print(x$regsum, digits=digits, ...)
	if(x$Warning != "")
		cat(x$Warning)
	if(length(x$Factors) > 0) {
		cat("\nFactor level information:\n")
		lapply(x$Factors, print)
	}
	G2 <- round(x$regsum$null.deviance - x$regsum$deviance, digits)
	df <- x$regsum$df[1]
	cat("\nLikelihood ratio test: ", G2, " on ", df-1,
			" degrees of freedom, p-value is ", round(1-pchisq(G2, df-1), digits),
			"\n\n", sep='')
	## Print the response matrix and the le Cressie and Houwelingen test
	cat("Response profile:\n")
	print(x$Profile, digits=digits, ...)
	cat("\n Goodness of fit tests\n")
	if(!is.null(x$leCessie))
		print(x$leCessie, digits=digits, ...)
	else
		cat("Le Cessie-Van Houwlingen test not computed.\n")
	if(!is.null(x$Hosmer))
		print(x$Hosmer, digits=digits, ...)
	else
		cat("Too few unique predcited values for Hosmer-Lemeshow Test\n")
	## Print the correct and concordant stats and the AUROC
	cat("\nPredictive power estimates:\n")
	## Print the R2 and adjusted R2
	R2 <- round(1 - x$object$deviance/x$object$null.deviance, digits)
	adjR2 <- round(1 - (x$object$aic - 2)/x$object$null.deviance, digits) # need to correct for intercept term
	cat("McFadden R-squared: ", R2, "\nadjusted R-squared: ", adjR2, "\n\n", sep="")
	if(!is.null(x$PctCorrect)) {
		cat("\nClassification table.\nPercent correct: (1 is sensitivity, 0 is specificity)\n")
		print(x$PctCorrect, digits=3)
		nConcord <-  sum(x$Concordance)
		cat("\nConcordance Index, based on",
				nConcord, "pairs\n")
		print(100*x$Concordance/nConcord, digits=digits, ...)
	}
	else
		cat("\nOther predictive power estimates cannot be computed for martix reponses.\n")
	print(x$roc)
	## Print the diagnostics
	cat("\nInfluence diagnostic test criteria:\n")
	print(x$crit.val, digits=digits, ...)
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
