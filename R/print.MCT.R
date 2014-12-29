#' Print Objects
#' 
#' Print the results of a multiple comparison test (\code{multicomp.test}).
#' 
#' 
#' @param x an object of class "MCT" from \code{multicomp.test}.
#' @param digits the number of significant digits to print numeric data.
#' @param \dots not used for method, required for other methods.
#' @return The object \code{x} is returned invisibly.
#' @note The printed output contains a description of the test, critical
#' values, the variables in the test, and two tables: the paired comparisons
#' and associations among the groups. The table of the paried comparisons shows
#' the groups in the comparison, the estimate of the difference between the
#' group means, the standard error of the difference, lower and upper
#' confidence intervals, and a flag that indicates if the confidence interval
#' excludes 0, which indicates wheter the difference is significantly different
#' from 0 at the user-specified value. The table of asociations shows the
#' group, the mean value of the response, the number of observations in the
#' group, and any number of columns names "A," "B," and so forth that represent
#' possible associations of the groups whare an "X" is present in the group.
#' @export
#' @method print MCT
print.MCT <- function(x, digits=4, ...) {
	## print function for objects of class MCT
	cat("\t", x$title, "\n")
	if(x$cv.method == "lsd")
		cat("Pairwise")
	else
		cat("Overall")
	cat(" error rate: ", x$alpha, "\nCritical value: ", round(x$crit.value, digits),
			" by the ", x$cv.method, " method\n\n", sep="")
	cat("Response variable: ", x$response, "\nGroup variable: ", x$groups,
			"\n\n", sep="")
	cat("Table of paired comparisons, ", round(1 - x$alpha, 4) * 100,
			" percent confidence intervals\n excluding 0 are flagged by *.\n", sep="")
	pmat <- format(signif(x$table, digits))
	pmat <- cbind(pmat, flag=ifelse(x$table[,3]*x$table[,4] > 0, "*", " "))
	print(pmat, quote=FALSE)
	cat("\nTable of associations among groups\n")
	pmat=cbind(Mean=signif(x$means, digits), Size=x$sizes, x$assoc)
	print(pmat, quote=FALSE)
	cat("\n")
	invisible(x)
}
