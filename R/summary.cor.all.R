#' Summarize Correlations
#' 
#' Extract a data frame that summarizes the correlations computed by
#' \code{cor.all}.
#' 
#' 
#' @param object an object created by the \code{cor.all} function.
#' @param p.adjust a character string describing the method to use to adjust
#' the p-values to account for m ultiple comparisons. See
#' \code{\link{p.adjust}} for the options and details.
#' @param variable if the name of a variable, then summarize only those
#' correlations of this variable with all others, otherwise summarize all
#' combinations.
#' @param \dots further arguments passed to or from other methods.
#' @return A data frame containing columns of the paired variables, the
#' correlation, the adjusted p-value, and the number of observations in the
#' correlation.
#' @seealso \code{\link{cor.all}}
#' @keywords htest
#' @examples
#' 
#' library(smwrData)
#' data(TNLoads)
#' # Extract only the correlations with the log of total nitrogen
#' summary(cor.all(TNLoads[, 1:5]), variable="LOGTN")
#' 
#' @export
#' @method summary cor.all
summary.cor.all <- function(object, p.adjust="holm", variable=NULL, ...) {
	NN <- nrow(object$estimates)
	RN <- matrix(rownames(object$estimates), nrow=NN, ncol=NN)
	CN <- matrix(rep(colnames(object$estimates), each=NN), ncol=NN)
	if(is.null(variable))
		sel <- lower.tri(RN)
	else {
		pck <- which(variable == rownames(object$estimates))
		if(length(pck) != 1)
			stop(variable, " not found in list of variables")
		sel <- RN == variable
		sel[pck, pck] <- FALSE
	}
	retval <- data.frame(Var1=RN[sel], Var2=CN[sel], Cor=object$estimates[sel], 
											 Pval=p.adjust(object$p.values[sel], method=p.adjust),
											 Counts=object$counts[sel], stringsAsFactors=FALSE)
	if(p.adjust != "none")
		names(retval)[4] <- paste("Pval", p.adjust, sep='.')
	names(retval)[3] <- switch(object$call.method,
														 pearson="Cor",
														 spearman="Rho",
														 kendall="Tau")
	return(retval)
}
