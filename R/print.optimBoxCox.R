#' Print Object
#' 
#' PRint the results of a multivariate unconditional Box-Cox transformation.
#' 
#' @param x an object of class "optimBoxCox" from \code{optimBoxCox}.
#' @param digits the number of significant digits to print numeric data.
#' @param \dots not used for method, required for other methods.
#' @return The object \code{x} is returned invisibly.
#' @note The printed output contains a table showing the power transformation
#'values and their standard errors.
#' @seealso \code{\link[smwrBase]{boxCox}}
#' @export
#' @method print optimBoxCox
print.optimBoxCox<-function(x, digits=4, ...){
	nc <- length(x$lambda)
	if(nc == 1L) {
		cat("Optimized Box-Cox Transformation to Normality\n\n")
	} else
		cat("Optimized Box-Cox Transformations to Multinormality\n\n")
	lambda<-x$lambda
	stderr<-x$stderr
	## round to .1
	rnd <- round(lambda, 1)
	table<-cbind(lambda,stderr, rnd)
	rownames(table)<-x$names
	colnames(table)<-c("Est. Lambda", "Std.Err.", "Rnd. Lambda")
	if(nc == 1L)
		rownames(table)<-""
	print(round(table,digits))
	invisible(x)
}
