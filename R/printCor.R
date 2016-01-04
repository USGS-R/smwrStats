#' Print Objects
#' 
#' Prints the results of a correlation test (\code{cor.all}) or correlation matrix
#'in a compact form indicating positive or negative correlation.
#' 
#' 
#' @param x an object of class "cor.all" from \code{cor.all} or the output from \code{cor}.
#' @param criterion a numeric value indicating the test criterion for showing a value. 
#' @return The object \code{x} is returned invisibly.
#' @note The printed output is a compressed table showing "+" where the value in \code{x}
#'is greater than \code{criterion}, "-" where the value in \code{x} is less than 
#'-\code{criterion}, "\" on the diagonal, if \code{x} is symmetric, and "." where \code{x} 
#'has a missing value, and " " otherwise.
#' @seealso \code{\link{cor.all}}
#' @examples
#' \dontrun{
#' library(smwrData)
#' data(TNLoads)
#' printCor(cor(TNLoads, method="spearman"), .5)
#' }
#' @export
printCor <- function(x, criterion=0.75) {
	## Get the matrix
	if(mode(x) == "list") {
		mat <- x$estimates
	} else {
		mat <- x
	}
	# Format the column names to a single column
	Names <- colnames(mat)
	if(is.null(Names)) {
		Names <- seq(ncol(mat))
	}
	Names <- format(Names, justify="right")
	ColNames <- do.call(cbind, strsplit(Names, split=""))
	# Format the row names to a single column
	Names <- rownames(mat)
	if(is.null(Names)) {
		Names <- seq(nrow(mat))
	}
	Names <- format(Names, justify="right")
	ColBuff <- rep(" ", nchar(Names[1L]))
	PrMat <- matrix(" ", ncol=ncol(mat), nrow=nrow(mat))
	PrMat[mat > criterion] <- "+"
	PrMat[mat < -criterion] <- "-"
	# if square and symmetric add the diagonal
	if(ncol(mat) == nrow(mat) && max(abs(mat - t(mat)), na.rm=TRUE) < 1.e-9) {
		PrMat[cbind(seq(nrow(mat)), seq(ncol(mat)))] <- "\\"
	}
	PrMat[is.na(mat)] <- "."
	# Print the results:
	for(i in seq(along=ColBuff)) { # just need the count
		cat(c(ColBuff, ColNames[i,], "\n"), sep="")
	}
	for(i in seq(nrow(PrMat))) {
		cat(c(Names[i], PrMat[i,], "\n"), sep="")
	}
	invisible(x)
}
