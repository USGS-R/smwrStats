#' Print Objects
#' 
#' Print the results of a correlation test (\code{cor.all}).
#' 
#' 
#' @param x an object of class "cor.all" from \code{cor.all}.
#' @param digits the number of significant digits to print numeric data.
#' @param lower logical, print only the lower triangular matrix? Otherwise,
#' print the full, square matrix.
#' @param \dots not used for method, required for other methods.
#' @return The object \code{x} is returned invisibly.
#' @note The printed output contains a description of the test; and 3 lines for
#' each comparison that is printed: the correlation statistic, the attained
#' p-value, and the number of pairs.
#' @export
#' @method print cor.all
print.cor.all <- function(x, digits = 4, lower=TRUE, ...) {
	##
	cat(paste("\n\t", x$method, "\n\n", sep = ""))
	cat("data: ", x$data.name, "\n\n")
	if(x$distribution == "lognormal") {
		cat("All data were log-transformed\n")
	} else if(x$distribution == "log1p") 
		cat("All data were log1p-transformed\n")
	est <- format(x$estimates, digits=digits)
	N <- ncol(est)
	val <- round(x$p.values, digits=digits)
	val0 <- (val == 0) & (abs(x$estimates) < 1)
	val <- format(val, width=digits+3, justify="right")
	if(any(val0)) {
		vallt <- paste("<0.", paste(rep(0, digits-1), collapse=''), "1", sep='')
		val[val0] <- vallt
	}
	val <- matrix(val, N, N) # Need to restore dim
	cnt <- format(x$counts)
	if(lower) {
		est <- est[-1, -N, drop=FALSE]
		val <- val[-1, -N, drop=FALSE]
		cnt <- cnt[-1, -N, drop=FALSE]
		N <- N - 1
		if(N > 1) {
			for(i in seq(1, N-1))
				for(j in seq(i+1, N)) {
					est[i, j] <- " "
					val[i, j] <- " "
					cnt[i, j] <- " "
				}
		}
	}
	outmat <- rbind(est[1,], val[1,], cnt[1,], ' ')
	if(N > 1)
		for(i in seq(2, N))
			outmat <- rbind(outmat, est[i,], val[i,], cnt[i,], ' ')
	## Make row names descriptive
	dn <- dimnames(est)
	dn[[1]] <- rep(dn[[1]], each=4)
	dn[[1]] <- switch(x$call.method,
										pearson=paste(dn[[1]], c("-cor", "-r!=0", "-N", "--"), sep=''),
										spearman=paste(dn[[1]], c("-rho", "-r!=0", "-N", "--"), sep=''),
										kendall=paste(dn[[1]], c("-tau", "-t!=0", "-N", "--"), sep=''))
	## Blank line--blank row name
	dn[[1]] <- sub(".*--$", "", dn[[1]])
	dimnames(outmat) <- dn
	print(outmat, quote=FALSE)
	invisible(x)
}
