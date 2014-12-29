#' All Subsets Regression
#' 
#' Create a table of the best subsets of explantory variables for a response
#' variable.
#' 
#' @include USGSwsStats-package.R
#' @param x matrix of candidate exmplanatory variables.
#' @param y the response variable.
#' @param wt the weight variable if needed.
#' @param nmax the maximum number of explanatory variables to include in the
#' largest model.
#' @param nbst the number of best models to determine for each subset size.
#' @param na.rm.x should missing values in x should be removed?
#' @return A data frame containing these columns: \item{model.formula}{the
#' subset model formula} \item{nvars}{the size (number of variables in the
#' subset model} \item{stderr}{the standard error of the subsbet model}
#' \item{R2}{the coefficient of determination for the subset model}
#' \item{adjr2}{the adjusted r-squared of the subbset model} \item{Cp}{Mallow's
#' Cp for the subset model} \item{press}{the press statistic for the subset
#' model}
#' @note This function is a wrapper for the \code{regsubsets} function in the
#' leaps package.
#' @seealso \code{\link{regsubsets}}
#' @references Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in
#' water resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' 
#' Mallow, C.L., 1973, Come comments of Cp: Technometrics, v. 15, p.
#' 661--675.\cr
#' 
#' Miller, A.J., 1990, Subset selection in regression in Monographs on
#' Statistics and Applied Probability 40: London, Chapman and Hall.\cr
#' @keywords models regression
#' @examples
#' 
#' # See the regression vignette for examples
#' .pager <- options("pager")
#' options(pager="console")
#' vignette(package="USGSwsStats")
#' options(.pager)
#' @importFrom leaps regsubsets
#' @export allReg
allReg <- function(x, y, wt=rep(1,nrow(x)), nmax=ncol(x), nbst=3,
                   na.rm.x=TRUE) {
	# Coding history:
	#    ????????? AVecchia Original coding
	#    2007Mar29 DLLorenz Modify argmuents for USGS library
	#    2007Mar30 DLLorenz More modifications
	#    2007May02 DLLorenz Changed column names in output
	#    2007Aug13 DLLorenz Added option to remove missing frm x
	#    2011Apr26 DLLorenz Conversion to R
	#    2011Oct25 DLLorenz Update argument documentation
	#    2014Dec22 DLLorenz Roxygen headers
  ##
  if(is.null(dimnames(y)))
    yname <- deparse(substitute(y))
  else
    yname <- dimnames(y)[[2]]
  ## Check for missing values in y
  sel <- !is.na(y)
  x <- as.matrix(x)
  if(na.rm.x)
    sel <- sel & apply(x, 1, function(xx) all(!is.na(xx)))
  y <- y[sel]
  x <- x[sel,]
  if(length(wt) > nrow(x))
    wt <- wt[sel] # protects against default
  xnames <- dimnames(x)[[2]]
  n <- nrow(x)
  p <- ncol(x)
  if(any(is.na(c(x, wt)))) {
    missings <- apply(x, 2, FUN= function(xx) sum(is.na(xx)))
    missings <- c(missings, sum(is.na(wt)))
    return(data.frame(vars=c(xnames, "wt"), number.missing=missings))
  }
  ##  Run regsubsets to select the best models:
  stepout <- summary(regsubsets(x, y, wt, method='exhaustive', nvmax=nmax, nbest=nbst))
  ##  Compute PRESS statistics
  which <- stepout$which[,-1] # drop intercept term
  size <- rowSums(which)
  ##  Run lsfit for each model to get the "hat" values and compute the PRESS
  press <- rep(0, length(size))
  for (i in seq(length(size))) {
    tmpout <- lsfit(x[,which[i,]], y, wt)
    tmpres <- tmpout$res*sqrt(wt)
    tmphat <- hat(x[,which[i,]])
    press[i] <- sum((tmpres/(1-tmphat))^2)
  }
  nmp <- n-size-1
  stderr <- stepout$rss / nmp
  xnames <- apply(which, 1, function(xx, xnames) paste(xnames[xx == 1], collapse=' + '), xnames=xnames)
  xnames <- paste(yname, xnames, sep=' ~ ')
  stuff.out <- data.frame(model.formula=xnames, nvars=size, stderr=sqrt(stderr),
                          R2=stepout$rsq*100, adjr2=stepout$adjr2*100,
                          Cp=stepout$cp, press=press, stringsAsFactors=FALSE)
  return(stuff.out)
}
