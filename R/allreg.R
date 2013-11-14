# All subsets regression problems
# Requires the leaps package.
#
# Coding history:
#    ????????? AVecchia Original coding
#    2007Mar29 DLLorenz Modify argmuents for USGS library
#    2007Mar30 DLLorenz More modifications
#    2007May02 DLLorenz Changed column names in output
#    2007Aug13 DLLorenz Added option to remove missing frm x
#    2011Apr26 DLLorenz Conversion to R
#    2011Oct25 DLLorenz Update argument documentation
#    2011Oct25          This version.
#

allReg <- function(x, y, wt=rep(1,nrow(x)), nmax=ncol(x), nbst=3,
                   na.rm.x=TRUE) {
  ## Arguments:
  ##  x (rectangular numeric) the candidate explanatory variables
  ##  y (numeric matrix or vector) the response variables
  ##  wt (numeric vector) any weighting variable
  ##  nmax (numeric scalar) the maximum number of explanatory variables to
  ##    include in the candidate models
  ##  nbst (numeric scalar) the number of best models to determine for each
  ##    subset size
  ##  na.rm.x (logical scalar) should missing values in x be removed
  ## 
  ## The output object contains the regression output. It is a data.frame
  ##     showing the variables included in the nbst models for each model
  ##     size along with ratings
  ##
  if(!require(leaps))
    stop("The leaps package is required for allreg")
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
