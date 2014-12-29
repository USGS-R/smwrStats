#' Variance Inflation Factors
#' 
#' Computes the variance inflation factor (Helsel and Hirsch, 2002) for each
#' variable in a linear regression fit.
#' 
#' 
#' @aliases vif vif.lm
#' @param model an object of class "lm." There is no default method for
#' \code{vif}.
#' @param \dots further arguments passed to or from other methods.
#' @return A named numeric vector containing the variance inflation factors for
#' each variable.
#' @references Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in
#' water resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' @keywords model
#' @export vif
vif <- function(model, ...) {
	# Coding history:
	#    2005Jul11 DLLorenz Original documented version
	#    2007Apr02 DLLorenz Added option for weighted regression
	#    2007May11 DLLorenz Modified for single explanatory variable
	#    2007Jun05 DLLorenz Bug fix.
	#    2008Jan04 DLLorenz Added computation for GLS model
	#    2011Jul27 DLLorenz Conversion to R
	#    2011Oct25 DLLorenz Update for package
	#    2012Dec12 DLLorenz Split into methods
	#    2013May13 DLLorenz Prevent VIF failure for intercept only model
	#    2014Dec29 DLLorenz Conversion to roxygen header
	UseMethod("vif")
}

#' @rdname vif
#' @export
#' @method vif lm
vif.lm <- function(model, ...) {
  ##
  ## beta needed for names only
  beta <- model$coef
  if(is.matrix(beta)) {
    beta <- beta[, 1]
  }
  na <- is.na(beta)
  beta <- beta[!na]
  lm.x <- model.matrix(model)[,-1, drop=FALSE]
  if(any(na))
    lm.x <- lm.x[, !na, drop = FALSE]
  if(ncol(lm.x) <= 1)
    VIFs <- 1
  else {
    if(is.null(model$weights))
      VIFs <- diag(solve(cor(lm.x)))
    else
      VIFs <- diag(solve(cov.wt(lm.x, wt=model$weights, cor=TRUE)$cor))
  }
  names(VIFs) <- names(beta)[-1]
  return(VIFs)
}
