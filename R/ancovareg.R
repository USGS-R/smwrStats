#' Diagnostics for Analysis of Covariance
#' 
#' Computes diagnostic statistics for an analysis of covariance (ANCOVA) with a
#' single factor variable.
#' 
#' @details The input model object (\code{object}) can be the complete ancova model
#' including all interaction terms or it can be any form of an ANCOVA model.
#' Most often, if it is not the complete ancova model, then \code{find.best}
#' should be \code{FALSE}.\cr The \code{find.best} option uses the \code{step}
#' function to select the "best" subset of terms in the model. In general, this
#' can be used to retain or drop significant interaction terms. It will not
#' look at individual factor levels in the model.
#' 
#' @param object the linear regression model object.
#' @param find.best select the "best" subset of terms in the model?
#' @param trace print the results of the selection process if \code{find.best}
#' is \code{TRUE}?
#' @return A list of class "ancovaReg" containing these components:
#' \item{aovtab}{the analysis of variance table from the original model}
#' \item{parmests}{a summary of the final \code{object}.} \item{vif}{a named
#' vector of variance inflation factors.} \item{diagstats}{a data.frame
#' containing the observed values, predicted values, residuals, standardized
#' residuals, studentized residuals, leverage, Cook's D, and dfits for each
#' observation.} \item{crit.val}{a named vector of the critical values for
#' leverage, Cook's D, and dfits.} \item{flagobs}{a logical vector indicating
#' which observations exceeded at least one of the critical values.}
#' \item{object}{the \code{lm} object.} \item{x}{the model matrix of
#' explanatory variables.} \item{factor.var}{the name of the factor variable}
#' \item{x.fr}{the model frame of explanatory variables.} If no factor variable
#' is found in the final model, either because one was not specified or it was
#' dropped from the model, then an object of class "multReg" is returned
#' instead. See \code{\link{multReg}} for details.
#' @note Objects of class "ancovaReg" have \code{print} and \code{plot}
#' methods.
#' @seealso \code{\link{lm}}, \code{\link{plot.ancovaReg}},
#' \code{\link{multReg}},
#' @references Draper, N.R. and Smith, H., 1998, Applied Regression Analysis,
#' (3rd ed.): New York, Wiley, 724 p.
#' @keywords models regression
#' @importFrom car Anova
#' @export ancovaReg
ancovaReg <- function(object, find.best=TRUE, trace=FALSE) {
	# Coding history:
	#    2012Jul20 DLLorenz Original coding
	#    2012Aug28 DLLorenz Change from diagPlot to plot
	#    2013Apr09 DLLorenz Added setGD to plot
	#    2013May10 DLLorenz changed anova to Anova
	#    2014Dec22 DLLorenz Roxygen header
  ##
  aovtab <- Anova(object) # Anova from car does type II
  if(find.best) {
    if(trace)
      cat("Stepwise elimintation of terms from original model\n")
    object <- step(object, direction="backward", trace=as.integer(trace))
  }
  ## Get the factor variable
  factor.var <- attr(object$terms, "dataClasses")
  factor.var <- names(factor.var)[factor.var == "factor"]
  regsum <- summary(object)
  vif <- vif(object)
  ## Compute the diagnostics using lsfit
  ## get x matrix
  if(is.null(object[['x']])) # to avoid conflict with the xlevels component.
    x <- model.matrix(object)[, -1L, drop=FALSE]
  else
    x <- object$x[, -1L, drop=FALSE]
  ## get y
  if(is.null(object$y))
    y <- model.extract(model.frame(object), "response")
  else
    y <- object$y
  respvar <- make.names(object$terms[[2L]])
  if(length(respvar) > 1L) # happens when function is used
    respvar <- paste(respvar, collapse='.')
  ## get weights
  if(is.null(object$weights))
    wt <- rep(1, length(y))
  else
    wt <- object$weights
  ##  Run lsfit and ls.diag to get the residual diagnostic stats we need:
  lsout <- lsfit(x, y, wt)
  lsdiag <- ls.diag(lsout)
  lev <- lsdiag$hat
  cooksd <- lsdiag$cooks
  std.res <- lsdiag$std.res
  dfits <- lsdiag$dfits
  stud.res <- lsdiag$stud.res
  ##  Compute critical values for diagnostic statistics and pick out offending
  ##  observations
  p <- ncol(x)+1
  n <- nrow(x)
  cvlev <- 3*p/n
  cvdfit <- qgrubbs(0.01, n)*sqrt(p/n)
  cvcook <- qf(.5,p+1,n-p)
  pck <- c(lev>cvlev | cooksd>cvcook | abs(dfits)>cvdfit)
  ##  Combine the diagnostic stats into a single data set, round it to 3
  ##  decimals, rename y, and combine the output into a list.
  yhat <- object$fitted.values
  diagstats <- data.frame(y=y, yhat=yhat, resids=lsout$residuals,
                          stnd.res=std.res, stud.res=stud.res,
                          leverage=lev, cooksD=cooksd, dfits=dfits)
  names(diagstats)[1L] <- respvar
  cvs <- c(cvlev,cvcook,cvdfit)
  names(cvs) <- c('leverage','cooksD','dfits')
  if(length(factor.var) > 0L) {
    stuff.out <- list(aovtab=aovtab, parmests=regsum, vif=vif,
                      diagstats=diagstats, crit.val=cvs, flagobs=pck,
                      object=object, x=x, factor.var=factor.var,
                      x.fr=model.frame(object)[, -1L, drop=FALSE]) # drop response
    oldClass(stuff.out) <- "ancovaReg"
  }
  else {
    stuff.out <- list(aovtab=aovtab, parmests=regsum, vif=vif,
                      diagstats=diagstats, crit.val=cvs, flagobs=pck,
                      object=object, x=x)
    oldClass(stuff.out) <- "multReg" # No factors, default to MLR
  }
  return(stuff.out)
}
