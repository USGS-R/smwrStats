#' Diagnostics for Linear Regression
#' 
#' Computes diagnostics for linear regression
#' 
#' 
#' @param object the linear regression model object
#' @return A object of class "multReg" having the following components:\cr
#' 
#' \item{aovtab}{ the analysis of variance table, using the type II sum of
#' squares. See \bold{Note}. } \item{parmests}{ a summary of \code{object}.  }
#' \item{vif}{ a named vector of variance inflation factors. }
#' \item{diagstats}{ a data.frame containing the observed values, predicted
#' values, residuals, standardized residuals, studentized residuals, leverage,
#' Cook's D, and dfits for each observation. } \item{crit.val}{ a named vector
#' of the critical values for leverage, Cook's D, and dfits.  See \bold{Note} }
#' \item{flagobs}{ a logical vector indicating which observations exceeded at
#' least one of the critical values. } \item{object}{ the \code{lm} object. }
#' \item{x}{ the model matrix of explanatory variables. }
#' @note The type II sum of squares are calculated according to the principle
#' of marginality, testing each term after all others, except ignoring the
#' term's higher-order relatives. This type sum of squares is useful for
#' assessing the overall marginal effect of each term in the model.\cr The
#' critical values for the test criteria are computed as: leverage,
#' \emph{3p/n}; Cook's D, median quantile for the \emph{F} distribution with
#' \emph{p+1} and \emph{n-p} degrees of freedonm; and dfits, the .01 quantile
#' of the \emph{grubbs} distribution for \emph{n} observations, where \emph{p}
#' is the number of parameters estiamted in the regression and \emph{n} is the
#' number of observations.\cr Objects of class "multReg" have \code{print} and
#' \code{plot} methods.
#' @seealso \code{\link{lm}}, \code{\link{plot.multReg}},
#' @references Draper, N.R. and Smith, H., 1998, Applied Regression Analysis,
#' (3rd ed.): New York, Wiley, 724 p.\cr
#' @keywords models regression
#' @export multReg
multReg <- function(object) {
	# Coding history:
	#    ????????? AVecchia Original coding for Env. Stats. class
	#    2007Apr04 DLLorenz Modified for library
	#    2007Jun05 DLLorenz Added span argument for plot function
	#    2008Apr21 DLLorenz changed comment from parests to parmests
	#    2008Apr24 DLLorenz Bug fix
	#    2008Aug14 DLLorenz Added Wooding's test to plot
	#    2008Oct20 DLLorenz Added PPCC test for normality
	#    2010Nov12 DLLorenz Conversion to R
	#    2011Apr26 DLLorenz Added polynomial nonlinearity test to partial plots
	#    2011Jun29 DLLorenz Bug fix and name change
	#    2011Oct25 DLLorenz Update for package
	#    2012Mar22 DLLorenz Bug fix in diagnostic plots for the correlogram
	#    2012Aug28 DLLorenz Change from diagPlot to plot
	#    2012Dec27 DLLorenz Minor tweaks, including integer indexes
	#    2013Jan08 DLLorenz Used Anova in car rather than anova for type II Sum Sq.
	#    2013Apr09 DLLorenz Added setGD to plot
	#    2013May10 DLLorenz Use na.omit to work with both na.omit and na.exclude in call
	#    2013Dec30 DLLorenz Added * to print of observations exceeding test criterion
	#    2014May14 DLLorenz Increased criterion for dfits
	#    2014Dec22 DLLorenz Roxzygen headers
	##
  ## Required USGS core functions: ppcc.test, vif
  ## Anova already imported from car by ancovareg
  aovtab <- Anova(object) # Anova from car does type II
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
  stuff.out <- list(aovtab=aovtab, parmests=regsum, vif=vif,
                    diagstats=diagstats, crit.val=cvs, flagobs=pck,
                    object=object, x=x)
  oldClass(stuff.out) <- "multReg"
  return(stuff.out)
}
