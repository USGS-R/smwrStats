#' Receiver Operator Characteristics for Logistic Regression
#' 
#' Computes the receiver operator characteristics (ROC) for a logistic
#' regression model.
#' 
#' 
#' @param object an object of class "glm."
#' @return An object of class "roc" haveing these components: \item{c.val}{the
#' area under the ROC curve} \item{table}{a 3-column matrix of the fitted
#' values, sensitivity and specificity}
#' @note Objects of class "roc" have \code{print} and \code{plot} methods.
#' @references Gonen, M., 2006, Receiver Operating Characteristics (ROC)
#' Curves: SAS Users Group International Conference, Paper 210-31 18 p.,
#' available at http://www2.sas.com/proceedings/sugi31/210-31.pdf, last
#' accessed Octocer 26, 2011.
#' @keywords regression
#' @export roc
roc <- function(object) {
	# Coding history:
	#    2008Apr03 DLLorenz Original Coding in plotLogistic()
	#    2009Feb11 DLLorenz Moved to this file
	#    2011Aug22 DLLorenz Conversion to R
	#    2011Oct25 DLLorenz Update for package
	#    2012Aug28 DLLorenz Change from diagPlot to plot
	#    2013Apr09 DLLorenz Added setGD to plot
	#    2014May19 DLLorenz fixed ranges
	#    2014Dec29 DLLorenz Convert to roxygen header
  ##
  fits <- fitted(object)
  ## Force to integer, just in case there is some off-1 values
  y <- as.integer(round(object$y,3))
  ## Sort y and fits by fits
  fits.order <- order(fits)
  y.sort <- y[fits.order]
  fits.sort <- fits[fits.order]
  sens <- 1-cumsum(y.sort)/sum(y.sort)
  spec <- cumsum(1-y.sort)/sum(1-y.sort)
  roc.mat <- do.call("rbind", by(cbind(fits=fits.sort, spec=spec, sens=sens),
                               spec, function(x) x[x$sens == max(x$sens),]))
  ## Compute c
  c.tbl <- table(outer(fits[y==1], fits[y==0], function(x,y) sign(x-y)))
  if(length(c.tbl) == 3) # ties in predicted values
    c.val <- (c.tbl[3] + .5*c.tbl[2])/sum(c.tbl)
  else # no ties
    c.val <- c.tbl[2]/sum(c.tbl)
  retval <- list(c.val=c.val, table=roc.mat)
  oldClass(retval) <- "roc"
  return(retval)
}
