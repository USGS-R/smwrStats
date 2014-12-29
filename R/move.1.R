#' Maintenance of Variance Extension, Type 1
#' 
#' Calculates the Maintenance Of Variance Extension, Type 1 (MOVE.1) for record
#' extension by fitting a Line of Organic Correlation (LOC) regression model.
#' 
#' If \code{distribution} is "normal," then the data in \code{x} and \code{y}
#' are assumed to have a bivariate normal distribution. Otherwise, they are
#' assumed to have a bivariate log-normal distribution and a logarithmic
#' transform is applied to both \code{x} and \code{y} before analysis. The
#' natural logarithm is used if \code{distribution} is "lognormal" and the
#' commmon logarithm is used if \code{distribution} is "commonlog."
#' 
#' @param formula a formula object with the response variable on the left of
#' the ~ operator and a single explantory variable on the right.
#' @param data the data frame containing the variables named in \code{formula}.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action a function that indicates what should happen when the data
#' contain missing values (NAs).  The default is set by the \code{na.action}
#' setting of \code{options} and is \code{na.fail} if that is unset. Ohter
#' possible options include \code{na.exclude} and \code{na.omit}.
#' @param distribution either "normal," "lognormal," or "commonlog" indicating
#' nature of the bivariate distribution, See \bold{Details}.
#' @return An object of class "move.1" having these components:
#' \item{coefficients}{ the intercept and slope of the line describing the fit.
#' } \item{na.action}{ a character string indicating the special handling of
#' NAs. } \item{R}{ Pearson's correlation coefficient. } \item{call}{ the
#' matched call to \code{move.1}. } \item{fitted.values}{ the fitted LOC values
#' for the response. } \item{residuals}{ a 2-column matrix containing the
#' signed distance from the predicted to to the corresponding \code{y} and
#' \code{x} values. } \item{x}{ the (possibly transformed) values for \code{x}.
#' } \item{y}{ the (possibly transformed) values for \code{y}. }
#' \item{x.stats}{ the mean and standard deviation of the (possibly
#' transformed) values for \code{x}. } \item{y.stats}{ the mean and standard
#' deviation of the (possibly transformed) values for \code{y}. }
#' \item{var.names}{ the names of \code{y} and \code{x}. } \item{model}{ the
#' model frame. }
#' @note Objects of class "move.1" have \code{print}, \code{predict}, and
#' \code{plot} methods.
#' @seealso \code{\link{predict.move.1}}, \code{\link{plot.move.1}}
#' @references Hirsch, R.M., 1982, A comparison of four streamflow record
#' extension techniques: Water Resources Research, v. 18, p. 1081--1088.
#' @keywords models regression
#' @examples
#' 
#' library(smwrData)
#' data(IonBalance)
#' # Build model for non missing Alkalinity
#' IB.move <- move.1(Anion_sum ~ Cation_sum, data=IonBalance, subset=abs(Pct_Diff) < 10) 
#' print(IB.move)
#' 
#' @export move.1
move.1 <- function(formula, data, subset, na.action, distribution="normal") {
	# Coding history:
	#    2006Dec05 DLLorenz Initial version
	#    2009Mar19 DLLorenz Dated version
	#    2010Dec08 DLLorenz Conversion to R
	#    2011Aug05 DLLorenz Added print and plot (diagPlot) files
	#    2011Oct25 DLLorenz Update for package
	#    2012Aug28 DLLorenz Change from diagPlot to plot
	#    2012Nov21 DLLorenz Added Q-Q plot to diagnostics
	#    2012Dec21 DLLorenz Change gaussian to normal
	#    2013Apr09 DLLorenz Added setGD to plot
	#    2014Dec22 DLLorenz Roxygen headers
  ##
  ## Process formula as with regular linear regression!
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$distribution <- m$modeluse.log10 <- NULL
  m$drop.unused.levels <- TRUE
  m[[1L]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  if(ncol(m) != 2L)
    stop("Must specify a single response and a single explanatory variable in the formula")
  Y <- m[, 1L]
  yname <- names(m)[1L]
  X <- m[, 2L]
  xname <- names(m)[2L]
  fit <- list(coefficients=double(2L))
  ## Extract missing value info
  fit$na.action <- attr(m, "na.action")
  ## Adjust for transforms
  distribution <- match.arg(distribution, c("normal", "lognormal", "commonlog"))
  if(distribution == "lognormal") {
    Y <- log(Y)
    X <- log(X)
    xname <- paste("log(", xname, ")", sep='')
    yname <- paste("log(", yname, ")", sep='')
  }
  else if(distribution == "commonlog") {
    Y <- log10(Y)
    X <- log10(X)
    xname <- paste("log10(", xname, ")", sep='')
    yname <- paste("log10(", yname, ")", sep='')
  }
  ybar <- mean(Y)
  yvar <- var(Y)
  xbar <- mean(X)
  xvar <- var(X)
  retcor <- cor.test(X, Y)
  fit$R <- retcor$estimate
  fit$p.value <- retcor$p.value
  fit$coefficients[2L] <- sqrt(yvar/xvar)
  if(fit$R < 0)
    fit$coefficients[2L] <-  - fit$coefficients[2L]
  fit$coefficients[1L] <- ybar - fit$coefficients[2L] * xbar
  names(fit$coefficients) <- c("(Intercept)", xname)
  ## Complete steps for the model
  fit$call <- call
  fit$fitted.values <- cbind(1, X) %*% as.matrix(fit$coefficients)
  yresiduals <- Y - fit$fitted.values
  ## Now reverse the sense of the fit to get X-residuals
  junk <- double(2L)
  junk[2L] <- 1/fit$coefficients[2L]
  junk[1L] <- xbar - junk[2L] * ybar
  xfitted <- cbind(1, Y) %*% as.matrix(junk)
  xresiduals <- X - xfitted
  fit$residuals <- cbind(yresiduals, xresiduals)
  colnames(fit$residuals) <- c(yname, xname)
  ## Pack the object with other necessary information
  fit$x <- X
  fit$y <- Y
  fit$xstats <- c(mean=xbar, sd=sqrt(xvar))
  fit$ystats <- c(mean=ybar, sd=sqrt(yvar))
  fit$var.names <- c(yname, xname)
  fit$model <- m
  ## set class to move.1
  oldClass(fit) <- "move.1"
  return(fit)
}

## The default extraction methods: coef, residuals, fitted work
