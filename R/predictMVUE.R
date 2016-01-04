#' Bias Corrected Predictions
#' 
#' Predicts bias-corrected expected mean response values from a log-transformed
#' regression model, using either the minimum variance unbiased estimate(MVUE),
#' Duan's smoothing estimate, or Ferguson's maximum likelihood estimate.
#' 
#' 
#' @aliases predictMVUE predictDuan predictFerguson
#' @param object an object of class "lm" on which to base the predicted values.
#' @param newdata an optional data.frame in which to look for variables with
#' which to predict.  If omitted, then the fitted values are used.
#' @param back.trans the back-transformation function. For common log
#' transforms, use \code{function(x) 10^x}.
#' @param Log10 is the transform of the response variable the common log?
#' @return A vector of predictions matching \code{newdata} or the model data.
#' @seealso \code{\link{lm}}
#' @references Bradu, D. and Mundlak, Y., 1970, Estimation in the lognormal
#' linear models: Journal of the American Statistical Association, v. 65, no.
#' 329, p. 198--211.
#' 
#' Duan, N., 1983, Smearing estimate: a nonparametric retransformation method:
#' Journal of the American Statistical Association, v. 78, p. 159--178.
#' 
#' Ferguson, R.I. 1986, River loads underestimated by rating curves: Water
#' Resources Research, v. 22, p 74--76.
#' 
#' Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in water resources:
#' U.S. Geological Survey Techniques of Water-Resources Investigations, book 4,
#' chap. A3, 522 p.
#' @keywords models regression
#' @examples
#' 
#' ## Generate random log-normal data and build the regression model
#' set.seed(111)
#' XX.df <- data.frame(x=sort(runif(32, 1, 5)), y=rlnorm(32, seq(1,2, length.out=32)))
#' XX.lm <- lm(log(y) ~ x, data=XX.df)
#' ## Compare the results for x=1:5
#' ## The simple back-transformed estimates
#' exp(predict(XX.lm, newdata=data.frame(x=1:5)))
#' ## The bias corrected estimates of the mean response
#' predictFerguson(XX.lm, newdata=data.frame(x=1:5))
#' predictDuan(XX.lm, newdata=data.frame(x=1:5))
#' predictMVUE(XX.lm, newdata=data.frame(x=1:5))
#' 
#' @useDynLib smwrStats phimvue
#' @export predictMVUE
predictMVUE <- function(object, newdata, Log10=FALSE) {
	# Coding history:
	#    2006Dec05 DLLorenz Original coding for Pat Rasmussen
	#    2007Jan26 DLLorenz Added to USGS library and added options
	#    2007Feb13 DLLorenz Made confidence interval future option
	#    2007May22 DLLorenz Added confidence interval
	#    2008Jan08 DLLorenz Renamed to predictMVUE to avoid confusion with methods
	#    2012Apr05 DLLorenz Conversion to R
	#    2012Jul31 DLLorenz Removed options related to ci.fit (see USGSqw for CI)
	#    2012Oct12 DLLorenz Bug fix in call to phimvue
	#    2014Dec29 DLLorenz Conversion to roxygen headers
	#
  Fact <- if(Log10) 2.30258509299405 else 1.
  if(missing(newdata))
    newdata <- eval(as.list(object$call)$data)
  firstguess <- predict(object, newdata, type = "response", se.fit = TRUE)
  firstguess$AW <- ((1 - (firstguess$se.fit/firstguess$residual.scale)^2) *
                    firstguess$df * (firstguess$residual.scale * Fact)^2)
  ## Compute the bias correction factor
  N <- length(firstguess$fit)
  firstguess$BCF <- .Fortran("phimvue",
                             as.double(firstguess$AW/4),
                             as.double(firstguess$df/2),
                             x = double(N),
                             Num = as.integer(N))$x
  retval <- exp(firstguess$fit * Fact) * firstguess$BCF
  return(retval)
}
