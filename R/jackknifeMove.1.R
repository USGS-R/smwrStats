#' Jackknife move.1
#' 
#' Computes the bias and variance of move.1 predictions using the jackknife estimator.
#'
#' If \code{type} is "response," then the predicted values are
#' back-transformed. Otherwise, the predicted values are computed directly from
#' the model equation.
#' 
#' @param object an object of class "move.1" on which to base the predicted
#' values.
#' @param newdata an optional data.frame in which to look for variables with
#' which to predict.  If omitted, then the calibration data are used; the response
#' values are sorted from smallest to largest and \code{type} is set to "link" if it
#' is not set to "response" in the call.
#' @param type the type of prediction ("response" or "link"). See
#' \bold{Details}.
#' @return A vector of predictions matching \code{newdata} or the model data.
#' @seealso \code{\link{move.1}}, \code{\link{predict.move.1}}
#' @references
#'Lorenz, D.L., 2015, smwrStats-an R package for analyzing hydrologic data, 
#'version 0.7.0: U.S. Geological Survey Open-File Report 2015-XXXX, XX p.
#' @keywords models regression
#' @examples
#' 
#' library(smwrData)
#' data(IonBalance)
#' # Build model for non missing Alkalinity
#' IB.move <- move.1(Anion_sum ~ Cation_sum, data=IonBalance, subset=abs(Pct_Diff) < 10) 
#' print(IB.move)
#' # Predict Anion_sum for missing Alkalinity
#' predict(IB.move, IonBalance[1, ])
#' 
#' @export
jackknifeMove.1 <- function(object, newdata, type = c("response", "link")) {
	if(missing(newdata)) {
		# construct the original data in order, sorted by x to predict sorted y
		newdata <- object$model
		newdata <- newdata[order(newdata[, 2]),]
		if(missing(type)) {
			type <- "link"
		}
	}
	est <- predict(object=object, newdata=newdata, type=type)
	N <- length(object$x)
	out <- matrix(0, ncol=N, nrow=length(est))
	# Make the jackknife estimate
	for(i in seq(N)) {
		out[,i] <- predict(update(object, subset=-i), newdata=newdata, type=type)
	}
	# Replace out with the scaled differences
	out <- (N-1)*(est - out)
	# Compute the bias and variance of the estimate (mean)
	bias <- -1/N * rowSums(out)
	var <- 1/(N*(N-1)) * (rowSums(out^2) - N * bias^2)
	retval <- cbind(Est=est, Bias=bias, Var=var)
	return(retval)
	
}
													 