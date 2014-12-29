#' Predict Maintenance of Variance Extension, Type 1
#' 
#' Predict new values from a maintenance of variance extension, type 1 (MOVE.1)
#' model fit.
#' 
#' If \code{type} is "response," then the predicted values are
#' back-transformed. Otherwise, the predicted values are computed directly from
#' the model equation.
#' 
#' @param object an object of class "move.1" on which to base the predicted
#' values.
#' @param newdata an optional data.frame in which to look for variables with
#' which to predict.  If omitted, then the fitted values are used.
#' @param type the type of prediction ("response" or "link"). See
#' \bold{Details}.
#' @param \dots not used, required for method function.
#' @return A vector of predictions matching \code{newdata} or the model data.
#' @seealso \code{\link{move.1}}
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
#' @method predict move.1
predict.move.1 <- function(object, newdata, type = c("response", "link")
													 , ...) {
	## Arguments:
	##  object (a move.1 object) the model from which predictions are made
	##  newdata (a data.frame) the new data, if missing use old data
	##  type (character scalar) what to predict
	##  ... (dots) not used, required for method function
	##
	if(missing(newdata))
		newdata <- object$model
	xlab <- attr(attr(object$model, "terms"), "term.labels")
	xindex <- newdata[, xlab] # xindex is simple vector
	ckdist <- FALSE
	if(!is.null(object$call$distribution)) {
		dist <- match.arg(object$call$distribution, c("normal", "lognormal",
																									"commonlog"))
		if(dist == "lognormal") {
			ckdist <- TRUE
			xindex <- log(xindex)
		}
		else if(dist == "commonlog") {
			ckdist <- TRUE
			xindex <- log10(xindex)
		}
	}
	out <- cbind(1, xindex) %*% object$coef
	type <- match.arg(type)
	if(type == "response" && ckdist) {
		if(dist == "commonlog")
			out <- 10^out
		else # Must be natural log
			out <- exp(out)
	}
	as.vector(out) # Strip matrix info
}
