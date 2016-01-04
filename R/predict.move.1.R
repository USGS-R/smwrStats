#' Predict Maintenance of Variance Extension, Type 1
#' 
#' Predicts new values from a maintenance of variance extension, type 1 (MOVE.1)
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
#' @param var.fit logical if \code{TRUE}, then compute the variance of the predicted
#' values. If \code{FALSE}, then the varainces are not computed.
#' @param \dots not used, required for method function.
#' @return If \code{var.fit} is \code{FALSE}, then a vector of predictions matching 
#'\code{newdata} or the model data. If \code{var.fit} is \code{TRUE}, then a data frame
#'containing the columns:
#'\item{fit}{the predicted values}
#'\item{var.fit}{the variance of the predicted values}
#'
#' @seealso \code{\link{move.1}}, \code{\link{jackknifeMove.1}}
#' @keywords models regression
#' @references
#'Lorenz, D.L., 2015, smwrStats-an R package for analyzing hydrologic data, 
#'version 0.7.0: U.S. Geological Survey Open-File Report 2015-XXXX, XX p.
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
predict.move.1 <- function(object, newdata, type = c("response", "link"),
													 var.fit=FALSE, ...) {
	##
	## Function to baack-compute the variance
	fixVar <- function(mu, sg2) {
		emu <- exp(mu + sg2/2)
		esg2 <- (exp(sg2) - 1)*exp(2*mu + sg2)
		retval <- esg2 + (exp(mu) - emu)^2
		return(retval)
	}
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
	if(var.fit) {
		Ns <- length(object$x)
		var <- (1-object$R^2)*object$ystats[2L]^2*(Ns-1)/(Ns-2)*
			(1/Ns + (xindex - object$xstats[1L])^2/(object$xstats[2L]^2*(Ns-1)))
	}
	type <- match.arg(type)
	if(type == "response" && ckdist) {
		if(dist == "commonlog") {
			if(var.fit) {
				var <- fixVar(log(10)*out, log(10)^2*var)
			}
			out <- 10^out
		} else { # Must be natural log
			if(var.fit) {
				var <- fixVar(out, var)
			}
			out <- exp(out)
		}
	}
	if(var.fit) {
		return(data.frame(fit=out, var.fit=var))
	} else {
		return(as.vector(out)) # Strip matrix info
	}
}
