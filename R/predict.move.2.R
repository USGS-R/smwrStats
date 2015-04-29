#' Predict Maintenance of Variance Extension, Type 2
#' 
#' Predict new values from a maintenance of variance extension, type 2 (MOVE.2)
#' model fit.
#' 
#' @details If \code{type} is "response," then the predicted values are
#' back-transformed. Otherwise, the predicted values are computed directly from
#' the model equation.
#' 
#' @param object an object of class "move.2" on which to base the predicted
#' values.
#' @param newdata an optional data.frame in which to look for variables with
#' which to predict.  If omitted, then the fitted values are used.
#' @param type the type of prediction ("response" or "link"). See
#' \bold{Details}.
#' @param \dots not used, required for method function.
#' @note If lag was set to a non-zero value in the call to \code{move.2}, then the
#'explanatory variable is lagged only when predictiong values from the calibration
#'data (\code{newdata} is not supplied.) This facilitates prediction of selected
#'statistics at the response site rather than the complete record.
#' @return A vector of predictions matching \code{newdata} or the model data.
#' @seealso \code{\link{move.2}}
#' @keywords models regression
#' @examples
#'\dontrun{
#'# See the vignette:
#'vignette("RecordExtension", package="smwrStats")
#'}
#' @export
#' @method predict move.2
predict.move.2 <- function(object, newdata, type = c("response", "link"), ...) {
  ## Set lag to 0 and reset only for prediction from original data
  lag <- 0
  if(missing(newdata)) {
    newdata <- object$data
    lag <- object$lag
  }
  type <- match.arg(type)
  xlab <- attr(attr(object$model, "terms"), "term.labels")
  xindex <- shiftData(newdata[, xlab], lag) # xindex is simple vector
  ckdist <- FALSE
  if(inherits(object$distribution, "character")) {
    dist <- object$distribution
    if(dist == "lognormal") {
      ckdist <- TRUE
      xindex <- log(xindex)
    } else if(dist == "commonlog") {
      ckdist <- TRUE
      xindex <- log10(xindex)
    } else if(dist == "log1p") {
    	ckdist <- TRUE
    	xindex <- log1p(xindex)
    }
    out <- cbind(1, xindex) %*% object$coef
    type <- match.arg(type)
    if(type == "response" && ckdist) {
      if(dist == "commonlog") {
        out <- 10^out
      } else if(dist == "log1p") {
      	out <- pmax(expm1(out), 0)
      } else # Must be natural log
        out <- exp(out)
    }
  } else { # distribution is from boxCox
    ## Convert predictor
    xindex <- boxCox(xindex, object$xtrans, object$xgm)
    out <- cbind(1, xindex) %*% object$coef
    if(type == "response") # Back transform
      out <- IboxCox(out, object$ytrans, object$ygm)
  }
  as.vector(out) # Strip matrix info
}
