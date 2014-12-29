#' Predict Values.
#' 
#' Predict new values from Sen slope (senSlope) model fit.
#' 
#' 
#' @param object an object of class "senSlope" on which to base the predicted
#' values.
#' @param newdata an optional data.frame in which to look for variables with
#' which to predict.  If omitted, then the fitted values are used.
#' @param \dots not used, required for method function.
#' @return A vector of predictions matching \code{newdata} or the model data.
#' @seealso \code{\link{senSlope}}
#' @keywords models regression
#' @export
#' @method predict senSlope
predict.senSlope <- function(object, newdata, ...) {
  ##
  if(missing(newdata))
    newdata <- object$model
  xlab <- attr(attr(object$model, "terms"), "term.labels")
  xindex <- newdata[, xlab] # xindex is simple vector
  ckdist <- FALSE
  out <- cbind(1, xindex) %*% object$coefficients
  as.vector(out) # Strip matrix info
}
