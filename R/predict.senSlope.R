predict.senSlope <- function(object, newdata, ...) {
  ## Arguments:
  ##  object (a senSlope object) the model from which predictions are made
  ##  newdata (a data.frame) the new data, if missing use old data
  ##  ... (dots) not used, required for method function
  ##
  if(missing(newdata))
    newdata <- object$model
  xlab <- attr(attr(object$model, "terms"), "term.labels")
  xindex <- newdata[, xlab] # xindex is simple vector
  ckdist <- FALSE
  out <- cbind(1, xindex) %*% object$coefficients
  as.vector(out) # Strip matrix info
}