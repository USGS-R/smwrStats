# correct log-transformed response variable using MVUE
#
# Coding history:
#    2006Dec05 DLLorenz Original coding for Pat Rasmussen
#    2007Jan26 DLLorenz Added to USGS library and added options
#    2007Feb13 DLLorenz Made confidence interval future option
#    2007May22 DLLorenz Added confidence interval
#    2008Jan08 DLLorenz Renamed to predictMVUE to avoid confusion with methods
#    2012Apr05 DLLorenz Conversion to R
#    2012Jul31 DLLorenz Removed options related to ci.fit (see USGSqw for CI)
#    2012Oct12 DLLorenz Bug fix in call to phimvue
#    2013Oct12          This version.

predictMVUE <- function(object, newdata, Log10=FALSE) {
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
