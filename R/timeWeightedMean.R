#' Time-Weighted Mean
#' 
#' Compute the mean weighted by the duration between values.
#' 
#' The values of \code{time} are expected to be in decimal format, where the
#' integer part indicates the period and the fractional part linearly
#' distributed through the period. For example, the year 2000 begins at 2,000.0
#' and July 2, 2000 is 2,000.5. The function \code{dectime} can be used to
#' convert \code{Date} or \code{POSIX} to decimal format. If \code{time} is of
#' class "Date," or any "POSIX" classs, it is converted using
#' \code{dectime}.\cr If \code{na.rm} is \code{FALSE} and \code{byeriod} is
#' \code{TRUE}, then missing values are returned only for periods that contain
#' a missing value (\code{NA}).
#' 
#' @param x a sequnece of numeric values whose mean is to be computed.
#' @param time the times associated with \code{x}. See \bold{Details}.
#' @param na.rm remove missing values before computing the mean?
#' @param fill expand the weights for the first and last observation to the end
#' of the period? Otherwise, the weight is based on the distance to the second
#' or next-to-last observation.
#' @param by.period computes wights by the periods defined by \code{time}.
#' @param excessive a warning is issued if the largest weight for any
#' observation exceeds \code{excessive} times the average weight.
#' @return If \code{by.period} is \code{TRUE}, then a vector with one entry per
#' period, otherwise the mean for the entire data set.
#' @seealso \code{\link{weighted.mean}}, \code{\link{dectime}}
#' @references Crawford,C.G., 2004, Sampling strategies for estimating acute
#' and chronic exposures of pesticides in streams: Journal of the American
#' Water Resources Association, v. 40, n. 2, p 485-502.
#' @keywords univar
#' @examples
#' \dontrun{
#' library(smwrData)
#' data(QW05078470)
#' with(QW05078470, timeWeightedMean(P00665, dectime(DATES)))
#' }
#' @export timeWeightedMean
timeWeightedMean <- function(x, time, na.rm=TRUE, fill=TRUE,
                           by.period=FALSE, excessive=2.5) {
	# Coding history:
	#    2010Aug17 DLLorenz Original
	#    2012Jun04 DLLorenz Conversion to R
	#    2014Dec29 DLLorenz Conversion to roxygen header
  ##
  ## time must be expressed in decimal time, force to be compliant
  time <- dectime(time)
  if(na.rm) {
    sel <- !(is.na(x) | is.na(time))
    x <- x[sel]
    time <- time[sel]
  }
  if(any(diff(time) < 0))
    stop('time must be sorted in increasing order')
  N <- length(x)
  if(N < 2)
    stop('requires more than a single observation')
  tmin <- floor(min(time))
  tmax <- ceiling(max(time))
  EMW <- (tmax - tmin)/N # expected mean weight
  if(by.period && !fill) { # force fill to be true and warn
    fill <- TRUE
    warning('fill forced to TRUE for by.period = TRUE')
  }
  if(fill) { # fill both ends to full period
    tmin <- 2*tmin - time[1] 
    tmax <- 2*tmax - time[N]
  }
  else { # extend ends to symmetric periods
    tmin <- 2*time[1] - time[2]
    tmax <- 2*time[N] - time[N-1]
  }
  cuttimes <- (c(tmin, time) + c(time, tmax))/2
  weights <- diff(cuttimes)
  ## Warn if excessive weight
  if((max(weights) / EMW) > excessive)
    warning('at least one weight is excessive (', excessive,
            ' * larger than expected)')
  if(by.period) {
    Seq <- unique(as.integer(time))
    retval <- sapply(Seq, FUN=function(i, x, ct) {
      ct <- pmin(i+1,pmax(ct,i))
      w <- diff(ct)
      weighted.mean(x, w=w)}, x=x, ct=cuttimes)
    names(retval) <- as.character(Seq) # Needed becuase Seq is integer
    return(retval)
  }
  ## otherwise just return the mean
  return(weighted.mean(x, w=weights))
}
