# calculate the time-weighted mean
#
# Coding history:
#    2010Aug17 DLLorenz Original
#    2012Jun04 DLLorenz Conversion to R
#    2012Jun05          This version.
#

timeWeightedMean <- function(x, time, na.rm=TRUE, fill=TRUE,
                           by.period=FALSE, excessive=2.5) {
  ## Arguments:
  ##  x, the variable to compute the mean
  ##  time, the date/time of the observation
  ##  na.rm, remove missings/
  ##  fill, extend the range to limits of the period
  ##  by.period, summarize by period (whole periods of time)
  ##  excessive, the limit for excessive weighting of any observation
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
