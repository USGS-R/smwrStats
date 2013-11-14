# Compute the (monthly) weighting vector
#
# Coding history:
#    2013Feb02 DLLorenz Original coding 
#    2013Apr02 DLLorenz Final tweaks for release, not exported

seasonalWave.wt <- function(loading, second.peak) {
  ## This is a seasonalWave support function.
  ## It computes the wtx (monthly application rates) vector.
  wtx <- c(rep(1, loading), rep(0, 12L - loading))
  if(!is.null(second.peak)) { # Add the second peak info
    lag <- second.peak$la  # Lag in months from peak
    if(is.null(lag))
      stop("second.peak must contain the \"lag\" component")
    if(lag + loading > 11)
      stop("second.peak lag plus loading must be less than 12 months")
    load <- second.peak$lo # Loading duration
    if(is.null(lag))
      stop("second.peak must contain the \"loading\" component")
    if(lag - load < 1)
      stop("second.peak loading must be shorter than the lag")
    wt <- second.peak$w    # weigth (0-1) scale
    if(is.null(wt))
      stop("second.peak must contain the \"weight\" component")
    if(any(wt > 1) || any(wt[1] <= 0))
      stop("second.peak weight must be between 0 and 1")
    ndx <- seq(load) + loading + lag - load
    wtx[ndx] <- wt
  }
  return(wtx)
}
