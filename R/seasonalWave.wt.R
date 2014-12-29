#' Compute Seasonal Wave Model
#' 
#' This is a support function for seasonalWave model (Vecchia and others,
#' 2008).
#' 
#' If \code{second.peak} is supplied, then it must be a list with these
#' components:\cr la, the lag (number of months) between the primary and second
#' peak;\cr lo, the loading duration (number of months) for the second peak
#' (must be less than \code{la};\cr w, the weigthing factor for the second peak
#' loading, relative to the primary peak. Must be greater than 0 and less than
#' or equal to 1.
#' 
#' @param loading the number of months of primary loading for the seasonal wave
#' model.
#' @param second.peak a list describing the loading characteristics of the
#' second peak. See \bold{Details}.  If \code{NULL}, then no second peak is
#' added to the seasonal wave model.
#' @return The weighting vector for the seasaon wave model. A vector of length
#' 12 that represents the relative loading to the system.
#' @note This function is support for the seasonalWave function and is not
#' intended to be called by the user.
#' @seealso \code{\link{seasonalWave}}
#' @references Vecchia, A.V., Martin, J.D., and Gilliom, R.J., 2008, Modeling
#' variability and trends in pesticide concentrations in streams: Journal of
#' the American Water Resources Association, v. 44, no. 5, p. 1308-1324.
#' @keywords manip
seasonalWave.wt <- function(loading, second.peak) {
	# Coding history:
	#    2013Feb02 DLLorenz Original coding 
	#    2013Apr02 DLLorenz Final tweaks for release, not exported
	#    2014Dec29 DLLorenz Conversion to roxygen header
	#
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
