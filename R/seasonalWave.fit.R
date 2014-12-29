#' Compute Seasonal Wave Model
#' 
#' This is a support function for seasonalWave model (Vecchia and others,
#' 2008).
#' 
#' 
#' @param cmax the time of the greatest peak value, expressed as a fraction of
#' the year.
#' @param wtx a vector of 12 monthly application rates.
#' @param pkt the numeric month of peak expected from the data in \code{wtx},
#' expressed as a fraction of the year.
#' @param phi the recession rate, in 1/years.
#' @return A vector of 361 values describing the annual seasonal wave.
#' @note This function is support for the seasonalWave function and is not
#' intended to be called by the user.
#' @author Dave Lorenz, original coding by Aldo Vecchia.
#' @seealso \code{\link{seasonalWave}}
#' @references Vecchia, A.V., Martin, J.D., and Gilliom, R.J., 2008, Modeling
#' variability and trends in pesticide concentrations in streams: Journal of
#' the American Water Resources Association, v. 44, no. 5, p. 1308-1324.
#' @keywords manip
seasonalWave.fit <- function(cmax, wtx, pkt, phi) {
	# Coding history:
	#    2007????? SVecchia Original Coding
	#    2007Aug21 DLLorenz Split into 2 functions and added to USGS library
	#    2007Aug22 DLLorenz Modified seasonalWave to take argument x and
	#                       return the wave values
	#    2008Sep10 DLLorenz Updated comments
	#    2009Jul22 DLLorenz Begin modification to allow any number of models
	#    2011May25 DLLorenz Begin Conversion to R and rename
	#    2011Jul06 DLLorenz Prep for package
	#    2012Aug11 DLLorenz Integer fixes
	#    2013Feb02 DLLorenz Modification to more easily facilitate new models
	#    2013Apr02 DLLorenz Final tweaks for release, not exported
	#    2014Dec29 DLLorenz Conversion to roxygen header
	#
  ## This is a seasonalWave support function.
  ## It computes seasonal wave with cmax=time of peak,
  ## wtx=monthly application rates, pkt=month of peak,
  ## phi=recession rate, in 12/months (actually 30-day months)
  del <- 1/12
  txx <- seq(0,1,1/360)
  sxx <- txx
  rho <- exp(-phi)
  z0xx <- rho^(txx)
  r12 <- rho^(-del*c(1:12))
  con <- wtx[1]*(r12[1L]-1)
  for(k in seq(2L,12L)) {con <- con+wtx[k]*(r12[k]-r12[k-1L])}
  con <- rho/(1-rho)*con
  z0xx <- z0xx*con
  zmat <- matrix(nrow=length(txx),ncol=12L)
  pckm <- matrix(nrow=length(txx),ncol=12L)
  ntot <- length(txx)
  pckm[,1] <- c(txx<=del)
  for (k in seq(2L,12L)) {
    pckm[,k] <- c(txx>(k-1)*del & txx<=k*del)
  }
  zmat[,1] <- rho^txx*(rho^(-replace(txx,txx>1/12,1/12))-1)
  for (k in seq(2L,12L)) {
    ztmp <- rep(0,ntot)
    for (j in seq(k-1L)) {
      ztmp[pckm[,j]] <- 0
    }
    ztmp[pckm[,k]] <- 1-rho^(txx[pckm[,k]]-(k-1)*del)
    if(k < 12L){ 
      for (j in seq((k+1), 12L)) {
        ztmp[pckm[,j]] <- rho^(txx[pckm[,j]]-del*k)-rho^(txx[pckm[,j]]-(k-1)*del)}
    }
    zmat[,k] <- ztmp
  }
  sst <- z0xx
  for (k in seq(12)) {
    sst <- sst+wtx[k]*zmat[,k]
  }
  sst <- sst/phi
  medxx <- (max(sst)+min(sst))/2
  rngxx <- 2*(max(sst)-medxx)
  sst <- (sst-medxx)/rngxx
######################
  if(cmax<=pkt) { txx2 <- (txx+1-pkt+cmax)}
  if(cmax>pkt)  {txx2 <- txx-pkt+cmax}
  txx2[txx2>1] <- txx2[txx2>1]-1
  otmp <- order(txx2)
  sst <- sst[otmp]
  ## Set last value to first--originally a dupe of the next-to-last value!
  sst[361L] <- sst[1L]
  sst
}
