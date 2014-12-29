#' Seasonal Wave
#' 
#' Compute a seasonal wave model to describe the variation of a constituent
#' over the course of a year. This model is particularly well suited to
#' describe the concentration of pesticides.
#' 
#' The seasonal wave model \eqn{W(t)} is expressed as a differential equation
#' \deqn{\frac{d}{dt} W(t) = \lambda(t) - \phi W(t) [0 \le t \le 1, }{ dW(t)/dt
#' = \lambda(t) - \phi W(t) [0 \le t \le 1, W(0)=W(1)]}\deqn{ W(0)=W(1)]}{
#' dW(t)/dt = \lambda(t) - \phi W(t) [0 \le t \le 1, W(0)=W(1)]}
#' \deqn{\lambda(t) = \sum_{k=1}^12 \omega_k I(\frac{k-1}{12} \le t < }{
#' \lambda(t) = \sum \omega_k I((k-1)/12 \le t < k/12)}\deqn{ \frac{k}{12})}{
#' \lambda(t) = \sum \omega_k I((k-1)/12 \le t < k/12)} where \eqn{\phi} is a
#' specified constant, which represents the rate of removal or chemical decay;
#' \eqn{[\omega_k, k=1,2,...,12]} are specified nonnegative constants, which
#' represent monthly input amounts; \eqn{I(.)} is the indicator function with
#' \eqn{I(.)=1} if \eqn{t} lies in the given interval and \eqn{I(.)=0}
#' otherwise; and \eqn{\lambda(t)} is the input amount at time \eqn{t}.\cr
#' 
#' The value of hlife is used to define the decay rate \eqn{\phi} in the
#' differential equation. The value of \eqn{\phi} is 12/\code{hlife} and
#' \eqn{W(t)} decays at a a rate of \eqn{exp(-\phi/12)} per month.\cr
#' 
#' The list for \code{second.peak} must contain these components:\cr la, the
#' lag from the primary peak to the second peak, in months. This is always the
#' time from the primary peak to the second peak, even if the second peak
#' occurs earlier in the year.\cr lo, the loading duration in months, this
#' value must be leas than \code{la}. w, the load scaling factor relative to
#' the primary peak. Must be greater than 0 and less than or equal to 1. For
#' practical pruposes, it should be greater than 0.5.
#' 
#' @param x a vector of decimal time representing dates. Missing values
#' (\code{NA}s) are permitted.
#' @param cmax the time of the greatest peak value, expressed as a fraction of
#' the year.
#' @param loading the number of months of contituent loading for the primary
#' peak.
#' @param hlife the half life of the decay rate, expressed in units of months.
#' This should be in the range of 1 though 4.
#' @param second.peak a list of the parameters for the second peak. See
#' \bold{Details}.
#' @return A vector expressing the expected variation of concentration for each
#' value in x; the values are scaled to a range of -0.5 to 0.5.
#' @author Dave Lorenz, original coding by Aldo Vecchia.
#' @seealso \code{\link{seasonalPeak}}, \code{\link{confirm.seasonalPeak}}
#' @references Vecchia, A.V., Martin, J.D., and Gilliom, R.J., 2008, Modeling
#' variability and trends in pesticide concentrations in streams: Journal of
#' the American Water Resources Association, v. 44, no. 5, p. 1308-1324.
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' # Selected single peak models
#' # 1 month loading, 1 month half-life
#' curve(seasonalWave(x, 3/12, 1, 1), 0, 1, n=361, xlab='fraction of year', ylab="W")
#' # 1 month loading, 3 month half-life
#' curve(seasonalWave(x, 3/12, 1, 3), 0, 1, n=361, add=TRUE, col="blue")
#' # 3 month loading, 2 month half-life
#' curve(seasonalWave(x, 3/12, 3, 2), 0, 1, n=361, add=TRUE, col="green")
#' # Add a second peak model
#' curve(seasonalWave(x, 3/12, 3, 2, second.peak=list(la=6, lo=2, w=.75) ), 0, 1, 
#'   n=361, add=TRUE, col="red")
#' }
#' 
#' @export seasonalWave
seasonalWave <- function(x, cmax, loading, hlife, second.peak=NULL) {
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
	#    2014Dec29 DLLorenz Conversion to roxygen header
  ##
  ## This is the user interface
  ##
  phi <- 12/hlife # compute directly rather than look up
  ## Construct the weighting (wtx) and peak timing of the loading model
  loading <- min(as.integer(loading), 9L)
  pkt <- loading/12
  ## The .wtx function computes the loading vector
  wtx <- seasonalWave.wt(loading, second.peak)
  ## The .fit function does the actual computation
  awave <- seasonalWave.fit(cmax, wtx, pkt, phi)
  x <- x - floor(x) # extract the decimal part
  ## Interpolate the values from the fit.
  return(approx(seq(0, 1, 1/360), awave, xout=x)$y)
}
