# Function to compute seasonal wave model
# Note! Whenever the model definitions change, confirm.seasonalPeak must
#   be updated
#
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
#    2012Aug11          This version
#

seasonalWave <- function(x, cmax, loading, hlife, second.peak=NULL) {
  ## Arguments:
  ##  x (numeric vector) the dates in decimal year format
  ##  cmax (numeric scalar) the fraction part of the year at which the
  ##    peak occurs.
  ##  loading (numeric scalar) the number of months of loading
  ##  hlife (numeric scalar) the half life (decay rate), must be 1-4
  ##  second.peak (NULL or list) indicating the parameters of the secondary
  ##    peak
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
