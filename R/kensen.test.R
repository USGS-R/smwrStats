# Kendall's tau used as a trend estimator (the t vector is time) along
#    with the Sen slope estimator.
#
# Coding History:
#    2000Oct27 JRSlack  Original coding.
#    2001Feb01 JRSlack  Additional input checking.
#    2003Apr04 JRSlack  Change variable x to y plus editorial cleanup.
#    2005Mar24 DLLorenz Fixed removal of NAs
#    2005May04 DLLorenz made class htest
#    2005Jul14 DLLorenz Bug fix
#    2006Apr10 DLLorenz Modifed to print slope and medians
#    2007Jul12 DLLorenz Removes 1:1 match between values and time
#    2007Aug27 DLLorenz Added S and VarS to return value
#    2007Sep05 DLLorenz Added protection against all ties in data
#    2008Sep17 DLLorenz Bugfix to computation on n
#    2011May02 DLLorenz Conversion to R (major changes needed)
#    2011Oct26 DLLorenz Renamed to indicate that a hypo test is done
#    2013Mar26 DLLorenz Added serial correction 
#

kensen.test <- function(y, t, n.min=10) {
  ## Arguments:
  ##  y (numeric vector) the observations for trend analysis
  ##  t (numeric vector) the time assocated with y
  ## Set names
  data.name <- paste(deparse(substitute(y)), "and", deparse(substitute(t)))
  ## Remove missing values and error checking.
  ## Make sure every data value has a time value.
  sel <- !is.na(y)
  y <- y[sel]
  t <- t[sel]
  if (any(is.na(t)))
    stop("Each data value must have a time value.")
  ## Make sure time is strictly increasing
  dift <- diff(t)
  if (any(dift <= 0))
    stop("Time vector must be strictly increasing.")
  ## Make sure we have enough data.
  n <- length(y)
  if (n < 3)
    stop("y and t should effectively be longer than 2.")  
  ## Calculate the statistics
  if(diff(range(y)) == 0) { # Protect against error message from cor.test
    retval <- list(statistic=0, parameter=NULL, p.value=1,
      estimate=c(tau=0), null.value=c(tau=0), 
      alternative="two.sided", method="", data.name="")
    class(retval) <- "htest"
    n.min <- Inf # Suppress adjustment
  }
  else
   retval <- cor.test(y, t, method='kendall', continuity=TRUE, exact=FALSE)
  retval$method <- "Kendall's tau with the Sen slope estimator"
  ## a fast, efficient way to compute slopes:
  slopes <- unlist(lapply(seq(along = y), function(i, y, t)
                                    ((y[i] - y[1:i]) / (t[i] - t[1:i])), y, t))
  ## The Sen slope estimator.
  sen.slope <- median(slopes,na.rm=TRUE)
  ## Median of the data values.
  median.data <- median(y)         
  ## Median of the time values.
  median.time <- median(t)
  ## A line representing the trend of the data then is given by
  ##
  ##    y(t) = sen.slope*(t-med.time)+med.data
  ##
  ##    This line has the slope of the trend and passes through
  ##       the point (t=med.time, y=med.data)
  ## Compute the coefficients for the line
  coef <- c(sen.slope*(-median.time)+median.data, sen.slope) # intercept and slope
  est <- c(sen.slope, median.data, median.time, n)
  names(est) <- c("slope", "median.data", "median.time", "n")
  ## Check if ajdustment can/should be made
  if(n.min < n && n > 9 && max(abs(dift - 1)) < 0.02) { # A small wiggle room
    Resids <- (y - coef[2L] * t)
    cor1 <- acf(Resids, lag.max=1L, plot=FALSE)$acf[-1L, 1L, 1L] # drop 0 lag
    ## Use eqn (7) Matalas & Langbein
    adjS <- 1 + 2 *(cor1^(n+1)-n*cor1^2+(n-1)*cor1)/(n*(cor1-1)^2)
    nstar <- n/adjS
    retval$p.value <- (1 - pnorm(abs(retval$statistic/sqrt(adjS))))*2
    retval$method <- "Kendall's tau (adjusted p-value) with the Sen slope estimator"
    est <- c(est, "n*"=nstar)
  }
  ## Return the statistics.
  retval$statistic <- retval$estimate # Replace T/z with tau
  retval$estimate <- est 
  retval$coef <- coef
  retval$data.name <- data.name
  return(retval)
}

