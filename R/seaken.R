# Seasonal Kendall test for trend with serial correlation correction.
#
# Usage:
#    seaken(series,nseas)
#       where: series is the equally spaced seasonal time series
#                 to examine for trend.
#              nseas is the number of seasons in the timeseries. This
#                 argument is optional and is assumed to 12 if absent.
#
# Limitations:
#    The number of time periods per year (nseas) must not be greater
#       then 52;i.e., the season can not be any shorter than one week.
#    The number of years may not be greater than 50.
#    Thus n must be no greater than 2600 (=52*50).
#
# References:
# 1) A Nonparametric Trend Test for Seasonal Data with Serial
#      Dependence, R.M.Hirsch & J.R.Slack, Water Resources Research,
#      Vol.20, No.6, Pages 727-732, June 1984.
#
# Coding history:
#    2000Oct18 JRSlack  Initial coding.
#    2000Dec22 JRSlack  Simplify the returned list and maintain mode single.
#    2002Feb20 JRSlack  Modified for S-PLUS 6 conventions.
#    2002Dec19 JRSlack  Remove upper limits on seasons and years.
#    2005Jul14 DLLorenz Date fix
#    2006Apr10 DLLorenz Modified to ouput htest and print slope and medians
#    2006Apr11 DLLorenz Added years seasons to data.name.
#    2006May26 DLLorenz Bug fix in retval
#    2012Feb09 DLLorenz Conversion to R
#    2012Aug06 DLLorenz Added series to return value and added "seaken" class
#    2012Dec17 DLLorenz Bug fix in call to set up axis
#    2012Dec17          This version.

seaken <- function(series, nseas=12) {
# Do error checking before calling seaken.
  Dname <- deparse(substitute(series))
  x <- as.single(series)
  n <- as.integer(length(series))
  ns <- as.integer(nseas)
  if (n/ns < 2) stop ("seaken requires at least 2 years of data.")
# Pad any trailing partial year to a full year. The original Fortran
#    code truncated the series using n<-floor(n/ns)*ns
  if (n%%ns != 0) {
    nfull <- ceiling(n/ns)*ns
    npo <- n+1
    x[npo:nfull] <- as.single(-99999.0)
    warning(paste("The original series of", n,
                  "values was padded with missing values to a length of", nfull,
                  "so that it contains an integral number of years."))
    n <- as.integer(nfull)
  }
  
  x[is.na(x)] <- as.single(-99999.0)   # Convert any NAs to -99999.0

# Call the Fortran seaken DLL
# x is the equally spaced seasonal time series to examine for trend
# n is the length of the timeseries
# ns is the number of seasons in the timeseries
# results is the statistics of the trend
  results <- .Fortran("seakenf", x , n , ns , as.single(vector("numeric", 9)))[[4]]
  ## Return the statistics.
  method <- "Seasonal Kendall with correlation correction"
  tau <- results[1]
  names(tau) <- "tau"
  zero <- 0
  names(zero) <- "slope"
  est <- c(results[4], results[6],  results[5])
  names(est) <- c("slope", "median.data", "median.time")
  ## use appropriate p.value
  if(n/ns < 10)
    p.value <- results[2]
  else
    p.value <- results[3]
  z <- list(method = method, data.name =
            paste(Dname, " (", n/ns, " years and ", ns, " seasons)", sep=""),
            nyears = n/ns, nseasons = ns, series=series,
            statistic=tau, p.value=p.value,
            p.value.raw = results[2], p.value.corrected = results[3],
            estimate=est, alternative = "two.sided", null.value = zero)
  oldClass(z) <- c("htest", "seaken")
  return(z)
}

seriesPlot.seaken <- function(x, # data
                              SeasonLine=list(name="", what="vertical", color="black"),
                              SeasonPoint=list(name="", what="points", symbol="circle", 
                                filled=TRUE, size=0.09, color="black"), # plot controls
                              yaxis.log=FALSE, yaxis.range=c(NA,NA), # y-axis controls
                              ylabels=7, xlabels, # labels
                              xtitle="",
                              ytitle=deparse(substitute(x)), # axis titles
                              caption="", # caption 
                              margin=c(NA, NA, NA, NA),  # margin controls
                              SeasonTrend=list(name="", what="lines",
                                color="green"), 
                              TrendLine=list(name="", what="lines",
                                color="blue"), ...) { # Trends
  if(dev.cur() == 1)
    setGD("SeriesPlot")
  ## Set reverse option for y-axis, needed as default and other defaults
  ytitle <- ytitle
  yaxis.rev <- FALSE
  xlabels <- x$nseasons
  TrendLine$slope <- x$estimate[1]
  if(is.list(ylabels))
    yax <- c(list(data=x$series, axis.range=yaxis.range, axis.log=yaxis.log,
                  axis.rev=yaxis.rev), ylabels)
  else
    yax <- list(data=x$series, axis.range=yaxis.range, axis.log=yaxis.log,
                axis.rev=yaxis.rev, axis.labels=ylabels)
  yax <- do.call("setAxis", yax)
  y <- yax$data
  yax <- yax$dax
  if(length(xlabels) == 1)
    xlabels <- seq(xlabels)
  xlabels <- as.character(xlabels)
  xax <- namePretty(xlabels, orientation="grid")
  ## set margins and controls
  margin.control <- setMargin(margin, yax)
  margin <- margin.control$margin
  right <- margin.control$right
  top <- margin.control$top
  left <- margin.control$left
  bot <- margin.control$bot
  par(mar=margin)
  ## Set up the defaults for the lines and explanation:
  SeasonLine <- setDefaults(SeasonLine, name="", what="vertical", color="black")
  if(SeasonLine$what == "vertical")
    type <- "h"
  else if(SeasonLine$what == "lines")
    type <- "l"
  else {
    warning('invalid value for what; set to "verical" in SeasonLine')
    SeasonLine$what <- "vertical"
    type <- "h"
  }
  SeasonPoint <- setPlot(SeasonPoint, name="", what="points", symbol="circle", 
                             filled=TRUE, size=0.09, color="black")
  SeasonTrend <- setDefaults(SeasonTrend, name="", what="lines", color="green")
  ## Replace NA with mean so that the complete plot is made
  N <- length(y)
  xseas <- rep(seq(along=xlabels), length.out=N)
  yst <- tapply(y, xseas, function(xx) ifelse(is.na(xx), mean(xx, na.rm=TRUE), xx))
  yst <- as.vector(do.call(rbind, yst))
  monthplot(yst, labels=xlabels, type="h", xlim=xax$range, xaxs="i", axes=FALSE,
       ylim=yax$range, yaxs="i", ylab="", xlab="", box=FALSE)
  if(SeasonPoint$what == "none")
    explan <- setExplan(setPlot(list(), name="", what=SeasonLine$what, type="solid",
                                width="standard", symbol="circle", filled=TRUE,
                                size=0.09, SeasonLine$color))
  else { # Set the explanation and draw the points
    explan <- setExplan(SeasonPoint)
    Nseas <- length(xlabels)
    xseq <- rep(seq(Nseas), N %/% Nseas) - 0.45 + seq(0,1, length.out=N) * 0.9
    points(xseq, y, pch=explan$current$pch, cex=explan$current$cex,
           col=explan$current$col, bg=explan$current$col)
  }
  box(lwd=frameWt())
  ## Label the axes
  renderY(yax, lefttitle=ytitle, left=left, right=right)
  renderX(xax, bottitle=xtitle, bottom=bot, top=top, caption=caption)
  ## Pack y into individual series
  y <- split(y, rep(xlabels, length.out=length(y)))
  retval <- (list(y=y, yaxis.log=yaxis.log, yaxis.rev=yaxis.rev,
                  xaxis.log=FALSE, explanation=explan, margin=margin))
  ## Add lines if requested
  ## Suppress any log transforms for next section
  retval$yaxis.log <- retval$xaxis.log <- FALSE
  explan <- retval$explanation
  par(lwd=stdWt())
  ## Add the trend lines if requested
  TrendLine <- setDefaults(TrendLine,name="", what="lines", color="blue")
  if(TrendLine$what != "none") {
    for(i in seq(along=xlabels)) {
      ymean <- mean(y[[xlabels[i]]], na.rm=TRUE)
      segments(i - .45, ymean - .45*TrendLine$slope,
               i + .45, ymean + .45*TrendLine$slope, col=TrendLine$color)
    }
    if(TrendLine$name != "") # add to explanation
      explan <- setExplan(setPlot(TrendLine), explan)
  }
  if(SeasonTrend$what != "none") {
    for(i in seq(along=xlabels)) {
      ydat <- y[[xlabels[i]]]
      ydiff <- outer(ydat, ydat, "-")
      xdiff <- outer(seq(along=ydat), seq(along=ydat), "-")
      slope <- median(ydiff[lower.tri(ydiff)] / xdiff[lower.tri(xdiff)], na.rm=TRUE)
      ymean <- mean(ydat, na.rm=TRUE)
      segments(i - .45, ymean - .45*slope,
               i + .45, ymean + .45*slope, col=SeasonTrend$color)
    }
    if(SeasonTrend$name != "") # add to explanation
      explan <- setExplan(setPlot(SeasonTrend), explan)
  }
  ## recover the log-transforms and explanation if necessary
  retval$yaxis.log <- yaxis.log
  retval$yax <- yax
  retval$xax <- xax
  retval$explanation <- explan
  invisible(retval)
}
