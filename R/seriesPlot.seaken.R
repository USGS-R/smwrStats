#' Series Plot
#' 
#' Produces a plot of time-series data by season that shows seasonal and overall
#' trends.
#' 
#' The argument \code{what} for \code{SeasonLine} must be either "lines" or
#' "vertical."  See \code{\link{monthplot}} for more information.\cr The
#' argument \code{what} for either \code{SeasonPoint}, \code{SeasonTrend}, or
#' \code{trendLine} may be set to "none" to suppress drawing of that feature.
#' 
#' @aliases seriesPlot seriesPlot.seaken
#' @param x data that can be treated as a regularly-spaced time series.
#' @param SeasonLine control parameters of the lines in the plot. See
#' \bold{Details}.
#' @param SeasonPoint control parameters of the points in the plot. See
#' \bold{Details}.
#' @param yaxis.log log-transform the y axis?
#' @param yaxis.range set the range of the y axis.
#' @param ylabels set the y-axis labels. See \code{\link{linearPretty}} for
#' details.
#' @param xlabels set the x-axis labels, may be a single numeric value
#' indicating the number of season in \code{x}. See \code{\link{namePretty}}
#' for details.
#' @param xtitle the x-axis title.
#' @param ytitle the y-axis title.
#' @param caption the figure caption.
#' @param margin the parameters of the margin of the plot area.
#' @param SeasonTrend control parameters of the trend line for each season. See
#' \bold{Details}.
#' @param TrendLine control parameters of the overall trend line compute by
#' \code{seaken}. See \bold{Details}.
#' @param \dots any additional arguments required for specific methods.
#' @return Information about the graph.
#' @seealso \code{\link{monthplot}}, \code{\link{seasonPlot}},
#' \code{\link{seaken}}
#' @keywords hplot
#' @export
#' @method seriesPlot seaken
seriesPlot.seaken <- function(x, # data
															SeasonLine=list(name="", what="vertical", color="black"),
															SeasonPoint=list(name="", what="points", symbol="circle", 
																							 filled=TRUE, size=0.09, color="black"), # plot controls
															yaxis.log=FALSE, yaxis.range=c(NA,NA), # y-axis controls
															ylabels=7, xlabels=x$seasonames, # labels
															xtitle="",
															ytitle=deparse(substitute(x)), # axis titles
															caption="", # caption 
															margin=c(NA, NA, NA, NA),  # margin controls
															SeasonTrend=list(name="", what="lines",
																							 color="brown"), 
															TrendLine=list(name="", what="lines",
																						 color="blue"), ...) { # Trends
	if(dev.cur() == 1)
		setGD("SeriesPlot")
	## Set reverse option for y-axis, needed as default and other defaults
	ytitle <- ytitle
	yaxis.rev <- FALSE
	if(is.null(xlabels)) {
		xlabels <- x$nseasons
	}
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
	xax <- namePretty(xlabels, orientation="grid", style="between")
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
	SeasonTrend <- setDefaults(SeasonTrend, name="", what="lines", color="brown")
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
	# Scale set to represent change over x-range of .9 to represent nyears of time
	scl <- 0.5/0.9*x$nyears 
	TrendLine <- setDefaults(TrendLine,name="", what="lines", color="blue")
	if(TrendLine$what != "none") {
		for(i in seq(along=xlabels)) {
			ymean <- mean(y[[xlabels[i]]], na.rm=TRUE)
			segments(i - .45, ymean - scl*TrendLine$slope,
							 i + .45, ymean + scl*TrendLine$slope, col=TrendLine$color, 
							 lwd=lineWt("bold"))
		}
		if(TrendLine$name != "") # add to explanation
			explan <- setExplan(setPlot(TrendLine, width="bold"), explan)
	}
	# Defaults set earlier
	if(SeasonTrend$what != "none") {
		for(i in seq(along=xlabels)) {
			ydat <- y[[xlabels[i]]]
			ydiff <- outer(ydat, ydat, "-")
			xdiff <- outer(seq(along=ydat), seq(along=ydat), "-")
			slope <- median(ydiff[lower.tri(ydiff)] / xdiff[lower.tri(xdiff)], na.rm=TRUE)
			ymean <- mean(ydat, na.rm=TRUE)
			segments(i - .45, ymean - scl*slope,
							 i + .45, ymean + scl*slope, col=SeasonTrend$color, 
							 lwd=lineWt("bold"))
		}
		if(SeasonTrend$name != "") # add to explanation
			explan <- setExplan(setPlot(SeasonTrend, width="bold"), explan)
	}
	## recover the log-transforms and explanation if necessary
	retval$yaxis.log <- yaxis.log
	retval$yax <- yax
	retval$xax <- xax
	retval$explanation <- explan
	invisible(retval)
}
