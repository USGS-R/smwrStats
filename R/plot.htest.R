#' Diagnostic Plots
#' 
#' Creates diagnostic plots for selected hypothesis tests.
#' 
#' The \code{kensen} method for \code{plot} graphs the y and t data with the
#' best fit line.
#' 
#' The \code{ppcc} method creates a single graph that shows the q-normal plot
#' with a line showing the fit.
#' 
#' @aliases plot.htest plot.kensen plot.ppcc
#' @param x an object having classes "htest" and some other class for which a
#' \code{plot} method exists.
#' @param which either "All" or a number indicating which plot, see
#' \bold{Details}.
#' @param set.up set up the graphics page?
#' @param \dots not used, required for other methods.
#' @return The object \code{x} is returned invisibly.
#' @keywords hplot
#' @export
#' @method plot htest
plot.htest <- function(x, which="All", set.up = TRUE, ...) {
  ## Coding history:
  ##    2013Jan06 DLLorenz Initial Coding and start ppw
  ##    2014Dec29 DLLorenz Conversion to roxygen headers
  NextMethod("plot")
}

#' @rdname plot.htest
#' @export
#' @method plot kensen
plot.kensen <- function(x, which="All", set.up = TRUE, ...) {
  ## Identify which plots to do: which is ignored
  ## 
  ## Set up graphics page
  if(set.up) 
    setGD("KenSen")
  ## Set up to do the plots, ignore which
  ## X and Y axis titles
	titles <-strsplit(x$data.name, " and ")[[1L]]
  timePlot(x$t, x$y, xtitle=titles[2L], ytitle=titles[1L])
  refLine(coefficients=x$coef)
  invisible(x)
}

#' @rdname plot.htest
#' @export
#' @method plot ppcc
plot.ppcc <- function(x, which="All", set.up = TRUE, ...) {
	## Identify which plots to do: which is ignored
	## 
	## Set up graphics page
	if(set.up) 
		setGD("PPCC")
	## Set up to do the plots, ignore which
	qqPlot(x$x, ytitle=x$data.name)
	invisible(x)
}
