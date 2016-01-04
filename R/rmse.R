#' Root-Mean-Squared and Relative Differences
#' 
#' Computes the root-mean-squared error (RMSE) of the difference between
#' observed values and the predicted values or the RMSE or relative percent
#' differences (RPD) between samples and duplicates.
#' 
#' 
#' @aliases rmse rmse.default rmse.lm rpd
#' @param x either a random vector an object for which a method exists.
#' @param y duplicate samples paired with \code{x}.
#' @param plotit logical, if \code{TRUE}, then create a Bland-Altman 
#'mean-difference plot (banld and Altman, 1986); otherwise no plot is created.
#' @param \dots arguments to be passed to or from methods.
#' @return For the \code{rmse} functions, a single value representing the
#' estimated RMSE. For \code{rpd}, the relative percent differences for each
#' paired sample and duplicate.
#' @note The definition for the RMSE of paired water-quality duplicates is
#' \deqn{ RMSE = \sqrt{\frac{\sum{(x_i - y_i)^2}}{2n}}}{RMSE = sqrt(sum of
#' squared differences/2 times the number of duplicates)} The definition for
#' RPD for paired water-quality duplicates is \deqn{ RPD = abs(x - y)/(x + y)/2
#' * 100}{RPD = abs(x - y)/(x + y)/2 * 100} Other disciplines may use different
#' equations.
#' @references 
#' Bland J.M. and Altman D.G., 1986 Statistical methods for assessing agreement 
#'between two methods of clinical measurement: Lancet, i, p. 307--310.\cr
#'
#' Clesceri, L.S., Greenberg, A.E., and Eaton, A.D., 1998, Standard
#'methods for the examination of water and wastewater, 20th edition:
#'Baltimore, Md, United Book Press, Inc., 1162 p.\cr

#'
#' Harvey, D., undated, Analytical chemistry 2.0: Analytical Sciences Digital
#'Library: online at URL:
#'http://www.asdlib.org/onlineArticles/ecourseware/Analytical%20Chemistry%202.0/Welcome.html\cr
#' 
#' Helsel, D.R., and Hirsch, R.M., 2002, Statistical methods in water
#'resources: U.S. Geological Survey Techniques of Water-Resources
#'Investigations, book 4, chap. A3, 522 p.\cr
#' @keywords univar
#' @examples
#' 
#' # Example 15.2 from Harvey.
#' dupX1 <- c(160, 196, 207, 185, 172, 133)
#' dupX2 <- c(147, 202, 196, 193, 188, 119)
#' rmse(dupX1, dupX2)
#' rpd(dupX1, dupX2)
#' 
#' @export rmse
rmse <- function(x, ...) {
	# Coding history:
	#   2009Oct02 DLLorenz Original Coding
	#   2012May24 DLLorenz Conversion to R
	#   2014Jul28 DLLorenz Bug fix to rpd.
	#   2014Dec29 DLLorenz Convert to roxygen headers
  UseMethod("rmse")
}

#' @name rmse
#' @export
#' @method rmse default
rmse.default <- function(x, y, ...) {
  ##
 sqrt(mean((x-y)^2)/2) # from Analytical Chemistry
}

#' @name rmse
#' @export
#' @method rmse lm
rmse.lm <- function(x, ...) {
  ## Arguments:
  ##  x, an object that inherit class lm, such as aov.
  ##
  rdf <- x$df.resid
  if(!is.null(x$weights))
    w <- x$weights
  else
    w <- 1.0
  sqrt(sum(w * x$residuals^2)/rdf)
}

#' @name rmse
#' @export rpd
rpd <- function(x, y, plotit=FALSE) {
  xname <- deparse(substitute(x))
  yname <- deparse(substitute(x))
  dif <- x - y
  summ <- x + y
  if(plotit){ # Create the Bland-Altman (Tukey m-d) plot
  	xyPlot(summ/2, dif, Plot=list(size=0.05), 
  				 xtitle=paste("Mean of", xname, "and", yname, sep=" "),
  				 ytitle=paste0(xname, " - ", yname))
  	refLine(horizontal=0)
  }
  ##
 abs(dif)/summ * 200 # from standard methods
}
