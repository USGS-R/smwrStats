#' Sample Quantiles
#' 
#' Computes sample quantiles corresponding to the given probabilities: 
#'method for "numeric" data. The smallest observation corresponds to a probability 
#'of 0 and the largest to a probability of 1. This method function is a simple 
#'wrapper for the default function that sets the default \code{type} to 2 for
#'numeric data.
#' 
#' @param x numeric vector whose sample quantiles are wanted.
#' @param probs numeric vector of probabilities with values in the range from 0
#' to 1.
#' @param na.rm remove missing values \code{NA}s before computation?
#' @param names include names of the probabilities, expressed as percentages?
#' @param type an interger between 1 and 9 that selects the method for
#' computing the quantile. See \bold{Note}.
#' @param \dots further arguments passed to or from other methods.
#' @return An optionally named vector corresponding to the quantiles of
#' \code{x} for the selected probabilities.
#' @note Helsel and Hirsh (2002) define the 75th percentile as "a value which
#' exceeds no more than 75 percent of the data and is exceeded by no more than
#' 25 percent of the data." This rule can be easily extended to other
#' percentiles. The selection of \code{type} equal to 2 ensures that this rule
#' is met for all data. The rule stated by Helsel and Hirsch is very useful for
#' an emprical description of the data, but Hyndman and Fan (1996) describe the
#' slection of \code{type} for other applications.
#' @seealso \code{\link{quantile.default}}
#' @references Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in
#' water resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' 
#' Hyndman, R.J. and Fan, Y. 1996, Sample quantiles in statistical packages:
#' American Statistician, v. 50, p. 361-365.
#' @keywords univar
#' @examples
#' 
#' # The default value for type (7) will compute a value that is exceeded by 4 values 
#' #  for a sample of size 14
#' quantile(seq(14), type=7)
#' # But 4/14 is greater than 25 percent. But setting type to 2 will result in
#' #  only 3 values that are larger than the computed 75th percentile.
#' quantile(seq(14))
#' 
#' @importFrom stats quantile
#' @export
#' @method quantile numeric
quantile.numeric <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE,
                             names = TRUE, type = 2, ...)
	# Coding History:
	#    2012Sep17 DLLorenz original Coding
	#    2014May08 DLLorenz update for R 3.1
	#    2014Dec29 DLLorenz Convert to roxygen headers
  ## Wrapper to force default type = 2 for all USGS applications
	NextMethod(, type=type)
