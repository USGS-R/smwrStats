# Qunatiles using type=2
#
# Coding History:
#    2012Sep17 DLLorenz original Coding
#    2014May08 DLLorenz update for R 3.1
#

quantile.numeric <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE,
                             names = TRUE, type = 2, ...)
  ## Wrapper to force default type = 2 for all USGS applications
	NextMethod(, type=type)
