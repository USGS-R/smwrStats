# Qunatiles using type=2
#
# Coding History:
#    2012Sep17 DLLorenz original Coding
#    2012Sep17          This version.
#

quantile.numeric <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE,
                             names = TRUE, type = 2, ...)
  ## Wrapper to force default type = 2 for all USGS applications
  quantile.default(x, probs = probs, na.rm = na.rm,
                   names = names, type = type, ...)
