# function to compute terms for predicting piecewise modeling of trends
#
# Coding history:
#    2013Aug15 DLLorenz Original coding
#

makepredictcall.trends <- function(var, call) {
  if (as.character(call)[1L] != "trends") 
    return(call)
  call$breaks <- attr(var, "breaks")
  call$boundary.breaks <- c(TRUE, TRUE) 
  return(call)
}