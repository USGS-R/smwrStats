# Generic confirm function
#
# Coding history:
#    2011Jul06 DLLorenz Original coding
#    2011Jul06          This version
#

confirm <- function(x, ...)
  UseMethod("confirm")

confirm.default <- function(x, ...) {
  warning(paste("No known confirm method for object of class ", class(x), sep=''))
  invisible(x)
}
