#' Print Objects
#' 
#' Prints the results of a MOVE.2 analysis (\code{move.2}).
#' 
#' 
#' @param x an object of class "move.2" from \code{move.2}.
#' @param digits the number of significant digits to print numeric data.
#' @param \dots additonal arguments for printing numeric values.
#' @return The object \code{x} is returned invisibly.
#' @note The printed output contains the call, the regression coefficients, and
#' selected statistics of the variables.
#' @export
#' @method print move.2
print.move.2 <- function(x, digits=4, ...) {
  ##
  cat("Call:\n")
  dput(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, digits=digits, ...)
  cat("\nStatistics of the variables:\nResponse (", x$var.names[1], "):\n",
      sep="")
  print(x$ystats.orig, digits=digits, ...)
  # This code and that for xstats nicely aligns the tables
  tmp <- capture.output(print(x$ystats, digits=digits, print.gap=2, ...))
  cat(" ",tmp[1L], "\n ", tmp[2L], "\n", sep="")
  cat("Predictor (", x$var.names[2], "):\n", sep="")
  print(x$xstats.orig, digits=digits, ...)
  tmp <- capture.output(print(x$xstats, digits=digits, print.gap=10, ...))
  cat("         ",tmp[1L], "\n         ", tmp[2L], "\n", sep="")
  cat("Correlation coefficient: ", round(x$R, digits), 
      "\n                p-value: ", round(x$p.value, digits), "\n", sep="")
  ## Need different summary of lengths 
  cat("\nConcurrent record length: ", x$N1,
      "\n  Extended record length: ", x$N2, "\n", sep="")
  invisible(x)
}
