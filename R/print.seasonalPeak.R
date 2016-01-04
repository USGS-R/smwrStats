#' Print Objects
#' 
#' Prints the results of a seasonal peak analysis (\code{seasonalPeak}).
#' 
#' 
#' @param x an object of class "seasonalPeak" from \code{seasonalpeak}.
#' @param digits the number of significant digits to print numeric data.
#' @param details logical, print the model details?
#' @param \dots not used for method, required for other methods.
#' @return The object \code{x} is returned invisibly.
#' @note The printed output contains the default time-of-peak value and
#' potential alternate values for unconfirmed peaks and the number of peaks,
#' timing of the primary peak, and optionally the model information for
#' confirmed peaks.
#' @export
#' @method print seasonalPeak
print.seasonalPeak <- function(x, digits=3, details=FALSE, ...) {
	# Coding history:
	#    2007Aug22 DLLorenz Original Coding ?
	#    2011May25 DLLorenz Begin Conversion to R and rename
	#    2013Apr02 DLLorenz Final tweaks for release
	#    2014Dec29 DLLorenz Convert to roxygen header
	#
  if(attr(x, "Confirmed")) {
    cat("Confirmed seasonal peak:\nNumber of peaks: ",
        attr(x, "NumPeaks"), "\n", sep="")
    cat("Time of peak:", round(as.double(x), digits=digits), "\n", sep=" ")
    if(details) {
      cat("\nLoading models:", attr(x, "Models"), "\n", sep=" ")
      if(attr(x, "NumPeaks") == 2L) {
        cat("Second peak models\n")
        print(attr(x, "Second.peak"))
      }
      if(!is.null(Hlives <- attr(x, "hlife")))
        cat("Half lives:", Hlives, sep=" ")
      cat("\n")
    }
    return(invisible(x))
  } # done if confirmed
  Extra <- attr(x, "Extra")
  Extra <- attr(x, "Points")$xout[Extra]
  cat("Unconfirmed seasonal peak:\nDefault value:",
      round(as.double(x), digits=digits), 
      "\nAlternate values:", round(Extra, digits=digits), "\n\n", sep=" ")
  invisible(x)
}
