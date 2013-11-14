# print a seasonalPeak
#
# Coding history:
#    2007Aug22 DLLorenz Original Coding ?
#    2011May25 DLLorenz Begin Conversion to R and rename
#    2013Apr02 DLLorenz  Final tweaks for release
#

print.seasonalPeak <- function(x, digits=3, details=FALSE, ...) {
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
