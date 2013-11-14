
print.senSlope <- function(x, digits=4, ...) {
  ## Coding history:
  ##    2013Apr15 DLLorenz Original Coding
  ##
  ## Arguments:
  ##  x (a senSlope object) the object to print
  ##  digits (integer scalar) the number of significant digits to print
  ##  ... (dots) any other arguments to print
  ##
  cat("Call:\n")
  dput(x$call)
  res <- sort(x$residuals) # get the actual residuals regardless of missing values
  cat("\nSmallest 5 residuals:\n")
  print(res[1:5], digits=digits, ...)
  cat("Largest 5 residuals:\n")
  print(rev(res)[1:5], digits=digits)
  cat("\nCoefficients:\n")
  coef <- x$coefficients
  names(coef) <- c("(Intercept)", x$var.name[2L])
  print(coef, digits=digits, ...)
  cat("Confidence intervals for the Sen slope:\n")
  print(x$slope.CI)
  if(!is.null(x$na.action)) {
    n.na <- length(x$na.action)
    if(n.na > 1)
      cat("  (", n.na, " observations deleted due to missing values)\n", sep="")
    else
      cat("  (", n.na, " observation deleted due to missing values)\n", sep="")
  }
  cat("\n")
  invisible(x)
}
