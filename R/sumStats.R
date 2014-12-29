#' Compute Summary Statistics
#' 
#' Create a dataset of specified summary statistics from a collection of data.
#' 
#' 
#' @param \dots any number of arguments, which can be of many different forms:
#' a dataset, selected columns of a dataset, vectors, and matrices. If a
#' dataset is supplied, then nonnumeric columns are removed before computing
#' statistics.
#' @param group list whose components are interpreted as categories, each of
#' the same length as the objects in \dots{} The e lements of the categories
#' define the position in a multi-way array corresponding to each observation.
#' Missing values (NAs) are allowed. The names of group are used as the names
#' of the columns in the output dataset.  If a vector is given, it will be
#' treated as a list with one component. The default is NULL, which indicates
#' no grouping variable.
#' @param Num a character string indicating the name of the column that
#' contains the number of nonmissing observations.
#' @param Stats a tagged list. An element in the list should have the name of
#' the target variable in the output data set, a nd it should be a function
#' that accepts the na.rm argument and computes a single value.  See
#' \bold{Notes} for commonly used functions. The default is the target columns
#' Mean and StdDev, which are computing using functions mean and stdev.
#' @param Probs vector of desired probability levels. Values must be between 0
#' and 1. Minimum returned for probs=0 and maximum returned for probs=1.
#' Default is c(0.0, 0.25, 0.50, 0.75, 1.0).
#' @param na.rm logical; if TRUE, then missing values are removed from each
#' column in \dots{} before computing the statistics
#' @return A data frame containing columns identifying each variable in
#' \dots{}, any grouping variables, and the requested statistics as named in
#' the call.
#' @note Commonly used functions referenced in the Stats arguments include
#' mean, sd, skew and var.\cr
#' 
#' The statistics requested by the Probs argument are computed by the quantile
#' function using \code{type}=2.
#' @seealso \code{\link{mean}}, \code{\link{sd}}, \code{\link{skew}},
#' \code{\link{var}}, \code{\link{quantile}},
#' @references Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in
#' water resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.
#' @keywords univar
#' @examples
#' 
#' ## Generate a random sample
#' set.seed(222)
#' XX.rn <- rexp(32)
#' sumStats(XX.rn)
#' 
#' @export sumStats
sumStats <- function(..., group=NULL, Num="Num",
                     Stats=list(Mean=mean, StdDev=sd),
                     Probs=(0:4)/4, na.rm=TRUE) {
	# Coding history:
	#    2009Mar04 DLLorenz Original Coding
	#    2009Mar17 DLLorenz Debug error in retgrp
	#    2009Nov04 DLLorenz Fix percentile labels
	#    2010Feb19 DLLorenz Fix for variable labels, added check for numerics
	#    2010Feb19          Added to the USGS library 4.0
	#    2011Aug09 DLLorenz Conversion to R
	#    2014Dec29 DLLorenz Conversion to roxygen header
  ##
  ## Note: only functions that take the na.rm can be passed to this
  ## function. 
  dots <- list(...)
  ndg <- length(dots)
  if(ndg == 1) {
    dotname <- deparse(substitute(...))
    dots <- dots[[1]]
    if(mode(dots) == "numeric") {
      dots <- list(dots)
      names(dots) <- dotname
    }
  }
  else { # multiple vectors were specified
    dotname <- as.list(match.call())
    ## Drop named components
    dotname$Num <- dotname$group <- dotname$Stats <- dotname$Probs <- dotname$na.rm <- NULL
    dotname <- sapply(dotname, deparse)
    names(dots) <- dotname[-1] # drop the call to sumStats
  }
  dots <- dots[sapply(dots, is.numeric)]
  ndg <- length(dots)
  ## Fix names of Probs so that they can be converted to a data.frame
  if(!is.null(Probs) && length(Probs) > 0)
    ProbNames <- paste("Pct", round(Probs * 100,
                                    if(length(Probs) > 1) 2 - log10(diff(range(Probs)))
                                    else 2), sep = ".")
  else
    ProbNames = ""
  if(!is.null(group)) {
    groups <- interaction(group, drop=TRUE)
    retval <- by(as.data.frame(dots), INDICES=groups, FUN=sumStats, Num=Num,
                 Stats=Stats, Probs=Probs, na.rm=na.rm)
    retgrp <- by(as.data.frame(group, stringsAsFactors=FALSE),
                 INDICES=groups, FUN=function(x, n) {
      xx <- x[1,,drop=FALSE]
      if(n > 1) {
        xx <- as.data.frame(lapply(xx, rep, times=n), stringsAsFactors=FALSE)
      }
      xx	
    }, n=ndg)
    retval <- do.call("rbind", retval)
    ## The by functions appears to strip a single column data frame of its class
    if(class(retgrp[[1]]) != 'data.frame') {
      retgrp <- as.matrix(unlist(retgrp))
      colnames(retgrp) <- "Group" # assign a simple name
    }
    else
      retgrp <- do.call("rbind", retgrp)
    retval <- cbind(retgrp, retval)
    ## Fix case where there is only one variable--results in Variable == dots
    if(ndg == 1 && !is.null(dotname))
      levels(retval$Variable) <- dotname
  }
  else { # no grouping
    retval <- lapply(dots, function(x, na.rm, Num, Stats, Probs, ProbNames) {
      retpart <- double()
      if(!is.null(Num) && Num != "") {
        retpart <- sum(!is.na(x))
        names(retpart) <- Num
      }
      if(!is.null(Stats) && length(Stats) > 0)
        retpart <- c(retpart, sapply(Stats, function(fcn, x, na.rm)
                                     fcn(x, na.rm=na.rm), x=x, na.rm=na.rm))
      if(!is.null(Probs) && length(Probs) > 0) {
        retprob <-  quantile(x, probs=Probs, na.rm=na.rm, type=2)
        names(retprob) <- ProbNames
        retpart <- c(retpart, retprob)
      }
      retpart
    } # end of function
                     ,na.rm=na.rm, Num=Num, Stats=Stats, Probs=Probs,
                     ProbNames=ProbNames)
    VarNames <- names(retval)
    retval <- as.data.frame(do.call("rbind", retval))
    if(!is.null(VarNames))
      retval <- cbind(Variable=VarNames, retval)
  } # end of else
  row.names(retval) <- seq(nrow(retval))
  return(retval)
}
