# create a data frame of summary statistics, with optional grouping vars
#
# Coding history:
#    2009Mar04 DLLorenz Original Coding
#    2009Mar17 DLLorenz Debug error in retgrp
#    2009Nov04 DLLorenz Fix percentile labels
#    2010Feb19 DLLorenz Fix for variable labels, added check for numerics
#    2010Feb19          Added to the USGS library 4.0
#    2011Aug09 DLLorenz Conversion to R
#    2012Feb03          This version.
#

sumStats <- function(..., group=NULL, Num="Num",
                     Stats=list(Mean=mean, StdDev=sd),
                     Probs=(0:4)/4, na.rm=TRUE) {
  ## Arguments:
  ##  ... (any number of vectors or a data.frame or a mnatrix) the data to
  ##    be summarized
  ##  group (a vector or data.frame or list) the data to form unique groups
  ##  Stats (tagged list) The output column name and the function to compute
  ##    the statistic
  ##  Probs (numeric vector or NULL) the arguemtn to quantile specifying the
  ##    desired probabilities.
  ##  na.rm (logical scalar) remove missing values?
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
