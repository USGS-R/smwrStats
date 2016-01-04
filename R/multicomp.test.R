#' Multiple Comparisons
#' 
#' Performs multiple comparison tests among groups of data. The tests may be
#'either parametric (Yandell, 1997), nonparametric (Higgins, 2004), or Dunn's 
#'nonparametric (Glantz, 2005).
#' 
#' @details The choices for \code{method} are "parametric," "nonparametric," and 
#'"dunn." If the \code{method} is "parametric," then the comparisons are based on 
#'the means and variances of the raw data and the valid choices for 
#'\code{critical.value} are "tukey" (default), "bonferroni," or "lsd." Otherwise, 
#'the comparisons are based on the ranks of the data. Valid choices for 
#'\code{critical.value} are "tukey" (default), "bonferroni," or "lsd" when 
#'\code{method} is "nonparametric" and "sidak" (default) or "bonferroni" 
#'when \code{method} is "dunn." The basic diffference between the default 
#'nonparametric method and Dunn's nonparametric method is in the handling of ties.
#' 
#' @param x the numeric vector of observations. Missing values (NAs) are
#'allowed and removed before the test is performed.
#' @param g any group vector for the observations. Missing values (NAs) are
#'allowed and removed before the test is performed.
#' @param method a character string describing the test. Only the first
#'character is necessary. See \bold{Details}.
#' @param critical.value a character string describing the method to use for
#'determining the critical value. Only the first character is necessary. See
#' \bold{Details}.
#' @param alpha the significance level of the test. See \bold{Note}.
#' @return An object of class MCT containing the following components:
#' \item{title}{ a description of the test. } \item{cv.method}{ the method used
#' to compute the critical value. } \item{alpha}{ the value of \code{alpha}.  }
#' \item{crit.value}{ the critical value for the pairwise comparisons. }
#' \item{response}{ the name of the response variable.  } \item{groups}{ the
#' name of the group variable.  } \item{means}{ the means for each group. }
#' \item{sizes}{ the number of observations in each group. } \item{table}{ the
#' table of the results of the pairwise comparisons. } \item{assoc}{ a data
#' frame containing the possible association for each group. }
#' @note All computations of the variance for unequal group sizes are based on
#' the harmonic mean as described in Yandell (1997). That adjustment is only
#' approximate when \code{critical.value} is "tukey" and \code{method} is
#' "parametric" but useful when the design is slightly unbalanced.\cr The
#' default nonparametric method \code{method} = "nonparametric" is only
#' assymptotically unbiased when some data are tied. For smaller data sets with
#' small numbers of ties, it may be preferable to use Dunn's nonparametric
#' method \code{method} = "dunn."
#' @references Glantz, S.A., 2005, Primer of biostatistics: McGraw Hill, New
#' York, 520 p.
#' 
#' Higgins, J.J., 2004, Introduction to modern nonparametric statistics: 
#'Pacific Grove, Calif., Brooks/Cole, 384 p.
#'
#' Yandell, B.S., 1997, Practical data analysis for designed experiments:
#'London, United Kingdom, Chapman & Hall, 437 p.
#' @keywords nonparametric htest
#' @export multicomp.test
multicomp.test <- function(x, g, method="parametric", critical.value="", alpha=0.05) {
	#   Coding History
	#    2007Apr24 DLLorenz Original Coding of multicompRank.test
	#    2011Aug06 DDLorenz Conversion to R
	#    2012Jan27 DLLorenz Begin conversion to a single MCT for all possibilities
	#                        becuase the basic computations are all the same.
	#    2012Jul24 DLLorenz Bug fix in lMat to prevent collapse to vector
	#    2012Sep21 DLLorenz Added Size to output of associations
	#    2013May15 DLLorenz Fix tbl assoc for inconsistent results
	#    2014Dec22 DLLorenz Roxygen headers
  ##
  ## valid methods: tukey, bonferroni, and lsd
  ## if lsd is selected, then error.type is comparison wise
  ##
  ## remove NAs from x and g
  xname=deparse(substitute(x))
  gname=deparse(substitute(g))
  sel <- !(is.na(x) | is.na(g))
  x <- x[sel]
  g <- g[sel, drop=T] # get rid of any unused factor levels
  ## Basic preliminary steps
  g <- as.factor(g)
  N <- length(x)
  Levels <- levels(g)
  Ngrp <- length(Levels)
  groupNum <- table(g)
  Ncomp <- choose(Ngrp, 2)
  method <- match.arg(method, c("parametric", "nonparametric", "dunn"))
  if(method == "parametric") {
    groupMean <- tapply(x, g, mean)
    df <- N - Ngrp
    overvar <- sum((x - groupMean[g])^2)/ df
    groupVar <- overvar/groupNum
    if(critical.value == "")
      critical.value <- "tukey"
    else
      critical.value <- match.arg(critical.value, c("tukey", "bonferroni", "lsd"))
    crit.val <- switch(critical.value,
                       tukey=qtukey(1 - alpha, Ngrp, df)/sqrt(2),
                       bonferroni=qt(1 - alpha/(2 * Ncomp), df),
                       lsd=qt(1 - alpha/2, df))
    title <- "Parametric Multiple Comparison Test"
  } else if(method == "nonparametric") {
    r <- rank(x)
    groupMean <- tapply(r, g, mean)
    if(any(duplicated(x))) # Correction as described in Higgins
      overvar <- var(r)
    else
      overvar <- N*(N+1)/12
    groupVar <- overvar/groupNum
    ## Set the critical value
    if(critical.value == "")
      critical.value <- "tukey"
    else
      critical.value <- match.arg(critical.value, c("tukey", "bonferroni", "lsd"))
    crit.val <- switch(critical.value,
                       tukey=qtukey(1 - alpha, Ngrp, Inf)/sqrt(2),
                       bonferroni=qnorm(1 - alpha/(2 * Ncomp)),
                       lsd=qnorm(1 - alpha/2))
    title <- "Nonparametric Multiple Comparison Test"
    xname <- paste("Rank of", xname, sep=" ")
  }
  else { # Dunn's
    r <- rank(x)
    groupMean <- tapply(r, g, mean)
    ## Adjust var for ties
    tieCor <- table(x)
    tieCor <- sum(tieCor^3 - tieCor)
    overvar <- N*(N+1)/12-tieCor/(12*(N-1))
    groupVar <- overvar/groupNum
    ## Set the critical value
    if(critical.value == "")
      critical.value <- "sidak"
    else
      critical.value <- match.arg(critical.value, c("sidak", "bonferroni"))
    if(critical.value == "sidak")
      alpha.adj <- (1 - (1-alpha)^(1/Ncomp))/2 # Sidak"s alpha level
    else
      alpha.adj <- alpha/Ncomp/2 # Bonferroni"s alpha level
    crit.val <- qnorm(1 - alpha.adj)
    title <- "Dunn's Nonparametric Multiple Comparison Test"
    xname <- paste("Rank of", xname, sep=" ")
  }
  ## Sort by means (largest to smallest for pretty output)
  orderMean <- order(-groupMean)
  groupMean <- groupMean[orderMean]
  groupVar <- groupVar[orderMean]
  Levels <- Levels[orderMean]
  ## Assemble tables
  diffMat <- outer(groupMean, groupMean, "-")
  stderrMat <- sqrt(outer(groupVar, groupVar, "+"))
  nameMat <- outer(Levels, Levels, paste, sep="-")
  sel <- lower.tri(nameMat)
  ## The pairwise comparison table, change sign of differences for visual effect
  tbl <- cbind(estimate=-diffMat[sel], stderr=stderrMat[sel],
               lower=-diffMat[sel] - stderrMat[sel]*crit.val,
               upper=-diffMat[sel] + stderrMat[sel]*crit.val)
  dimnames(tbl)[[1]] <- t(nameMat)[sel]
  ## The group association table
  tMat <- abs(diffMat)/stderrMat < crit.val # these are the groups
  ## Now construct the letter matrix
  lMat <- matrix(" ", Ngrp, Ngrp)
  lMat[lower.tri(lMat, TRUE)] <- ifelse(tMat[lower.tri(tMat, TRUE)], "X", " ")
   ## Check for consistent groups (mostly result from wide ranges in N)
  OK <- TRUE
  for(i in seq(3L, Ngrp)) { # First possible mismatch
    X <- FALSE
    for(j in seq(1L, i)) {
      if(!X && lMat[i, j] == "X")
        X <- TRUE
      else if(X && lMat[i, j] ==  " ") {
        OK <- FALSE
        for(jj in seq(1L, j))
          if(lMat[i, jj] == "X")
            lMat[j, jj] <- " "
      }
    }
  }
  ## Collapse columns that match
  toKeep <- rep(TRUE, Ngrp)
  if(OK) {
    for(i in seq(Ngrp-1L)) {
      j <- i+1L
      jseq <- seq(j, Ngrp)
      if(all(lMat[jseq,i] == lMat[jseq,j]))
        toKeep[j] <- FALSE
    }
  } else { # Need to loop back
    for(j in seq(2L ,Ngrp)) {
      for(i in seq(1L, j-1L)) {
        jseq <- seq(j, Ngrp)
        if(all(lMat[jseq,i] == lMat[jseq,j]))
         toKeep[j] <- FALSE
      }
    }
  }
  lMat <- lMat[, toKeep, drop=FALSE]
  dimnames(lMat) <- list(Levels , LETTERS[seq(ncol(lMat))])
  as.data.frame(lMat)
  ## Construct the return value
  ## Fix the match for sizes
  Order <- match(names(groupMean), names(groupNum))
  retval <- list(title=title, cv.method=critical.value, alpha=alpha, crit.value=crit.val,
                 response=xname, groups=gname, means=groupMean, 
  							 sizes=groupNum[Order], table=tbl, assoc=lMat)
  oldClass(retval) <- "MCT"
  return(retval)
}
