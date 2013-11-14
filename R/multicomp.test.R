# perform multicomp test 
#
#   Coding History
#    2007Apr24 DLLorenz Original Coding of multicompRank.test
#    2011Aug06 DDLorenz Conversion to R
#    2012Jan27 DLLorenz Begin conversion to a single MCT for all possibilities
#                        becuase the basic computations are all the same.
#    2012Jul24 DLLorenz Bug fix in lMat to prevent collapse to vector
#    2012Sep21 DLLorenz Added Size to output of associations
#    2013May15 DLLorenz Fix tbl assoc for inconsistent results
#

multicomp.test <- function (x, g, method="parametric", critical.value="", alpha=0.05) {
  ## Arguments:
  ##  x (numeric vector) The response data
  ##  g (grouplike vector) The groups for x
  ##  method (character scalar) the method to use
  ##  critical.value (character scalar) the method to compute the critical value
  ##  alpha (numeric scalar) the probability level
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
      overvar <- N*(N-1)/12
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
  retval <- list(title=title, cv.method=critical.value, alpha=alpha, crit.value=crit.val,
                 response=xname, groups=gname, means=groupMean, sizes=groupNum,
                 table=tbl, assoc=lMat)
  oldClass(retval) <- "MCT"
  return(retval)
}

print.MCT <- function(x, digits=4, ...) {
  ## print function for objects of class MCT
  cat("\t", x$title, "\n")
  if(x$cv.method == "lsd")
    cat("Pairwise")
  else
    cat("Overall")
  cat(" error rate: ", x$alpha, "\nCritical value: ", round(x$crit.value, digits),
      " by the ", x$cv.method, " method\n\n", sep="")
  cat("Response variable: ", x$response, "\nGroup variable: ", x$groups,
      "\n\n", sep="")
  cat("Table of paired comparisons, ", round(1 - x$alpha, 4) * 100,
      " percent confidence intervals\n excluding 0 are flagged by *.\n", sep="")
  pmat <- format(signif(x$table, digits))
  pmat <- cbind(pmat, flag=ifelse(x$table[,3]*x$table[,4] > 0, "*", " "))
  print(pmat, quote=FALSE)
  cat("\nTable of associations among groups\n")
  pmat=cbind(Mean=signif(x$means, digits), Size=x$sizes, x$assoc)
  print(pmat, quote=FALSE)
  cat("\n")
  invisible(x)
}
