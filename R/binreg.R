# diagnostics for logistic regression
#
# Coding history:
#    2009Jan26 DLLorenz Original Coding and start of tweaks
#    2010Jan29 DLLorenz Modified to allow for matrix responses
#    2011Aug22 DLLorenz Conversion to R
#    2011Oct25 DLLorenz Update for package
#    2012Aug28 DLLorenz Change from diagPlot to plot
#    2013Apr09 DLLorenz Added setGD to plot
#

binaryReg <- function(object, lc.max=1000) {
  ## Argument:
  ##  object (a glm model object) the logistic regression model
  ##
  ## Get the model frame and do some preliminary processing
  obj.frame <- model.frame(object)
  resp.var <- obj.frame[[1]]
  resp.val <- object$y
  ## The response profile describe how the response data were coded
  if(class(resp.var) == 'matrix')
    resp.profile <- quantile(resp.val)
  else {
    resp.profile <- table(resp.var, resp.val)
    resp.profile <- cbind(expand.grid(dimnames(resp.profile)),
                          Counts=as.vector(resp.profile))
    resp.profile <- resp.profile[resp.profile$Counts > 0,]
    names(resp.profile)[1:2] <- c(names(obj.frame[1]), "Response")
    rownames(resp.profile) <- as.character(seq(nrow(resp.profile)))
  }
  ## Compute the summary and Wald stats are added to the coefs matrix
  obj.sum <- summary(object)
  ## Check to make sure that contrasts were specified if any factors
  classes <- sapply(obj.frame, class)
  factors <- classes[-1] %in% c('factor', 'ordered')
  Warning=''
  Factors=list()
  if(any(factors)) {
    classes <- names(classes[-1])[factors] # keep only those that are factors
    if(is.null(object$call$contrasts)) {
      if(length(classes) > 1)
        Warning=paste("\nWarning:\n", paste(classes, collapse=', '),
          "were included as factor explanatory variables, but no contrasts specified\n")
      else
        Warning=paste("\nWarning:\n", classes,
          "was included as a factor explanatory variable, but no contrasts specified\n")
    } # end of if null contrasts
    else {
      contrs <- names(object$call$contrasts)
      noclasses <- classes[!classes %in% contrs]
      if(length(noclasses) > 1)
        Warning=paste("\nWarning:\n", paste(noclasses, collapse=', '),
          "were included as factor explanatory variables, but no contrasts specified\n")
      else if(length(noclasses) == 1)
        Warning=paste("\nWarning:\n", noclasses,
          "was included as a factor explanatory variable, but no contrasts specified\n")
      ## No warning if length is 0, create class info
      for(i in classes[classes %in% contrs]) {
        info <- eval(call(object$call$contrasts[[i]], levels(obj.frame[[i]])))
        info.col <- colnames(info)
        if(length(info.col) ==0)
          info.col <- as.character(seq(ncol(info)))
        colnames(info) <- paste(i, info.col, sep='')
        Factors[[i]] <- info
      } # end of processing factors
    } # end of else null contrasts
  } # end of if any factors
  ## do the le Cessie and Houwelingen and Hosmer-Lemeshow tests
  if(class(resp.var) == 'matrix' || length(resp.var) > lc.max)
    leCessie <- NULL
  else
    leCessie <- leCessie.test(object)
  ## H-L only if there are more than 15 unique predicted values
  resp.fit <- fitted(object)
  if(length(unique(resp.fit)) > 15)
    HL <- hosmerLemeshow.test(object)
  else
    HL <- NULL
  ## Compute the Percent correct and concordant stats and AUROC
  if(class(resp.var) == 'matrix') {
    PctCorrect <- NULL
    Concord <- NULL
  }
  else {
    n <- length(resp.fit)
    PctCorrect <- c("1" = 100*sum(resp.fit[resp.val == 1] > .5)/sum(resp.val),
                    "0" = 100*sum(resp.fit[resp.val == 0] < .5)/sum(1-resp.val))
    Concord <- table(sign(outer(resp.fit[resp.val == 1],
                                resp.fit[resp.val == 0], '-')))
    names(Concord) <- c("Discordant", "Tied", "Concordant")[match(names(Concord), c("-1", "0",  "1"))]
  }
  rocout <- roc(object)
  ## Compute some influence stats using the method in H-L and Frank Harrell
  ##  get x, y, and weights
  x <- model.matrix(object)
  y <- fitted(object)*(1-fitted(object))
  wt <- object$weights
  ## Run lsfit and ls.diag to get the residual diagnostic stats we need:
  lsout <- lsfit(x, y, wt, intercept=FALSE)
  lsdiag <- ls.diag(lsout)
  lev <- lsdiag$hat
  cooksd <- lsdiag$cooks
  dfits <- lsdiag$dfits
  ## Compute critical values for diagnostic statistics and pick out offending
  ##  observations
  p <- ncol(x)+1
  n <- nrow(x)
  cvlev <- 3*p/n
  cvdfit <- 2*sqrt(p/n)
  cvcook <- qf(.5,p+1,n-p)
  pck <- c(lev>cvlev | cooksd>cvcook | abs(dfits)>cvdfit)
  ## Combine the diagnostic stats into a single data set, round it to 3
  ##  decimals, rename y, and combine the output into a list.
  yhat <- object$fitted.values
  diagstats <- data.frame(y=resp.var,
                          yhat=yhat, resids=resid(object, 'response'),
                          deviance.res=resid(object, 'deviance'),
                          pearson.res=resid(object, 'pearson'),
                          leverage=lev, cooksD=cooksd, dfits=dfits)
  respvar <- make.names(object$terms[[2]])
  if(length(respvar) > 1) # happens when function is used
    respvar <- paste(respvar, collapse='.')
  names(diagstats)[1] <- respvar
  cvs <- c(cvlev,cvcook,cvdfit)
  names(cvs) <- c('leverage','cooksD','dfits')
  ## pack it up
  retval <- list(regsum=obj.sum, Warning=Warning, Factors=Factors,
                 Profile=resp.profile, Hosmer=HL,
                 leCessie=leCessie, PctCorrect=PctCorrect,
                 Concordance=Concord, roc=rocout,
                 diagstats=diagstats,
                 crit.val=cvs, flagobs=pck, object=object)
  oldClass(retval) <- "binaryreg"
  return(retval)
}

print.binaryreg <- function(x, digits=4, ...) {
  ## Arguments:
  ##  x (binaryreg object) the object to print
  ##  digits (integer scalar) how many digits to use when prining
  ##   ... (dots) any additional arguments to print
  ##
  ## Print the summary, any warning, the factor information, and G2
  print(x$regsum, digits=digits, ...)
  if(x$Warning != "")
    cat(x$Warning)
  if(length(x$Factors) > 0) {
    cat("\nFactor level information:\n")
    lapply(x$Factors, print)
  }
  G2 <- round(x$regsum$null.deviance - x$regsum$deviance, digits)
  df <- x$regsum$df[1]
  cat("\nLikelihood ratio test: ", G2, " on ", df-1,
      " degrees of freedom, p-value is ", round(1-pchisq(G2, df-1), digits),
      "\n\n", sep='')
  ## Print the response matrix and the le Cressie and Houwelingen test
  cat("Response profile:\n")
  print(x$Profile, digits=digits, ...)
  cat("\n Goodness of fit tests\n")
  if(!is.null(x$leCessie))
    print(x$leCessie, digits=digits, ...)
  else
    cat("Le Cessie-Van Houwlingen test not computed.\n")
  if(!is.null(x$Hosmer))
    print(x$Hosmer, digits=digits, ...)
  else
    cat("Too few unique predcited values for Hosmer-Lemeshow Test\n")
  ## Print the correct and concordant stats and the AUROC
  cat("\nPredictive power estimates:\n")
  ## Print the R2 and adjusted R2
  R2 <- round(1 - x$object$deviance/x$object$null.deviance, digits)
  adjR2 <- round(1 - (x$object$aic - 2)/x$object$null.deviance, digits) # need to correct for intercept term
  cat("McFadden R-squared: ", R2, "\nadjusted R-squared: ", adjR2, "\n\n", sep="")
  if(!is.null(x$PctCorrect)) {
    cat("\nClassification table.\nPercent correct: (1 is sensitivity, 0 is specificity)\n")
    print(x$PctCorrect, digits=3)
    nConcord <-  sum(x$Concordance)
    cat("\nConcordance Index, based on",
        nConcord, "pairs\n")
    print(100*x$Concordance/nConcord, digits=digits, ...)
  }
  else
    cat("\nOther predictive power estimates cannot be computed for martix reponses.\n")
  print(x$roc)
  ## Print the diagnostics
    cat("\nInfluence diagnostic test criteria:\n")
  print(x$crit.val, digits=digits, ...)
  if(any(x$flagobs)) {
    cat("\tObservations exceeding at least one test criterion\n")
    print(x$diagstats[x$flagobs,], digits=digits, ...)
  }
  else
    cat("\tNo observations exceeded any test criteria\n")
  invisible(x)
}

plot.binaryreg <- function(x, which='All', set.up=TRUE, bandw=0.3, ...) {
  ## Arguments:
  ##  x (binaryreg object) the object to plot
  ##  which (character or integer) which graphs to plot
  ##  bandw (numeric scalar) bandwidth for kernel smoothing for H-L graph
  ##
  ## Set up graphics page
  if(set.up) {
    setGD("LOGISTIC")
  }
  ## Set up to do all 5 plots
  doPlot <- rep(TRUE,5) 
  if(is.numeric(which)) {
    if(min(which) > 0) # select which to plot
      doPlot[seq(5)[-which]] <- FALSE
    else # which not to plot
      doPlot[-which] <- FALSE
  }
  ## Anything else produces all plots
  ## 
  ## Last are the L-W-Y partial plots
  if(doPlot[5]) {
    x.resid <- resid(x$object, 'response')
    x.coefs <- coef(x$object)
    x.mat <- model.matrix(x$object)
    for(i in seq(2, length(x.coefs))) {
      xs <- x.mat[,i]
      ## Skip if binary
      if(length(unique(xs)) < 3) next
      rtoplot <- x.resid[order(xs)]
      xs <- sort(xs)
      rtoplot <-  cumsum(rtoplot)/sqrt(length(xs))
      xyPlot(xs, rtoplot, Plot=list(what='stairstep'),
             xtitle=names(x.coefs[i]),
             ytitle='L-W-Y Cumulative Residuals', yaxis.range=c(-.6, .6))
      refLine(horizontal=0)
    }
  }
  ## ROC is plot # 4
  if(doPlot[4])
    plot(x$roc, set.up=FALSE)
  ## Third plot--classification graph
  if(doPlot[3] && !is.null(x$PctCorrect)) {
    pd.0 <- density(x$object$fitted.values[x$object$y == 0], bw=0.25, kern='tri')
    pd.1 <- density(x$object$fitted.values[x$object$y == 1], bw=0.25, kern='tri')
    ylim <- c(0, max(pd.0$y, pd.1$y))*1.1
    AA <- xyPlot(pd.0$x, pd.0$y, Plot=list(name="specificity", color='darkblue'),
                 xaxis.range=c(0,1), yaxis.range=ylim,
                 xtitle='Predicted', ytitle='Relative Density', margin=c(NA, NA, 1.2, NA))
    AA <- addXY(pd.1$x, pd.1$y, Plot=list(name="sensitivity", color='darkgreen'), current=AA)
    addExplanation(AA, where='ur', title='')
    refLine(vertical=0.5)
    addTitle(paste("specificity:", round(x$PctCorrect[["0"]], 3),
                "  sensitivity:", round(x$PctCorrect[["1"]], 3), sep=' '), Bold=FALSE)
  }
  ## Second Plot, overall fit with H-L
  if(doPlot[2] && !is.null(x$Hosmer)) {
    fits <- fitted(x$object)
    y <- x$object$y
    xyPlot(fits, y, Plot=list(what='points', filled=FALSE),
           xtitle='Fitted Values', ytitle='Observed values',
           xlabels=list(labels=5, extend=5), ylabels=list(labels=5, extend=5),
           margin=c(NA, NA, 1.2, NA))
    p2.smo <- ksmooth(fits, y, bandwidth=bandw, kernel='normal')
    addXY(range(fits), range(fits), Plot=list(what='lines'))
    addXY(p2.smo$x, p2.smo$y, Plot=list(what='lines'))
    HLx <- x$Hosmer$estimate[,2] / x$Hosmer$estimate[,1]
    HLy <- x$Hosmer$estimate[,3] / x$Hosmer$estimate[,1]
    addXY(HLx, HLy, Plot=list(what='points', filled=TRUE))
    addTitle(paste("Hosmer-Lemeshow Test, p-value =", round(x$Hosmer$p.value,4)), Bold=FALSE)
  }
  ## First plot (last on page), leCessie
  if(!is.null(x$leCessie)) {
    plot(x$leCessie, which=1, set.up=FALSE)
    addTitle(paste("le Cessie-van Houwelingen Test, p-value =", round(x$leCessie$p.value,4)),
    				 Bold=FALSE)
  }
  invisible(x)
}
