# R function for doing multiple regression diagnostics
#
# Coding history:
#    ????????? AVecchia Original coding for Env. Stats. class
#    2007Apr04 DLLorenz Modified for library
#    2007Jun05 DLLorenz Added span argument for plot function
#    2008Apr21 DLLorenz changed comment from parests to parmests
#    2008Apr24 DLLorenz Bug fix
#    2008Aug14 DLLorenz Added Wooding's test to plot
#    2008Oct20 DLLorenz Added PPCC test for normality
#    2010Nov12 DLLorenz Conversion to R
#    2011Apr26 DLLorenz Added polynomial nonlinearity test to partial plots
#    2011Jun29 DLLorenz Bug fix and name change
#    2011Oct25 DLLorenz Update for package
#    2012Mar22 DLLorenz Bug fix in diagnostic plots for the correlogram
#    2012Aug28 DLLorenz Change from diagPlot to plot
#    2012Dec27 DLLorenz Minor tweaks, including integer indexes
#    2013Jan08 DLLorenz Used Anova in car rather than anova for type II Sum Sq.
#    2013Apr09 DLLorenz Added setGD to plot
#    2013May10 DLLorenz Use na.omit to work with both na.omit and na.exclude in call
#    2013Dec30 DLLorenz Added * to print of observations exceeding test criterion
#    2014May14 DLLorenz Increased criterion for dfits

multReg <- function(object) {
  ## Arguments:
  ##  object (lm object) the regression object for the diagnostics
  ##
  ## Required USGS core functions: ppcc.test, vif
  ##
  ## Output value is a list with the following components:
  ##      aovtab     (anova table)
  ##      parmests    (parameter estimates, std.errors, etc.)
  ##      vif        (variance inflation factors)
  ##      diagstats  (residual diagnostic statistics (cooks D, leverage, etc) 
  ##      crit.val   (critical values for diagnostic statistics)
  ##      flagobs    ("flagged" observations that one or more diagnostic 
  ##                   statistics exceeding the critical value)
  ##      object     (the original regression object)
  ##      x          (the predictor matrix, without the intercept column)
  ##
  aovtab <- Anova(object) # Anova from car does type II
  regsum <- summary(object)
  vif <- vif(object)
  ## Compute the diagnostics using lsfit
  ## get x matrix
  if(is.null(object[['x']])) # to avoid conflict with the xlevels component.
    x <- model.matrix(object)[, -1L, drop=FALSE]
  else
    x <- object$x[, -1L, drop=FALSE]
  ## get y
  if(is.null(object$y))
    y <- model.extract(model.frame(object), "response")
  else
    y <- object$y
  respvar <- make.names(object$terms[[2L]])
  if(length(respvar) > 1L) # happens when function is used
    respvar <- paste(respvar, collapse='.')
  ## get weights
  if(is.null(object$weights))
    wt <- rep(1, length(y))
  else
    wt <- object$weights
  ##  Run lsfit and ls.diag to get the residual diagnostic stats we need:
  lsout <- lsfit(x, y, wt)
  lsdiag <- ls.diag(lsout)
  lev <- lsdiag$hat
  cooksd <- lsdiag$cooks
  std.res <- lsdiag$std.res
  dfits <- lsdiag$dfits
  stud.res <- lsdiag$stud.res
  ##  Compute critical values for diagnostic statistics and pick out offending
  ##  observations
  p <- ncol(x)+1
  n <- nrow(x)
  cvlev <- 3*p/n
  cvdfit <- qgrubbs(0.01, n)*sqrt(p/n)
  cvcook <- qf(.5,p+1,n-p)
  pck <- c(lev>cvlev | cooksd>cvcook | abs(dfits)>cvdfit)
  ##  Combine the diagnostic stats into a single data set, round it to 3
  ##  decimals, rename y, and combine the output into a list.
  yhat <- object$fitted.values
  diagstats <- data.frame(y=y, yhat=yhat, resids=lsout$residuals,
                          stnd.res=std.res, stud.res=stud.res,
                          leverage=lev, cooksD=cooksd, dfits=dfits)
  names(diagstats)[1L] <- respvar
  cvs <- c(cvlev,cvcook,cvdfit)
  names(cvs) <- c('leverage','cooksD','dfits')
  stuff.out <- list(aovtab=aovtab, parmests=regsum, vif=vif,
                    diagstats=diagstats, crit.val=cvs, flagobs=pck,
                    object=object, x=x)
  oldClass(stuff.out) <- "multReg"
  return(stuff.out)
}

plot.multReg <- function(x, which='All', set.up=TRUE, span=1.0, ...) {
  ## Arguments:
  ##  x (multreg object) the object to plot
  ##  which (character or integer) which graphs to plot
  ##  span (numeric scalar) smoothing parameter for loess.smooth
  ##  set.up and  ... (dots) unused, required for method function
  ## 
  ## Identify which plots to do:
  ## 1 Fitted - Actual
  ## 2 Fitted - Residual
  ## 3 S-L
  ## 4 correlogram (only if dates are available)
  ## 5 Q - normal
  ## 6 Influence
  ## 7 Outliers
  ## 8 Residual dependence plots
  ##
  ## Set up graphics page
  if(set.up)
    setGD("MULTREG")
  ## Set up to do all plots
  doPlot <- rep(TRUE, 8L)
  do8 <- FALSE
  if(is.numeric(which)) {
    if(min(which) > 0) # select which to plot
      doPlot[seq(8L)[-which]] <- FALSE
    else # which not to plot
      doPlot[-which] <- FALSE
  }
  else if(is.character(which) && which[1] != "All") {
    doPlot[1:7] <- FALSE
    xnames <- which
    do8 <- TRUE
  }
  ## Anything else produces all plots
  ## Final plot (8) residuals vs predictors
  ## Extract weights and other data needed later
  Fits <- na.omit(fitted(x$object)) # needed for later plots, but without missing values
  wt.showCD <- if(is.null(x$object$weights))
    rep(1,length(Fits))
  else
    x$object$weights
  ## Protect aginst intercept only fit
  if(diff(range(Fits)) < 1.e-15) {
    doPlot[-5L] <- FALSE
    doPlot[5L] <- TRUE
  }
  if(doPlot[8L]) {
    xpred <- x$x
    if(!do8) # get all explanatory variable names
      xnames <- dimnames(xpred)[[2L]]
    ## Residual dependence plots
    if(do8 || length(xnames) > 1L) {
      for(i in xnames) {
        ## Protect against possible binary predictors as created by factors
        if(length(unique(xpred[,i])) < 3L)
          warning(i, " appears to be created by a factor variable, not plotted")
        else {
          xyPlot(xpred[,i], x$diagstats$stud.res,
                 Plot=list(what="points", size=0.05),
                  xtitle=i, ytitle="Studentized Residual",
                 margin=c(NA, NA, 1.5, .5))
          if(span > 0) {
            smo <- loess.smooth(xpred[,i], x$diagstats$stud.res, family="sym", span=span)
            addXY(smo$x, smo$y)
          }
          refLine(horizontal=0, Plot=list(what="lines", width="standard", type="dashed"))
          ## The p-value of the second order fit on the residuals exactly matches
          ## the p-value of adding the second order term to the regression
          nl.p <- summary(lm(x$diagstats$resids ~ poly(xpred[,i], 2),
                             weights=wt.showCD, model=TRUE), FALSE)$coefficients[3,4]
          addTitle(Main=paste("Second order polynomial test for linearity: p=",
                     round(nl.p, 4L), sep=""), Bold=FALSE)
        }
      }
    }
  } # end of Residual dependence plots
  ## Studentized vs fitted
  if(doPlot[7L]) {
    tval <- qt(.975, x$object$df.residual - 1)
    ylim=max(tval, abs(range(x$diagstats$stud.res))) * 1.05
    xyPlot(Fits, x$diagstats$stud.res,
           Plot=list(what="points", size=0.04),
           yaxis.range=c(-ylim, ylim),
           xtitle="Fitted", ytitle="Studentized Residual")
    refLine(horizontal=0, Plot=list(what="lines", width="standard", type="dashed"))
    refLine(horizontal=c(-tval, tval),
            Plot=list(what="lines", width="color", type="dashed", color="red"))
  }
  ## Show the effect that each flagged point has on the regression
  ## (an influence plot)
  Act <- Fits + na.omit(residuals(x$object)) # needed later
  if(doPlot[6L]) {
    xyPlot(Fits, Act,
           Plot=list(what="points", size=0.02),
           xtitle="Fitted", ytitle="Actual")
    refLine(coefficients=c(0,1), Plot=list(what="lines", width="standard", type="dashed"))
    if(is.null(x$object$na.action))
      fg.showCD <- which(x$flagobs)
    else 
      fg.showCD <- which(naresid(x$object$na.action, x$flagobs))
    for(j in seq(along=fg.showCD)) {
      i <- fg.showCD[j]
      ## generate a random color
      Col <- sprintf("#%02x%02x%02x", as.integer(runif(1, 10, 245)),
                     as.integer(runif(1, 10, 245)),
                     as.integer(runif(1, 10, 245)))
      refLine(coefficients=lsfit(Fits[-i], Act[-i], wt.showCD[-i])$coef,
              Plot=list(what="lines", color=Col, width="color"))
      addXY(Fits[i], Act[i],
            Plot=list(what="points", color=Col, size=0.04))
      labelPoints(Fits[i], Act[i], dimnames(x$x)[[1L]][i],
                  dir="NE", size=12, color=Col)
    }
  } # end of Influence plot
  ## Q-normal of standardized residuals (H&H criterion 4)
  if(doPlot[5L]) {
    qqPlot(x$diagstats$stnd.res, Plot=list(size=0.05),
             yaxis.log=FALSE, ylabels=7,
             ytitle="Standardized Residual",
             xtitle="Standard Normal Quantiles",
             margin=c(NA, NA, 2.4, NA))
    refLine(coefficients=c(0,1))
    PPCC <- ppcc.test(x$diagstats$stnd.res)$p.value
    addTitle(Main=paste("PPCC test for normality: p=", round(PPCC,4), sep=""), Bold=FALSE)
  }
  ## for the next plots use Pearson residuals, which are weighted
  Res <- na.omit(residuals(x$object, type="pearson"))
  ## If possible, plot a correlogram--requires finding 1 datelike column in
  ## the data
  if(doPlot[4L] && !is.null(x$object$call$data)) {
    data <- x$object$call$data
    if(is.name(data))
      data <- try(eval(data))
    if(class(data) != "try-error") {
      anyDate <- which(sapply(data, isDateLike))
      if(length(anyDate) == 1L) { # if more than 1, skip it
        Date <- dectime(data[rownames(data) %in% rownames(xpred),anyDate])
        corGram(Date, Res)
      }
    }
  }
  ## Add details of call on regression model to next plots
  Mod <- format(x$object$call$formula)
  ## 3rd plot, S-L
  RSE <- rmse(x$object)
  if(doPlot[3L]) {
    xyPlot(Fits, sqrt(abs(Res)),
           Plot=list(what="points", size=0.05),
           xtitle="Fitted",
           ytitle=as.expression(substitute(sqrt(abs(YL)),
               list(YL = as.name("Residuals")))),
           margin=c(NA, NA, 2.4, NA))
    if(span > 0) {
      smo <- loess.smooth(Fits, sqrt(abs(Res)), span=span)
      addXY(smo$x, smo$y)
    }
    ## 0.82218 is the expected value of the sqrt(abs(x)) for a normal dist.:
    ## integrate(function(x) sqrt(abs(x))*dnorm(x), -Inf, Inf)
    refLine(horizontal=0.82218*sqrt(RSE), Plot=list(what="lines", width="standard", type="dashed"))
    Woodings <- cor.test(Fits, abs(Res), method="s", exact=FALSE)
    addTitle(Main=paste("Woodings test for heteroscedasticity: p=",
               round(Woodings$p.value,4), sep=""), Bold=FALSE)
  } # end of S-L 
  ## 2nd plot response vs. fit
  if(doPlot[2L]) {
    xyPlot(Fits, Res,
           Plot=list(what="points", size=0.05),
           xtitle="Fitted",
           ytitle="Residuals")
    if(span > 0) {
      smo <- loess.smooth(Fits, Res, span=span)
      addXY(smo$x, smo$y)
    }
    refLine(horizontal=0, Plot=list(what="lines", width="standard", type="dashed"))
  }
  ## First plot is actual vs fitted, with regression details
  if(doPlot[1L]) {
    xyPlot(Fits, Act,
           Plot=list(what="points", size=0.05),
           xtitle=paste("Fitted:", Mod, sep=" "),
           ytitle="Response")
    if(span > 0) {
      smo <- loess.smooth(Fits, Act, span=span)
      addXY(smo$x, smo$y)
    }
    refLine(coefficients=c(0,1), Plot=list(what="lines", width="standard", type="dashed"))
    ## Add some details, regression eqn and RSE
    Eqn <- x$object$coef
    names(Eqn)[1L] <- ""
    Eqn <- paste(as.character(round(Eqn, 3)), names(Eqn), sep=" ")
    Eqn <- paste(Eqn, collapse=" + ")
    Resp <- as.expression(substitute(hat(R), 
        list(R=as.name(deparse(x$object$call$formula[[2]])))))
    Model <- expression()
    RSE <- signif(RSE, 3)
    legend("topleft", legend=c(
                        as.expression(substitute(hat(R) == EQN, 
                            list(R=as.name(deparse(x$object$call$formula[[2]])),
                                 EQN=Eqn))),
                        paste("Residual Standard Error: ", RSE, sep="")), bty="n")
  }
  invisible(x)
}

print.multReg <- function(x, digits=3, ...) {
  ## Arguments:
  ##  x (binaryreg object) the object to print
  ##  digits (integer scalar) how many digits to use when prining
  ##   ... (dots) not used, required for method function
  ##
  print(x$parmests, digits=digits, signif.stars=FALSE)
  ## Add model comparison stats:
  cat("press: ", signif(press(x$object), digits),
    "\n  AIC: ", signif(AIC(x$object), digits),
    "\n  BIC: ", signif(BIC(x$object), digits),
    "\n\n", sep="")
  if(length(x$vif) > 1L) {
    print(x$aovtab, digits=digits, signif.stars=FALSE) # really only needed for MLR
    cat("\nVariance inflation factors\n")
    namvif <- format(names(x$vif), justify="left")
    valvif <- format(round(x$vif, 2), justify="right")
    for(i in seq(along=x$vif))
    	cat(namvif[i], " ", valvif[i], "\n", sep="")
  }
  cat("\nTest criteria\n")
  print(x$crit.val, digits=digits)
  if(any(x$flagobs)) {
    cat("\tObservations exceeding at least one test criterion\n")
    dstats <- format(x$diagstats[x$flagobs,] , digits=digits)
    ## Append * to each value that exceeds its criterion
    dstats$leverage <- paste(dstats$leverage, 
                             ifelse(x$diagstats[x$flagobs, "leverage"] > x$crit.val[1L],
                                    "*", " "), sep="")
    dstats$cooksD <- paste(dstats$cooksD, 
                             ifelse(x$diagstats[x$flagobs, "cooksD"] > x$crit.val[2L],
                                    "*", " "), sep="")
    dstats$dfits <- paste(dstats$dfits, 
                             ifelse(abs(x$diagstats[x$flagobs, "dfits"]) > x$crit.val[3L],
                                    "*", " "), sep="")
    print(dstats)
  }
  else
    cat("\tNo observations exceeded any test criteria\n")
  invisible(x)
}
