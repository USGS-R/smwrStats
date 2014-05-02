# R function for doing ancova model selection and diagnostics
#
# Coding history:
#    2012Jul20 DLLorenz Original coding
#    2012Aug28 DLLorenz Change from diagPlot to plot
#    2013Apr09 DLLorenz Added setGD to plot
#    2103May10 DLLorenz changed anova to Anova
#    

ancovaReg <- function(object, find.best=TRUE, trace=FALSE) {
  ## Arguments:
  ##  object (lm object) the regression object for the diagnostics
  ##  find.best (logical) use step to find the "best" subset
  ##
  ## Required USGS core functions: ppcc.test, vif
  ##
  ## Output value is a list with the following components:
  ##      aovtab     (anova table of the complete model)
  ##      parmests    (parameter estimates, std.errors, etc.)
  ##      vif        (variance inflation factors)
  ##      diagstats  (residual diagnostic statistics (cooks D, leverage, etc) 
  ##      crit.val   (critical values for diagnostic statistics)
  ##      flagobs    ("flagged" observations that one or more diagnostic 
  ##                   statistics exceeding the critical value)
  ##      object     (the final regression object)
  ##      x          (the predictor matrix, without the intercept column)
  ##
  aovtab <- Anova(object) # Anova from car does type II
  if(find.best) {
    if(trace)
      cat("Stepwise elimintation of terms from original model\n")
    object <- step(object, direction="backward", trace=as.integer(trace))
  }
  ## Get the factor variable
  factor.var <- attr(object$terms, "dataClasses")
  factor.var <- names(factor.var)[factor.var == "factor"]
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
  cvdfit <- 2*sqrt(p/n)
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
  if(length(factor.var) > 0L) {
    stuff.out <- list(aovtab=aovtab, parmests=regsum, vif=vif,
                      diagstats=diagstats, crit.val=cvs, flagobs=pck,
                      object=object, x=x, factor.var=factor.var,
                      x.fr=model.frame(object)[, -1L, drop=FALSE]) # drop response
    oldClass(stuff.out) <- "ancovaReg"
  }
  else {
    stuff.out <- list(aovtab=aovtab, parmests=regsum, vif=vif,
                      diagstats=diagstats, crit.val=cvs, flagobs=pck,
                      object=object, x=x)
    oldClass(stuff.out) <- "multReg" # No factors, default to MLR
  }
  return(stuff.out)
}

plot.ancovaReg <- function(x, which='All', set.up=TRUE, span=0.8, ...) {
  ## Arguments:
  ##  x (ancovareg object) the object to plot
  ##  which (character or integer) which graphs to plot
  ##  set.up (logical) call setPage?
  ##  span (numeric scalar) smoothing parameter for loess.smooth
  ## 
  ## Identify which plots to do:
  ## 1 Fitted - Actual
  ## 2 Fitted - Residual
  ## 3 S-L; colors by factor.var
  ## 4 correlogram (only if dates are available)
  ## 5 Q - normal and boxplot by factor.var
  ## 6 Influence
  ## 7 Outliers
  ## 8 Residual dependence plots; colors by factor.var
  ##
  ## Set up graphics page if requested
  if(set.up)
    setGD("ANCOVA")
  ## Set up to do all plots
  doPlot <- rep(TRUE, 8) 
  if(is.numeric(which)) {
    if(min(which) > 0) # select which to plot
      doPlot[seq(8)[-which]] <- FALSE
    else # which not to plot
      doPlot[-which] <- FALSE
  }
  ## Anything else produces all plots
  ## Final plot (8) residuals vs predictors
  ## Extract weights and other data needed later
  Fits <- na.omit(fitted(x$object)) # needed for later plots
  wt.showCD <- if(is.null(x$object$weights))
    rep(1,length(Fits))
  else 
    x$object$weights
  xpred <- x$x.fr # Need the frame to preserve the factor
  if(doPlot[8L]) {
    xnames <- dimnames(xpred)[[2L]]
    xnames <- xnames[xnames != x$factor.var]
    ## Residual dependence plots
    if(length(xnames) > 0L) {
      for(i in xnames) {
        ## Need to protect against things like fourier and quadratic, which result in matrix entries
        if(is.matrix(xpred[, i])) {
          for(j in seq(ncol(xpred[,i]))) {
            xtopl <- xpred[, i][, j]
            AA <- colorPlot(xtopl, x$diagstats$stud.res,
                        color=xpred[[x$factor.var]],
                        Plot=list(what='points', size=0.05),
                        xtitle=paste(i, colnames(xpred[, i])[j], sep=""), 
                        ytitle='Studentized Residual',
                        margin=c(NA, NA, 2.4, .5))
            if(span > 0) {
              smo <- loess.smooth(xtopl, x$diagstats$stud.res, family='sym', span=span)
              addXY(smo$x, smo$y)
            }
            refLine(horizontal=0, Plot=list(what='lines', width='standard',
                           type='dashed'))
            ## The p-value of the second order fit on the residuals exactly matches
            ## the p-value of adding the second order term to the regression
            nl.p <- summary(lm(x$diagstats$resids ~ poly(xtopl, 2),
            									 weights=wt.showCD, model=TRUE), FALSE)$coefficients[3L,4L]
            addTitle(Main=paste("Second order polynomial test for linearity: p=",
            										round(nl.p,4), sep=''), Bold=FALSE)
            ## Must be after Title
            addExplanation(AA, "ul")
          }
        } else {
          AA <- colorPlot(xpred[,i], x$diagstats$stud.res,
                        color=xpred[[x$factor.var]],
                        Plot=list(what='points', size=0.05),
                        xtitle=i, ytitle='Studentized Residual',
                        margin=c(NA, NA, 2.4, .5))
          if(span > 0) {
            smo <- loess.smooth(xpred[,i], x$diagstats$stud.res, family='sym', span=span)
            addXY(smo$x, smo$y)
          }
          refLine(horizontal=0, Plot=list(what='lines', width='standard',
                             type='dashed'))
          ## The p-value of the second order fit on the residuals exactly matches
          ## the p-value of adding the second order term to the regression
          nl.p <- summary(lm(x$diagstats$resids ~ poly(xpred[,i], 2),
          									 weights=wt.showCD, model=TRUE), FALSE)$coefficients[3L,4L]
          addTitle(Main=paste("Second order polynomial test for linearity: p=",
          										round(nl.p,4), sep=''), Bold=FALSE)
          ## Must be after title
          addExplanation(AA, "ul")
        }
      }
    }
  } # end of Residual dependence plots
  ## Studentized vs fitted
  if(doPlot[7L]) {
    tval <- qt(.975, x$object$df.residual-1)
    ylim=max(tval, abs(range(x$diagstats$stud.res))) * 1.05
    xyPlot(Fits, x$diagstats$stud.res,
           Plot=list(what='points', size=0.04),
           yaxis.range=c(-ylim, ylim),
           xtitle='Fitted', ytitle='Studentized Residual')
    refLine(horizontal=0, Plot=list(what='lines', width='standard', type='dashed'))
    refLine(horizontal=c(-tval, tval),
            Plot=list(what='lines', width='color', type='dashed', color='red'))
  }
  ## Show the effect that each flagged point has on the regression
  ## (an influence plot)
  Act <- Fits + residuals(x$object) # needed later
  if(doPlot[6L]) {
    xyPlot(Fits, Act,
           Plot=list(what='points', size=0.02),
           xtitle="Fitted", ytitle="Actual")
    refLine(coefficients=c(0,1), Plot=list(what='lines', width='standard', type='dashed'))
    if(is.null(x$object$na.action))
      fg.showCD <- which(x$flagobs)
    else 
      fg.showCD <- which(naresid(x$object$na.action, x$flagobs))
    for(j in seq(along=fg.showCD)) {
      i <- fg.showCD[j]
      ## generate a random color
      Col <- sprintf('#%02x%02x%02x', as.integer(runif(1, 10, 245)),
                     as.integer(runif(1, 10, 245)),
                     as.integer(runif(1, 10, 245)))
      refLine(coefficients=lsfit(Fits[-i], Act[-i], wt.showCD[-i])$coef,
              Plot=list(what='lines', color=Col, width='color'))
      addXY(Fits[i], Act[i],
            Plot=list(what='points', color=Col, size=0.04))
      labelPoints(Fits[i], Act[i], dimnames(x$x)[[1L]][i],
                  dir='NE', size=12, color=Col)
    }
  } # end of Influence plot
  ## Q-normal of standardized residuals (H&H criterion 4)
  if(doPlot[5L]) {
    ## First is box plot by factor level
    boxPlot(x$diagstats$stnd.res, group=x$x.fr[[x$factor.var]],
            Box=list(type="tukey"))
    probPlot(x$diagstats$stnd.res,
             yaxis.log=FALSE, ylabels=7,
             ytitle='Standardized Residual',
             margin=c(NA, NA, 2.4, NA), mean=0, sd=1)
    PPCC <- ppcc.test(x$diagstats$stnd.res)$p.value
    addTitle(Main=paste("PPCC test for normality: p=", round(PPCC,4), sep=''), Bold=FALSE)
  }
  ## for the next plots use Pearson residuals, which are weighted
  Res <- residuals(x$object, type="pearson")
  ## If possible, plot a correlogram--requires finding 1 datelike column in
  ## the data
  if(doPlot[4L] && !is.null(x$object$call$data)) {
    data <- x$object$call$data
    if(is.name(data))
      data <- try(eval(data))
    if(class(data) != "try-error") {
      anyDate <- which(sapply(data, isDateLike))
      if(length(anyDate) == 1) { # if more than 1, 
        Date <- dectime(data[[anyDate]])
        if(!is.null(skips <- x$object$na.action)) {
          if(attr(skips, "class") == "omit")
            Date <- Date[-skips] # Must remove if class is omit 
        }
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
           Plot=list(what='points', size=0.05),
           xtitle="Fitted",
           ytitle=as.expression(substitute(sqrt(abs(YL)),
               list(YL = as.name("Residuals")))),
           margin=c(NA, NA, 2.4, NA))
    if(span > 0) {
      smo <- loess.smooth(Fits, sqrt(abs(Res)), span=span)
      addXY(smo$x, smo$y)
    }
    refLine(horizontal=0.82218*sqrt(RSE), Plot=list(what='lines', width='standard', type='dashed'))
    Woodings <- cor.test(Fits, abs(Res), method='s', exact=FALSE)
    addTitle(Main=paste("Woodings test for heteroscedasticity: p=",
               round(Woodings$p.value,4), sep=''), Bold=FALSE)
  } # end of S-L 
  ## 2nd plot response vs. fit
  if(doPlot[2L]) {
    xyPlot(Fits, Res,
           Plot=list(what='points', size=0.05),
           xtitle="Fitted",
           ytitle="Residuals")
    if(span > 0) {
      smo <- loess.smooth(Fits, Res, span=span)
      addXY(smo$x, smo$y)
    }
    refLine(horizontal=0, Plot=list(what="lines", width="standard", type="dashed"))
  }
  ## First plot is actual vs partial fitted, with regression details
  if(doPlot[1L]) {
    ## Construct fitted values using only the continuous variables
    Levels <- levels(xpred[[x$factor.var]])
    ## This plot requires the input data set
    data <- x$object$call$data
    if(is.null(data)) {
      warning("Plot # 1 requires the original data")
      return(invisible(x))
    }
    else if(is.name(data)) {
      Par.fr <- try(eval(data))
      if(class(data) == "try-error") {
        warning("Plot # 1 requires the original data")
        return(invisible(x))
      }
      ## We must protect against finding a symbol with the same name as an object in this function!
      ## Why here and not for doPlot4?
      if(is.name(Par.fr)) {
        Par.fr <- as.character(Par.fr)
        Par.fr <- get(Par.fr, pos=1L)
      }
    }
    else # Must be a dataset
      Par.fr <- data
    ## This makes a dummy data set for factor = first level
    Par.fr[[x$factor.var]] <- Levels[1L]
    ## Make predictions and remove missings to match data
    Par.fits <- predict(x$object, Par.fr)
    if(!is.null(toDel <- x$object$na.action))
      Par.fits <- Par.fits[-toDel]
    AA <- colorPlot(Par.fits, Act, color=xpred[[x$factor.var]],
                    Plot=list(what="points", size=0.05),
                    xtitle=paste("Partial Fitted: ", x$factor.var, " dropped",
                      sep=""),
                    ytitle="Response")
    Colors <- setColor(Levels) # The default for colorPlot
    for(i in seq(along=Levels)) {
      Par <- lm(Act ~ Par.fits, subset=xpred[[x$factor.var]] == Levels[i])
      refLine(coefficients=coefficients(Par),
              Plot=list(what="lines", width="color", color=Colors[i]),
              xrange=range(Par.fits[xpred[[x$factor.var]] == Levels[i]]))
    }
    addExplanation(AA, where="lr")
    ## Add some details, regression eqn and RSE
    Eqn <- x$object$coef
    names(Eqn)[1L] <- ""
    Eqn <- paste(as.character(round(Eqn, 3)), names(Eqn), sep=' ')
    Eqn <- paste(Eqn, collapse=' + ')
    Resp <- as.expression(substitute(hat(R), 
        list(R=as.name(deparse(x$object$call$formula[[2L]])))))
    Model <- expression()
    RSE <- signif(RSE, 3)
    legend("topleft", legend=c(
                        as.expression(substitute(hat(R) == EQN, 
                            list(R=as.name(deparse(x$object$call$formula[[2L]])),
                                 EQN=Eqn))),
                        paste("Residual Standard Error: ", RSE, sep='')), bty='n')
  }
  invisible(x)
}

print.ancovaReg <- function(x, digits=3, ...) {
  ## Arguments:
  ##  x (ancovareg object) the object to print
  ##  digits (integer scalar) how many digits to use when prining
  ##   ... (dots) not used, required for method function
  ##
  cat("\nOriginal model\n")
  print(x$aovtab, digits=digits, signif.stars=FALSE)
  cat("\n\n\nFinal model\n")
  print(x$parmests, digits=digits, signif.stars=FALSE)
  if(length(x$vif) > 1) {
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
    print(x$diagstats[x$flagobs,], digits=digits)
  }
  else
    cat("\tNo observations exceeded any test criteria\n")
  invisible(x)
}
