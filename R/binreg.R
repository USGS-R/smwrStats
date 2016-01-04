#' Diagnostics for Logistic Regression
#' 
#' Computes diagnostic statistics for logistic regression.
#' 
#' Because the le Cessie test is very slow to compute for many observations,
#' the test is not performed if there are more than \code{lc.max} observations.
#' 
#' @param object the logistic regression model object.
#' @param lc.max the lmaximum number of observations for the le Cessie test.
#' @return A list of class "binaryreg" containing these components:
#' \item{regsum}{the output from \code{summary(object)}} \item{Warning}{any
#' warnings relevant to the model} \item{Factors}{information about factor
#' explanatory variables} \item{Profile}{summary information about the coding
#' of the response variable} \item{Hosmer}{output from the Hosmer-Lemeshow test
#' on \code{object}} \item{leCessie}{output from the le Cessie-van Houwelingen
#' test on \code{object}} \item{PctCorrect}{the classification table}
#' \item{Concordance}{the concordance table} \item{roc}{output from the
#' receiver operating characteristics test on \code{object}} \item{diagstats}{a
#' data frame containing the response variable, the predicted response
#' probability, the response residual, the deviance residuals, the Pearson
#' residuals, the leverage, the value of Cook's D, and the dfits value}
#' \item{crit.val}{the critical values for leverage, Cook's D, and dfits}
#' \item{flagobs}{a logical value indicating which observaitons exceeded any
#' one of the critical values} \item{object}{the \code{object}}
#' @note Logistic regression can be very useful alternative method for heavily
#' censored water-quality data.\cr The critical values for the test criteria
#' are computed as: leverage, \emph{3p/n}; Cook's D, median quantile for the
#' \emph{F} distribution with \emph{p+1} and \emph{n-p} degrees of freedonm;
#' and dfits, the .01 quantile of the \emph{grubbs} distribution for \emph{n}
#' observations, where \emph{p} is the number of parameters estiamted in the
#' regression and \emph{n} is the number of observations.\cr Objects of class
#' "binaryreg" have \code{print} and \code{plot} methods.
#' @seealso \code{\link{roc}}, \code{\link{leCessie.test}},
#' \code{\link{hosmerLemeshow.test}}
#' @references Harrell, F.E., Jr., 2001, Regression modeling strategies with
#' applications to linear models, logistic regression and survival analysis:
#' New York, N.Y., Springer, 568 p.\cr
#' 
#' Helsel, D.R., and Hirsch, R.M., 2002, Statistical methods in water
#' resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' 
#' McFadden, D., 1974, Conditional logit analysis of qualitative choice
#' behavior: p. 105-142 in Zarembka, P. (ed.), Frontiers in Econometrics.
#' London, Academic Press, 252 p.
#' @keywords models regression
#' @export binaryReg
binaryReg <- function(object, lc.max=1000) {
	# Coding history:
	#    2009Jan26 DLLorenz Original Coding and start of tweaks
	#    2010Jan29 DLLorenz Modified to allow for matrix responses
	#    2011Aug22 DLLorenz Conversion to R
	#    2011Oct25 DLLorenz Update for package
	#    2012Aug28 DLLorenz Change from diagPlot to plot
	#    2013Apr09 DLLorenz Added setGD to plot
	#    2014May19 DLLorenz Modified criteria, and flag print
	#    2014Dec22 DLLorenz Roxygen header
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
  cvdfit <- qgrubbs(0.01, n)*sqrt(p/n)
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
