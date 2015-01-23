#' Maintenance of Variance Extension, Type 2
#' 
#' Calculates the Maintenance Of Variance Extension, Type 2 (MOVE.2) for record
#' extension by fitting a Line of Organic Correlation (LOC) regression model.
#' 
#' @details MOVE.2 has a necessarily predefined method for missing values---the
#'response variable is assumed to contain missing values and they are the values
#'to be estimated by the model equation. For the function \code{move.2}, missing 
#'values in the explanatory variable are excluded from the computations.
#'
#'If \code{distribution} is "normal," then the data in the explanatory variable 
#'and the response variable are assumed to have a bivariate normal distribution. 
#'Otherwise, they are assumed to have a bivariate log-normal distribution and a 
#'logarithmic ttransform is applied to both the explanatory variable and the 
#'response variable before analysis. The natural logarithm is used if 
#'\code{distribution} is "lognormal" and the commmon logarithm is used if 
#'\code{distribution} is "commonlog." Alternatively, the output from \code{optimBoxCox}
#'that contains both the response and explanatory variables can be supplied to
#'transform those variables by other than a logarithmic transform.
#' 
#' @param formula a formula object with the response variable on the left of
#'the ~ operator and a single explantory variable on the right.
#' @param data the data frame containing the variables named in \code{formula}.
#' @param subset an optional vector specifying a subset of observations to be
#'used in the fitting process.
#' @param distribution either "normal," "lognormal," "commonlog," or an object
#'of class "optimBoxCox" indicating the nature of the bivariate distribution, 
#'See \bold{Details}.
#' @param lag the number of days to account for tiem of travel between the explanatory
#'and response sites. If the explanatory site is upstream, then lag can be positive,
#'otherwise, lag can be negative to account for the travel time between the sites.
#' @return An object of class "move.2" having these components:
#'\item{coefficients}{ the intercept and slope of the line describing the fit.}
#'\item{R}{ Pearson's correlation coefficient.}
#'\item{p.value}{ the p-value from the correlation test, given a two-sided alternate hypothesis.}
#'\item{call}{ the matched call to \code{move.1}.}
#'\item{fitted.values}{ the fitted LOC values for the response.} 
#'\item{residuals}{ a 2-column matrix containing the signed distance from the 
#'predicted to to the corresponding the response variable and the explanatory variable values.}
#'\item{x}{ the (possibly transformed and lagged) values for the explanatory variable.} 
#'\item{y}{ the (possibly transformed) values for the response variable. }
#'\item{lag}{ the value of the \code{lag} argument. }
#'\item{xstats}{ the mean and standard deviation of the (possibly
#'transformed) values for the explanatory variable.}
#'\item{ystats}{ the mean and standard deviation of the (possibly 
#'transformed) values for the response variable.}
#'\item{var.names}{ the names of the response variable and the explanatory variable.} 
#'\item{model}{ the model frame.}
#'\item{data}{ the data frame supplied in \code{data}.}
#'\item{distribution}{ the value supplied in \code{distribution}.}
#' @note Objects of class "move.2" have \code{print}, \code{predict}, and
#' \code{plot} methods.
#' @seealso \code{\link{predict.move.2}}, \code{\link{plot.move.2}}, 
#' \code{\link{optimBoxCox}}
#' @references Hirsch, R.M., 1982, A comparison of four streamflow record
#'extension techniques: Water Resources Research, v. 18, p. 1081--1088.\cr
#'Moog, D.B., Whiting, P.J., and Thomas, R.B., 1999, Streamflow record extension using
#'power transformations and applicaitons to sediment transport: Water Resources
#'Research, v. 35, p 243--254.
#' @keywords models regression
#' @examples
#'\dontrun{
#'# See the vignette:
#'vignette("RecordExtension", package="smwrStats")
#'}
#' @export move.2
move.2 <- function(formula, data, subset, distribution="normal", lag=0) {
  ##
  ## Process formula as with regular linear regression!
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$distribution <- m$lag <- NULL
  m$drop.unused.levels <- TRUE
  m$na.action <- "na.pass"
  m[[1L]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  if(ncol(m) != 2L)
    stop("Must specify a single response and a single explanatory variable in the formula")
  Yall <- m[, 1L]
  yname <- names(m)[1L]
  Xall <- shiftData(m[, 2L], lag)
  xname <- names(m)[2L]
  # Extract the matching data and the extended x data (subscript 2 in Hirsch)
  good <- complete.cases(Xall,  Yall)
  X <- Xall[good]
  Y <- Yall[good]
  X2 <- Xall[!is.na(Xall) & !good]
  fit <- list(coefficients=double(2L))
  ## Extract missing value info
  fit$na.action <- attr(m, "na.action")
  ## Adjust for transforms
  if(inherits(distribution, "character")) {
    distribution <- match.arg(distribution, c("normal", "lognormal", "commonlog"))
    addTrans <- FALSE
    if(distribution == "lognormal") {
      Y <- log(Y)
      X <- log(X)
      X2 <- log(X2)
      xname <- paste("log(", xname, ")", sep='')
      yname <- paste("log(", yname, ")", sep='')
    } else if(distribution == "commonlog") {
      Y <- log10(Y)
      X <- log10(X)
      X2 <- log10(X2)
      xname <- paste("log10(", xname, ")", sep='')
      yname <- paste("log10(", yname, ")", sep='')
    }
    ## Otherwise normal
  } else { # Must be the output from optimBoxCox
    ## Transform according to boxCox
    ## Find the corresponding name
    xpick <- which(distribution$names == xname)
    ypick <- which(distribution$names == yname)
    xtrans <- round(distribution$lambda[xpick], 1L)
    xgm <- signif(distribution$gm[xpick], 4)
    xname <- paste("boxCox(", xname, ",", xtrans, ",", xgm, ")", sep='')
    X <- as.numeric(boxCox(X, xtrans, xgm))
    X2 <- as.numeric(boxCox(X2, xtrans, xgm))
    ytrans <- round(distribution$lambda[ypick], 1L)
    ygm <- signif(distribution$gm[ypick], 4)
    yname <- paste("boxCox(", yname, ",", ytrans, ",", ygm, ")", sep='')
    Y <- as.numeric(boxCox(Y, ytrans, ygm))
    addTrans <- TRUE
  }
  ybar <- mean(Y)
  yvar <- var(Y)
  xbar <- mean(X)
  xvar <- var(X)
  retcor <- cor.test(X, Y)
  fit$R <- retcor$estimate
  fit$p.value <- retcor$p.value
  # MOVE.2 code here
  N1 <- length(X)
  if(!(N2 <- length(X2))) {
    warning("There is no extended period for the explanatory variable (N2=0)!")
  } else {
    ybar1 <- ybar
    yvar1 <- yvar
    xbar1 <- xbar
    xvar1 <- xvar
    xbar2 <- mean(X2)
    xvar2 <- var(X2)
    xbar <- mean(c(X, X2))
    xvar <- var(c(X, X2))
    # Compute the adjusted mean of ybar
    ybar <- ybar1 + N2/(N1 + N2) * fit$R * sqrt(yvar1/xvar2) * (xbar2 - xbar1)
    alpha2 <- N2*(N1 - 4)*(N1 - 1)/((N2 - 1) * (N1 - 3) * (N1 -2))
    R2 <- fit$R^2
    yvar <- ((N1 - 1)*yvar1 + 
               (N2 - 1)*R2*yvar1/xvar1*xvar2 +
               (N2 - 1)*alpha2*(1 - R2)*yvar1 +
               N1*N2/(N1 + N2)*R2*yvar1/xvar1*(xbar2 - xbar1)^2)/(N1 + N2 - 1)
  }
  fit$coefficients[2L] <- sqrt(yvar/xvar)
  if(fit$R < 0)
    fit$coefficients[2L] <-  - fit$coefficients[2L]
  fit$coefficients[1L] <- ybar - fit$coefficients[2L] * xbar
  names(fit$coefficients) <- c("(Intercept)", xname)
  ## Complete steps for the model
  fit$call <- call
  fit$fitted.values <- cbind(1, X) %*% as.matrix(fit$coefficients)
  yresiduals <- Y - fit$fitted.values
  ## Now reverse the sense of the fit to get X-residuals
  junk <- double(2L)
  junk[2L] <- 1/fit$coefficients[2L]
  junk[1L] <- xbar - junk[2L] * ybar
  xfitted <- cbind(1, Y) %*% as.matrix(junk)
  xresiduals <- X - xfitted
  fit$residuals <- cbind(yresiduals, xresiduals)
  colnames(fit$residuals) <- c(yname, xname)
  ## Pack the object with other necessary information
  fit$x <- X
  fit$y <- Y
  fit$lag <- lag
  fit$xstats <- c(mean=xbar, sd=sqrt(xvar))
  fit$ystats <- c(mean=ybar, sd=sqrt(yvar))
  fit$var.names <- c(yname, xname)
  fit$N1 <- N1
  fit$N2 <- N2
  fit$model <- m
  fit$data <- data
  fit$distribution <- distribution # Required to evaluate if not character
  if(addTrans) {
    fit$xtrans <- xtrans
    fit$xgm <- xgm
    fit$ytrans <- ytrans
    fit$ygm <- ygm
  }
  ## set class to move.2
  oldClass(fit) <- "move.2"
  return(fit)
}

## The default extraction methods: coef, residuals, fitted work

