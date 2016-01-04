#' The le Cessie-van Houwelingen Test
#' 
#' Performs the le Cressie-van Houwelingen test for goodness-of-fit for a
#' logistic regression model.
#' 
#' If \code{bandwidth} is missing, then the mean distance between observations
#' is used.
#' 
#' @param object an object of class "glm" on which to perform the test.
#' @param bandwidth the bandwidth for smoothing the residuals.
#' @param newterms any new variables to add to the model. Expressed as the
#' right-hand side of a formula.
#' @return An object of class "lecessie" having these components:
#' \item{method}{ a description of the method.  } \item{statistic}{ the test
#' statistic.  } \item{parameters}{ the degrees of freedom of the chi-squared
#' test.  } \item{p.value}{ the attained p-level of the test statistic.  }
#' \item{data.name}{ the name of \code{object}.  } \item{alternative}{ the
#' alternate hypothesis--"some lack of fit."  } \item{estimate}{ the estimated
#' values for the test.  } \item{object}{ the original \code{object}.  }
#' \item{target.object}{ the \code{object} with any added \code{newterms}.  }
#' \item{bandwidth}{ the bandwidth used for smoothing the residuals.  }
#' \item{max.distance}{ the maximum distance between observations.  }
#' \item{smoothed.residuals}{ the smoothed residuals.  }
#' \item{distance.matrix}{ a matrix of the distances between observations.  }
#' \item{hat}{ the hat matrix.  }
#' @note The null hypothesis is "no lack of fit." Rejection of the null
#' hypothesis indicates "some lack of fit."
#' @seealso \code{\link{binaryReg}}
#' @references le Cessie, S. and van Houwelingen, H.C., 1995, Testing the fit
#' of a regression model via score tests in random effects models: Biometrics,
#' v. 51, p 600-614.
#' @keywords htest
#' @export leCessie.test
leCessie.test <- function(object, bandwidth, newterms) {
	# Coding history:
	#    2009Jan25 DLLorenz Original Coding and begin of modifications
	#    2011Aug22 DLLorenz Conversion to R
	#    2011Oct25 DLLorenz Update for package
	#    2012Aug28 DLLorenz Renamed to correct spelling
	#    2012Aug28 DLLorenz Change from diagPlot to plot
	#    2013Apr09 DLLorenz Added setGD to plot
	#    2013Sep26 DLLorenz Slight improvement to speed for large models
	#    2014Dec22 DLLorenz Roxygen headers
  ##
  ## Note there are very minor differences between this version and the SAS
  ## code by Saskia le Cessie.
  ## The errors begin with very small differences in the computed distances
  ## and propogate from there.
  ##
  ## Create the categorical distance function
  categorical <- function(x) {
    x <- as.matrix(x)
    retval <- apply(x, 2, function(y) {
      nrcats <- length(unique(y))
      ## allow for the case where the number of categories is 1
      if(nrcats == 1)
        disty <- rep(0, length(y)*(length(y)-1))
      else
        disty <- (dist(y, 'manhattan') != 0) * nrcats/(nrcats - 1)
    })
    retval <- rowSums(retval)
    return(retval)
  }
  ## Keep the original object first 
  object.orig <- object
  Dname <- as.character(object$call$formula)
  Dname <- paste(Dname[2], Dname[1], Dname[3]) # get it in the correct order
  ## Extract the components and compute the smoothed residuals and raw Q
  fits <- fitted(object)
  resids <- resid(object, 'response')
  N <- length(resids)
  if(!missing(newterms))
    object <- update(object, formula=newterms)
  obj.frame <- model.frame(object)
  if(class(obj.frame[[1]]) == 'matrix')
    stop("Cannot perform test where response is matrix")
  ## Compute distances
  distances <- lapply(obj.frame[,-1, drop=FALSE], function(x) {
    if(is.numeric(x))
      return(unclass(0.5*(dist(scale(x)))^2))
    if(class(x) == 'factor')
      return(unclass(categorical(x)))
    ## Must be a structure created from a function that creates a matrix
    return(unclass(0.5*dist(scale(x[[1]]))^2))
  })
  dist.mat <- matrix(0, N, N)
  dist.mat[lower.tri(dist.mat)] <- sqrt(rowSums(as.data.frame(distances)))
  dist.mat <- dist.mat + t(dist.mat)
  if(missing(bandwidth))
    bandwidth <- mean(dist.mat)
  R.raw <- pmax(1 - dist.mat/bandwidth, 0)
  smoothres <- t(resids) %*% R.raw
  Q.raw <- sum(smoothres * resids)
  ## Extract components and compute corrected estimated value of Q
  obj.mat <- model.matrix(object.orig)
  mu2 <- fits * (1 - fits)
  V <- diag(mu2)
  Vx <- diag(V)
  obj.hat <- Vx*obj.mat %*% solve(t(obj.mat) %*% (Vx * obj.mat)) %*% t(obj.mat)
  R.cor <- (diag(N) - obj.hat) # set up to re-use below
  R.cor <- R.cor %*% R.raw %*% R.cor
  E.Q <- sum(diag(R.cor) * mu2)
  ## Extract components and compute corrected std.err of Q
  mu4 <- mu2*(1-3*mu2)
  VarQ1 <- sum(diag(R.cor)^2*(mu4 - 3*mu2^2))
  R.tmp <- R.cor * rep(mu2, each=nrow(R.cor))
  VarQ2 <- 2*sum(diag(R.tmp %*% R.tmp))
  VarQ <- VarQ1 + VarQ2
  Test <- Q.raw * 2 * E.Q / VarQ
  df <- 2 * E.Q^2 / VarQ
  ## make the return value able to be printed with print.htest
  names(Test) <- "Chisq"
  names(df) <- "df"
  retval <- list(statistic=Test, parameters=df,
                 p.value=1 - pchisq(Test, df),
                 estimate=c(Q=Q.raw, 'E[Q]'=E.Q, 'se[Q]'=sqrt(VarQ)),
                 method="le Cessie-van Houwelingen GOF test",
                 data.name=Dname,
                 alternative="Some lack of fit\nnull hypothesis: No lack of fit",
                 ## end of htest material
                 object=object.orig, target.object=object,
                 bandwidth=bandwidth, max.distance=max(dist.mat),
                 smoothed.residuals=as.vector(smoothres),
                 distance.matrix=dist.mat,
                 hat=obj.hat)
  oldClass(retval) <- "lecessie"
  return(retval)
}
