#' Multivariate Unconditional Box-Cox Transformations
#'
#' Compute Box-Cox transformations that maximize the log likelihood
#'of the transformations.
#'
#' @param X a data frame or matrix of the data to find the optimized Box-Cox
#'transforms to produce multivariate normality. Can also be a numeric vector for 
#'a simple Box-Cox transform to normality.
#' @param start a numeric vector of length matching the number of columns in \code{X}
#'to provide starting values for the Box-Cox transforms.
#' @note The maximum likelihood estimate of the Box-Cox transformations corresponds to the
#'minimum determinant of the variance-covariance matrix of the transformed \code{X}. The
#'methodology is described in Andrews and others (1971).
#' @return An object of class "optimBoxCox" having these components:
#'\item{start}{ the starting values for the Box-Cox transformations.}
#'\item{criterion}{ the log-likelihood of the Box-Cox transformations.}
#'\item{names}{ the names of the columns.}
#'\item{lambda}{ the values of the Box-Cox transformations.}
#'\item{stderr}{ the standard errors of the Box-Cox transformations.}
#'\item{return.code}{ the convergence value returned by \code{optim}.}
#'\item{gm}{ the geometric means of the data in \code{X}.}
#'\item{data}{ the data in \code{X} with missing values removed.}
#' @seealso \code{\link[smwrBase]{boxCox}}, \code{\link[stats]{optim}}
#' @references
#'Andrews, D.F., Gnanadesikan, R., and Warner, J.L., 1971, Transformations
#'of multivariate data: Biometrics, v. 27, p. 825--840.
#' @export
optimBoxCox<-function(X, start=NULL) {
  ## The function to compute the negative of the log likelihood of the transformation
  negLogL<-function(X, lambda, gm) {
    for (j in 1:ncol(X)){
      X[,j]<-boxCox(X[,j],lambda[j],gm[j])
    }
    (nrow(X)/2)*log(((nrow(X)-1)/nrow(X))*det(var(X)))
  }
  ## Begin
  X <- na.omit(data.matrix(X))
  nc <- ncol(X)
  if(any(X <= 0)) stop("All values must be > 0")
  gm <- apply(X, 2, function(x) exp(mean(log(x))))
  ## This is a reasonable starting value for lambda
  if (is.null(start)) {
    start <- double(nc)
    for (j in 1:nc) {
      start[j] <- 1 - skew(X[,j])/2
    }
  }
  ## Different arguments for single/multivariate
  if(nc == 1L) {
    method <- "Brent"
	lower <- -100
	upper <- 100
  } else {
	method <- "Nelder-Mead"
	lower <- -Inf
	upper <- Inf
  }
  ## Ready, set, go.
  ret<-optim(start, negLogL, hessian=TRUE, method=method, X=X, gm=gm, 
             lower=lower, upper=upper)
  retval <- list(start=start, criterion=ret$value, names=colnames(X),
                 lambda=ret$par, stderr=sqrt(diag(solve(ret$hessian))),
                 return.code=ret$convergence, gm=gm, data=X)
  if(retval$return.code != 0) 
        warning(paste("Convergence failure: return code =",
                      retval$return.code))
  # Convert to log-likelihood
  retval$criterion <- -retval$criterion
  class(retval)<-"optimBoxCox"
  retval
}

