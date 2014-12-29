#' The Hosmer-Lemeshow Test
#' 
#' Perform the Hosmer-Lemeshow test for goodness-of-fit for a logistic
#' regression model.
#' 
#' 
#' @param object an object of class "glm" on which to perform the test.
#' @param groups the number of groups to use for the test.
#' @return An object of class "htest" having these components: \item{method}{ a
#' description of the method.  } \item{statistic}{ the test statistic.  }
#' \item{p.value}{ the attained p-level of the test statistic.  }
#' \item{data.name}{ the name of \code{object}.  } \item{alternative}{ the
#' alternate hypothesis---"some lack of fit."  } \item{estimate}{ a data frame
#' of the size, expected value, and actual counts in each group. If the model
#' has a single explanatory variable, then the mean value is included as column
#' 4.  }
#' @note The null hypothesis is "no lack of fit." Rejection of the null
#' hypothesis indicates "some lack of fit."
#' @seealso \code{\link{binaryReg}}
#' @references Hosmer, D.W. and Lemeshow, S., 1980, Goodness-of-fit tests for
#' the multiple logistic regression model: Communications in Statistics ---
#' Theory and Methods, v. 9, p. 1043--1069
#' @keywords htest
#' @export hosmerLemeshow.test
hosmerLemeshow.test <- function(object, groups=10) {
	# Coding History
	#    2007May29 DLLorenz Initial version
	#    2008Apr03 DLLorenz Begin bug fixes and tweaks
	#    2011Aug22 DLLorenz Conversion to R
	#    2011Oct25 DLLorenz Update for package
	#    2014Dec22 DLLorenz Roxygen headers
  ##
	ExpClasses <- attributes(terms(object))$dataClasses[-1L]
	if(length(ExpClasses) == 1L) {
		if(!(ExpClasses %in% c("numeric", "integer")))
			stop("Hosmer-Lemeshow test not appropriate for a single, non-numeric explanatory variable.")
		ExpVar <- object$model[, 2L]
	}
  Dname <- as.character(object$call$formula)
  Dname <- paste(Dname[2], Dname[1], Dname[3]) # get it in the correct order
  Test <- na.omit(fitted(object))
  Groups <- cut(rank(Test), groups, factor.result=TRUE,
                labels=paste("group", seq(groups)))
  Sizes <- table(Groups)
  Resp <- object$y
  Expect <- tapply(Test, Groups, mean) * Sizes
  Observed <- tapply(Resp, Groups, sum)
  Crit <- sum((Observed - Expect)^2/(Expect*(1 - Expect/Sizes)))
  names(Crit) <- "Chi-square"
  P.val <- 1 - pchisq(Crit, groups-2)
  names(groups) <- "Number of groups"
  ## Return an object of class htest
	if(length(ExpClasses) > 1L) {
		retval <- list(method="Hosmer-Lemeshow goodness of fit test",
									 statistic=Crit,
									 parameters=groups,
									 p.value=P.val,
									 data.name=Dname,
									 alternative="Some lack of fit\nnull hypothesis: No lack of fit",
									 estimate=data.frame(Size=as.vector(Sizes),
									 							 Expected=as.vector(round(Expect, 3)),
									 							 Counts=as.vector(Observed)))
	} else {
		ExpVar <- tapply(ExpVar, Groups, mean)
		estimate=data.frame(Size=as.vector(Sizes),
												Expected=as.vector(round(Expect, 3)),
												Counts=as.vector(Observed), 
												V4=as.vector(ExpVar))
		names(estimate)[4L] <- names(ExpClasses)
		retval <- list(method="Hosmer-Lemeshow goodness of fit test",
									 statistic=Crit,
									 parameters=groups,
									 p.value=P.val,
									 data.name=Dname,
									 alternative="Some lack of fit\nnull hypothesis: No lack of fit",
									 estimate=estimate)
	}
	oldClass(retval) <- "htest"
	return(retval)
}
