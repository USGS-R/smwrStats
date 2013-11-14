# Compute the Hosmer-Lemeshow test for logistic regression
#
# Coding History
#    2007May29 DLLorenz Initial version
#    2008Apr03 DLLorenz Begin bug fixes and tweaks
#    2011Aug22 DLLorenz Conversion to R
#    2011Oct25 DLLorenz Update for package
#    2012Feb03          This version.
#

hosmerLemeshow.test <- function(object, groups=10) {
  ## Arguments:
  ##  object (a glm model object) the logistic regression model
  ##  groups (integer scalar) the number of groups to use for the test.
  ##
  Dname <- as.character(object$call$formula)
  Dname <- paste(Dname[2], Dname[1], Dname[3]) # get it in the correct order
  Test <- fitted(object)
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
  ## This should return an object of class htest
  retval <- list(method="Hosmer-Lemeshow goodness of fit test",
                 statistic=Crit,
                 parameters=groups,
                 p.value=P.val,
                 data.name=Dname,
                 alternative="Some lack of fit\nnull hypothesis: No lack of fit",
                 estimate=cbind(Size=Sizes,
                   Expected=round(Expect, 3),
                   Counts=Observed))
  oldClass(retval) <- "htest"
  return(retval)
}
