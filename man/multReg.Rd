\name{multReg}
\alias{multReg}
\title{
Diagnostics for Linear Regression
}
\description{
Computes diagnostics for linear regression
}
\usage{
multReg(object)
}
\arguments{
  \item{object}{
the linear regression model object
}
}
\value{
A object of class "multReg" having the following components:\cr

 \item{aovtab}{
the analysis of variance table, using the type III sum of squares.
}
 \item{parmests}{
a summary of \code{object}. 
}
 \item{vif}{
a named vector of variance inflation factors.
}
 \item{diagstats}{
a data.frame containing the observed values, predicted values, residuals, 
standardized residuals, studentized residuals, leverage, Cook's D, and dfits for each observation.
}
 \item{crit.val}{
a named vector of the critical values for leverage, Cook's D, and dfits.
}
 \item{flagobs}{
a logical vector indicating which observations exceeded at least one of the critical values.
}
 \item{object}{
the \code{lm} object.
}
 \item{x }{
the model matrix of explanatory variables.
}
}
\references{
Draper, N.R. and Smith, H., 1998, Applied Regression Analysis, (3rd ed.): New York, Wiley, 724 p.\cr
}

\note{
The type II sum of squares are calculated according to the principle of marginality, testing each term 
after all others, except ignoring the term's higher-order relatives. This type sum of squares is useful
for assessing the overall marginal effect of each term in the model.\cr
Objects of class "multReg" have \code{print} and \code{plot} methods.
}

\seealso{
\code{\link{lm}}, \code{\link{plot.multReg}},
}
\keyword{models}
\keyword{regression}
