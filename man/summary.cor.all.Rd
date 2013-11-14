\name{summary.cor.all}
\alias{summary.cor.all}
\title{
Summarize Correlations
}
\description{
Extract a data frame that summarizes the correlations computed by \code{cor.all}.
}
\usage{
\method{summary}{cor.all}(object, p.adjust = "holm", variable = NULL, ...)
}
\arguments{
  \item{object}{
an object created by the \code{cor.all} function.
}
  \item{p.adjust}{
a character string describing the method to use to adjust the p-values to account for m ultiple comparisons.
See \code{\link{p.adjust}} for the options and details.
}
  \item{variable}{
if the name of a variable, then summarize only those correlations of this variable with all others, otherwise
summarize all combinations.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\value{
A data frame containing columns of the paired variables, the correlation, the adjusted p-value, and the number of
observations in the correlation.
}
\seealso{
\code{\link{cor.all}}
}
\examples{
library(USGSwsData)
data(TNLoads)
# Extract only the correlations with the log of total nitrogen
summary(cor.all(TNLoads[, 1:5]), variable="LOGTN")
}
\keyword{htest}

