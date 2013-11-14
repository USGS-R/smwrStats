\name{Grubbs}
\alias{pgrubbs}
\alias{qgrubbs}
\title{
The Grubbs Distribution
}
\description{
The one-sided critical value and the attained signicance level for the Grubbs test for an outlier in a sample from Normal distribution.
}
\usage{
pgrubbs(G, N)
qgrubbs(alpha, N)
}
\arguments{
  \item{G}{
the maximum or minimum scaled difference from the mean.
}
  \item{N}{
the number of values in the sample. 
}
  \item{alpha}{
the significance level.
}
}
\value{
The function \code{pgrubbs} returns the attained p-value, not the cumulative distribution, for the maximum scaled distance from the mean.
the funciton \code{qgrubbs} returns the one-sided critical value for the maximum scaled difference from the mean.
}
\references{
Grubbs, F., 1969, Procedures for Detecting Outlying Observations in Samples, Technometrics, v. 11, no. 1, pp. 1-21. 
}
\seealso{
\code{\link{grubbs.test}}
}
\examples{
# The difference is due to rounding errors
pgrubbs(c(.9, .95, .99), 32)
qgrubbs(c(5.905348, 5.483234, 5.159097), 32)
}
\keyword{distribution}

