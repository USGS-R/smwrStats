\name{qtiles.CI}
\alias{qtiles.CI}
\title{
Quantiles with Confidence Limits
}
\description{
Compute sample quantiles and confidence limts for specified probabilities.
}
\usage{
qtiles.CI(x, probs = 0.5, CI = 0.9, bound=c("two.sided", "upper", "lower"),
  na.rm = TRUE)
}
\arguments{
  \item{x}{
numeric vector to compute the sample quantiles.
}
  \item{probs}{
numeric vector of desired probabilities with values between 0 and 1.
}
  \item{CI}{
the minimum desired confidence interval for each level specifed in probs. 
}
  \item{bound}{
a character string indicating the desired bounds, "two.sided" means the two-sided
interval, "upper" means the upper bound of the interval, and "lower" means the
lower bound of the interval. Only a single character is needed. The lower confidence
limit is \code{-Inf} when \code{bound} is "upper" and the upper confidence limit is
\code{Inf} when \code{bound} is "lower."
}
  \item{na.rm}{
logical; if TRUE, then missing values are removed before computation.
}
}
\value{
A matrix of sample quantiles, the lower confidence limit, the upper 
confidence limit, and the probability represented by the 
confidence interval corresponding to the probs levels in the sorted x data. 
}
\references{
Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in water resources: 
U.S. Geological Survey Techniques of Water-Resources Investigations, book 4, chap. A3, 522 p.\cr
}

\seealso{
\code{\link{quantile}}
}
\examples{
## Generate a random sample
set.seed(222)
XX.rn <- rexp(32)
qtiles.CI(XX.rn, probs=c(.25, .5, .75))
}

\keyword{univar}

