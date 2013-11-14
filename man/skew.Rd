\name{skew}
\alias{skew}
\title{
Skewness
}
\description{
Compute the skewness statistic.
}
\usage{
skew(x, na.rm = TRUE, method = "fisher")
}
\arguments{
  \item{x}{
any numeric vector.
}
  \item{na.rm}{
logical; if TRUE, then remove missing value before computation.
}
  \item{method}{
  the method to use for computing the skew, must be "fisher" or "moments."
}
}
\value{
a single value representing the skewness of the data in \code{x}.
}
\references{
Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in water resources: 
U.S. Geological Survey Techniques of Water-Resources Investigations, book 4, chap. A3, 522 p.\cr
}

\examples{
skew(c(1.0, 1.2, 1.5, 1.9, 2.5))
}
\keyword{univariate}

