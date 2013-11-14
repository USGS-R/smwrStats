\name{print.lecessie}
\alias{print.lecessie}
\title{
Print Objects
}
\description{
Print the results of a le Cessie-van Houwelingen test (\code{leCessie.test}).
}
\usage{
\method{print}{lecessie}(x, digits = 4, ...)
}
\arguments{
  \item{x}{
an object of class "lecessie" from \code{leCessie.test}.
}
  \item{digits}{
the number of significant digits to print numeric data.
}
  \item{\dots}{
not used for method, required for other methods.
}
}
\value{
The object \code{x} is returned invisibly.
}
\note{
The printed ouput is very similar to the printed ouput for class "htest," the output
from a hypothesis test, but includes the additional output summarizing the distance
between observations, which can be useful for comparing the reults with different
bamdwidth settings.
}
