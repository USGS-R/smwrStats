\name{plot.htest}
\alias{plot.htest}
\alias{plot.kensen}
\alias{plot.ppcc}
\title{Plot}
\usage{
\method{plot}{htest}(x, which = "All", set.up = TRUE, ...)

\method{plot}{kensen}(x, which = "All", set.up = TRUE, ...)

\method{plot}{ppcc}(x, which = "All", set.up = TRUE, ...)
}
\arguments{
  \item{x}{an object having classes "htest" and some other
  class for which a \code{plot} method exists.}

  \item{which}{either "All" or a number indicating which
  plot, see \bold{Details}.}

  \item{set.up}{set up the graphics page?}

  \item{\dots}{not used, required for other methods.}
}
\value{
The object \code{x} is returned invisibly.
}
\description{
Create diagnostic plots for selected hypothesis tests.
}
\details{
The \code{kensen} method for \code{plot} graphs the y and t
data with the best fit line.

The \code{ppcc} method creates a single graph that shows
the q-normal plot with a line showing the fit.
}
\keyword{hplot}

