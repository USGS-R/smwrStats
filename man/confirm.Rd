\name{confirm}
\alias{confirm}
\alias{confirm.default}
\title{
Confirm an Analysis
}
\description{
Review and accept the results of an analysis
}
\usage{
confirm(x, ...)
\method{confirm}{default}(x, ...)
}
\arguments{
  \item{x}{
the object to be confirmed
  }
  \item{\dots}{
additional arguments required for specific methods
  }
}
\value{
The object returned depends on the specific method.
}
\note{
The default method simply returns the object and issues a warning.
}
\keyword{manip}
