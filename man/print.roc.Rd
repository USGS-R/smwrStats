\name{print.roc}
\alias{print.roc}
\title{
Print Objects
}
\description{
Print the results of a reciever operator characteristics analysis (\code{roc}).
}
\usage{
\method{print}{roc}(x, digits = 3, ...)
}
\arguments{
  \item{x}{
an object of class "roc" from \code{roc}.
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
The area under the reciever operator characteristics (ROC) curve is printed.
}
