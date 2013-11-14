\name{vif}
\alias{vif}
\alias{vif.lm}
\title{
Variance Inflation Factors
}
\description{
Computes the variance inflation factor (Helsel and Hirsch, 2002) for each variable in a linear regression fit. 
}
\usage{
vif(model, ...)
\method{vif}{lm}(model, ...)
}
\arguments{
  \item{model}{
an object of class "lm." There is no default method for \code{vif}.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\value{
A named numeric vector containing the variance inflation factors for each variable.
}
\references{
Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in water resources: 
U.S. Geological Survey Techniques of Water-Resources Investigations, book 4, chap. A3, 522 p.\cr
}
\keyword{model}

