\name{print.ancovaReg}
\alias{print.ancovaReg}
\title{
Print Objects
}
\description{
Print the results of a analysis of covariance (\code{ancovaReg}).
}
\usage{
\method{print}{ancovaReg}(x, digits = 3, ...)
}
\arguments{
  \item{x}{
an object of class "ancovaReg" from \code{ancovaReg}.
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
The printed output contains the ANOVA table for the orignal models, the regression summary
for the final model, variance inflation factors for each explanatory variable in the final model,
and selected test criteria with observations that exceed one of more of the criteria.
}
