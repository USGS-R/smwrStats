\name{press}
\alias{press}
\title{
Prediction Error Sum of Squares
}
\description{
Compute the prediction error sum of squares statistic (PRESS) (Helsel and Hirsch, 2002) for a linear regression model. 
}
\usage{
press(model)
}
\arguments{
  \item{model}{
an object of class "lm" or the output from \code{lsfit}.
}
}
\value{
The prediction error sum of squares statistic.
}
\references{
Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in water resources: 
U.S. Geological Survey Techniques of Water-Resources Investigations, book 4, chap. A3, 522 p.\cr
}
\seealso{
\code{\link{multReg}}, \code{\link{lm}}
}
\keyword{regression}

