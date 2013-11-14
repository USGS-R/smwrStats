\name{trends}
\alias{trends}
\title{
Basis for Piecewise Linear Trends
}
\description{
Generate a basis matrix for piecewise linear modeling of trends.
}
\usage{
trends(x, breaks, boundary.breaks = c(FALSE, FALSE), steps)
}
\arguments{
  \item{x}{
a vector of dates/times, assumed to be in dectime format. Missing values are
permitted and result in corresponding missing values in the output.
}
  \item{breaks}{
a vector of breakpoints in the linear trends.
}
  \item{boundary.breaks}{
a logical vector of length 2, indicating whether the breaks include the range of the data. 
If the first is TRUE, then the first break denotes the beginning of a trend. 
If the last is TRUE, then the last break denotes the end of a trend. Otherwise, the trends
begin at the \code{floor} of the first value of \code{x} and end at the \code{ceiling} of 
the last value of \code{x}.
}
  \item{steps}{
a vector indicating any step trends.
}
}

\value{
A matrix of dimension length of x by number of linear trend and step trends. 
The breaks are included as an attribute.
}

\note{
Each trend is 0 prior to the break, and then increases at the rate of 1 per unit and 
maintain that maximum value after the break. The regression coefficient then reprents the
trend as a rate. A step trend is 0 before the step and 1 after it. 
}

\seealso{
\code{\link{floor}}, \code{\link{ceiling}}, \code{\link{curvi}}
}
\examples{
# model two piecewise linear trends from 2000 to 2004, with a break at 2001 and 2003
trends(2000 + seq(0,20)/5, breaks=c(2001, 2003))
}
\keyword{model}

