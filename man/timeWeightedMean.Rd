\name{timeWeightedMean}
\alias{timeWeightedMean}
\title{
Time-Weighted Mean
}
\description{
Compute the mean weighted by the duration between values.
}
\usage{
timeWeightedMean(x, time, na.rm = TRUE, fill = TRUE, by.period = FALSE, excessive = 2.5)
}
\arguments{
  \item{x}{
a sequnece of numeric values whose mean is to be computed.
}
  \item{time}{
the times associated with \code{x}. See \bold{Details}.
}
  \item{na.rm}{
remove missing values before computing the mean?
}
  \item{fill}{
expand the weights for the first and last observation to the end of the period? Otherwise,
the weight is based on the distance to the second or next-to-last observation.
}
  \item{by.period}{
computes wights by the periods defined by \code{time}.
}
  \item{excessive}{
a warning is issued if the largest weight for any observation exceeds \code{excessive} times the average weight.
}
}
\details{
The values of \code{time} are expected to be in decimal format, where the integer part
indicates the period and the fractional part linearly distributed through the period. For example,
the year 2000 begins at 2,000.0 and July 2, 2000 is 2,000.5. The function \code{dectime}
can be used to convert \code{Date} or \code{POSIX} to decimal format. If \code{time} is of
class "Date," or any "POSIX" classs, it is converted using \code{dectime}.\cr
If \code{na.rm} is \code{FALSE} and \code{byeriod} is \code{TRUE}, then missing values are
returned only for periods that contain a missing value (\code{NA}).
}
\value{
If \code{by.period} is \code{TRUE}, then a vector with one entry per period, otherwise the
mean for the entire data set.
}
\references{
Crawford,C.G., 2004, Sampling strategies for estimating acute and chronic exposures of 
pesticides in streams: Journal of the American Water Resources Association, v. 40, n. 2, p 485-502.
}
\seealso{
\code{\link{weighted.mean}}, \code{\link{dectime}}
}
\examples{
library(USGSwsData)
data(QW05078470)
with(QW05078470, timeWeightedMean(P00665, dectime(DATES)))
}
\keyword{univar}
