\name{seasonalWave.wt}
\alias{seasonalWave.wt}
\title{
Compute Seasonal Wave Model
}
\description{
This is a support function for seasonalWave model (Vecchia and others, 2008).
}
\usage{
seasonalWave.wt(loading, second.peak)
}
\arguments{
  \item{loading}{
the number of months of primary loading for the seasonal wave model.
}
  \item{second.peak}{
a list describing the loading characteristics of the second peak. See \bold{Details}. 
If \code{NULL}, then no second peak is added to the seasonal wave model.
}
}
\details{
If \code{second.peak} is supplied, then it must be a list with these components:\cr
la, the lag (number of months) between the primary and second peak;\cr
lo, the loading duration (number of months) for the second peak (must be less than \code{la};\cr
w, the weigthing factor for the second peak loading, relative to the primary peak. Must be
greater than 0 and less than or equal to 1.
}
\value{
The weighting vector for the seasaon wave model. A vector of length 12 that represents the
relative loading to the system.
}
\references{
Vecchia, A.V., Martin, J.D., and Gilliom, R.J., 2008, Modeling variability and trends in 
pesticide concentrations in streams: Journal of the American Water Resources Association, 
v. 44, no. 5, p. 1308-1324.
}
\note{
This function is support for the seasonalWave function and is not intended to be called by the user.
}
\seealso{
\code{\link{seasonalWave}}
}
\keyword{manip}
