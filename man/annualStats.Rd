\name{annualStats}
\alias{annualStats}
\title{
Daily Value Statistics
}
\description{
Retrieve and compute annual summary statistics for USGS daily data.
}
\usage{
annualStats(gage, dtype = "sw", year = "water", begin.date = "", end.date = "", 
 param = NULL, complete = TRUE, statistic = mean, ...)
}
\arguments{
  \item{gage}{
a single USGS station identifier as a character string.
}
  \item{dtype}{
the data type, either "sw" for surface water, or "gw" for groundwater.
}
  \item{year}{
the kind of year, either "water" for water year, ending on September 30, or "calendar" for the calendar year.
}
  \item{begin.date}{
the date to use for the earliest value. If, "", then retrieve beginning with the first record in the database. 
}
  \item{end.date}{
the date to use for the latest value. If "", then retrieve the most recent values in the database.
}
  \item{param}{
the parameter code to retrieve. See \bold{Details}.
}
  \item{complete}{
logical: require complete years (no missing values) and remove provisional data? Otherwise, accept incomplete years.
}
  \item{statistic}{
the name of the statistical function to use. Must return a vector, may be longer than 1, like \code{quantile}.
}
  \item{\dots}{
additional arguments to \code{statistic}.
}
}
\details{
The value for \code{param} must match a valid 5-digit USGS parameter code. The default value for dtype = "sw" is "00060" 
(streamflow) and for dtype = "gw" is "72019" (water level below land surface).
}
\value{
A data frame contain columns for the agency code, station identifier, the year, the parameter code and the statistics.
}

\seealso{
\code{\link{readNWIS}}
}
\examples{
\dontrun{
# Get the annual mean flow for USGS streamgage 0533000
annualStats("05330000")
}
}

\keyword{IO}

