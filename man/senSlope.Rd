\name{senSlope}
\alias{senSlope}
\title{
Compute the Sen Slope
}
\description{
Computes the Sen slope, confidence intervals, and an intercept for paired data.
}
\usage{
senSlope(formula, data, subset, na.action, intercept = "Ac", CI = 0.95)
}
\arguments{
  \item{formula}{
a model formula with exactly one explanatory variable.
}
  \item{data}{
the data.
}
  \item{subset}{
any descriptiopn to subset \code{data}.
}
  \item{na.action}{
the function to handle missing values.
}
 \item{intercept}{
a character string indicating the method to compute the intercept. See \bold{Details}.
}
  \item{CI}{
the desired confidence interval for the slope.
}
}
\details{
The argument \code{intercept} may be either "Ac" or "A1m." If it is "Ac," then
the intercept is computed from the median of \code{y} and \code{x}, also known
as the Conover method. If it is "A1m," then the intercept is chosen so that the
median of the residuals is zero.
}
\value{
An object of class "senSlope" with these components:
  \item{call}{the matched call.}
  \item{coefficients}{the intercept ans Sen slope.}
  \item{slope.CI}{the lower and upper confidence limits of the Sen slope.}
  \item{residuals}{the residuals of the regression.} 
  \item{fitted.values}{the fitted values.}
  \item{na.action}{information about any missing values.}
  \item{x}{the explanatory variable.}
  \item{y}{the response variable.}
  \item{var.names}{the response and explanatory variable names.}
  \item{model}{the model frame.}
}
\references{
Dietz, E.J., 1989, Teaching regression in a nonparametric statistics course:
The American Statstician, v. 43, p. 35--40\cr

Helsel, D.R., and Hirsch, R.M., 2002, Statistical methods in water resources: 
U.S. Geological Survey Techniques of Water-Resources Investigations, book 4, chap. A3, 522 p.\cr

Sen, P.K., 1968, Estimates of the regression coefficient based on Kendall's tau:
Journal of the American Statistical Association, v. 63 p. 1379--1389\cr
}

\seealso{
\code{\link{kensen.test}}, \code{\link{serial.test}}
}
\examples{
\dontrun{
library(USGSwsData)
data(SaddlePeaks)
senSlope(Flow ~ Year, data=SaddlePeaks)
}
}
\keyword{regression}
\keyword{nonparametric}
\keyword{robust}
