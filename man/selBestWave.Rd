\name{selBestWave}
\alias{selBestWave}
\title{
Select the "Best" Seasonal Wave
}
\description{
Select the "best" parameters for a seasonal wave fit.
}
\usage{
selBestWave(formula, data, dec.time, wave.list, exhaustive = FALSE, Regression = lm, Test = AIC, ...)
}
\arguments{
  \item{formula}{
the formula describing the model without a \code{seasonalWave} term.
}
  \item{data}{
the data.frame that contians the veriable specified in the formula.
}
  \item{dec.time}{
a character string of the name of the column in \code{data} in decimal time format.
}
  \item{wave.list}{
an object of class "seasonalPeak" (confirmed) describing the timing of the peak and potential candidate models.
}
  \item{exhaustive}{
logical; if TRUE, then do a fairly complete search for the timing of the peak, otherwise
accept the timing specified in \code{wave.list}.
}
  \item{Regression}{
the regression function.
}
  \item{Test}{
the function to perform the comparison test among all of the candidate models.
}
  \item{\dots}{
any additional arguments to \code{Regression}.
}
}
\details{
For logistic regression, use \code{Regression}=\code{glm}, \code{Test}=\code{deviance},
\code{family=binomial(link="logit")}.
}
\value{
A 4- or 7-column matrix of the "best" models. The columns are the timing of the peak (\code{Cmax}),
the primary peak loading (\code{Loading}), the half life (\code{Hlife}), and the 
test score (\code{Test}). If the model is a two-peak model, then additional columns
(\code{la}, \code{lo}, and \code{w}) describeing the second peak are included.
}
\references{
Vecchia, A.V., Martin, J.D., and Gilliom, R.J., 2008, Modeling variability and trends in 
pesticide concentrations in streams: Journal of the American Water Resources Association, 
v. 44, no. 5, p.1308--1324.
}
\examples{
## See the SeasonalWave demo
}
\keyword{models}
\keyword{regression}
