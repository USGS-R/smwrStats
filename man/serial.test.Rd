\name{serial.test}
\alias{serial.test}
\title{
Test for Serial Correlation
}
\description{
Perform either of two nonparametric tests (Wilcoxon Rank-Sum test or runs test) for serial correlation. 
}
\usage{
serial.test(x, method = "wilcoxon")
}
\arguments{
  \item{x}{
the numeric vector of observations. Missing values (NAs) are not allowed.
}
  \item{method}{
the string "wilcoxon" or "runs," depending on which test should be used. 
Only the first character is necessary. See \bold{Details}. 
}
}
\details{
If the method is "wilcoxon," then the Wilcoxon Rank-Sum test (Wilcoxon, 1945) 
is performed to compare the distribution of positive Z values to the distribution of negative Z values.\cr

If the method is "runs," then the runs test (Wald and Wolfowitz, 1940) is 
performed on sequences of positive and negative values of Z.\cr

Z is the product of sequential values of x after 
removing any trend and normalizing so that the median of Z values is 0.
}
\value{
An object of class htest containing the following components: 
  \item{statistic}{
the score for either test.
}
  \item{p.value}{
the two-sided p-value of the statistic.
}
  \item{Z}{
the raw data of the analysis. 
}
  \item{alternative}{
a string "two.sided" indicating the hypothesis.
}
  \item{method}{
a description of the method. 
}
  \item{data.name}{
a string containing the actual name of the input data. 
}
}
\references{
Dufour, J.M., 1981, Rank test for serial dependence: Journal of the Time Series Analysis, v. 2, no. 3, p. 117--128.\cr

Wald, A., and Wolfowitz, J., 1940, On a test whether two samples are from the same population: 
Annals of Mathematical Statistics, v. 11, p. 147--162.\cr

Wilcoxon, F., 1945, Individual comparisons by ranking methods: Biometrics, v. 1, p. 80--83.\cr

}

\note{
The null hypothesis is that the data are uncorrelated.
}

\keyword{nonparametric}
\keyword{htest}
