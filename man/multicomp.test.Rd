\name{multicomp.test}
\alias{multicomp.test}
\title{
Multiple Comparisons
}
\description{
Perform mutliple comparison tests among groups of data. the tests may
be either parametric, nonparametric, or Dunn's nonparametric (Glantz, 2005).
}
\usage{
multicomp.test(x, g, method="parametric", critical.value="", alpha=0.05)
}
\arguments{
  \item{x}{
the numeric vector of observations. Missing values (NAs) are allowed and removed before the test is performed.
}
  \item{g}{
any group vector for the observations. Missing values (NAs) are allowed and removed before the test is performed.
}
  \item{method}{
a character string describing the test.
Only the first character is necessary. See \bold{Details}. 
}
  \item{critical.value}{
a character string describing the method to use for determining the critical value.
Only the first character is necessary. See \bold{Details}. 
}
  \item{alpha}{
the significance level of the test. See \bold{Note}. 
}

}
\details{
If the method is "parametric," then the comparisons are based on the means and variances of the raw data.
 The valid choices for \code{critical.value} are "tukey" (default), "bonferroni," or "lsd."\cr
Otherwise, the comparisons are based on the ranks of the data. Valid choices for \code{critical.value} are 
 "tukey" (default), "bonferroni," or "lsd" when \code{method} is "nonparametric" and "sidak" (default)
 or "bonferroni" when \code{method} is "dunn."\cr
The basic diffference between the default nonparametric method and Dunn's nonparametric method is the handling of
ties.
}
\value{
An object of class MCT containing the following components: 
  \item{title}{
a description of the test.
}
  \item{cv.method}{
the method used to compute the critical value.
}
  \item{alpha}{
the value of \code{alpha}. 
}
  \item{crit.value}{
the critical value for the pairwise comparisons.
}
  \item{response}{
the name of the response variable. 
}
  \item{groups}{
the name of the group variable. 
}
  \item{means}{
the means for each group.
}
  \item{sizes}{
the number of observations in each group.
}
  \item{table}{
the table of the results of the pairwise comparisons.
}
 \item{assoc}{
a data frame containing the possible association for each group.
}
}
\note{
All computations of the variance for unequal group sizes are based on the harmonic mean as
 described in Yandell (1997). That adjustment is only approximate when \code{critical.value}
 is "tukey" and \code{method} is "parametric" but useful when the design is slightly unbalanced.\cr
The default nonparametric method \code{method} = "nonparametric" is only assymptotically unbiased when some
 data are tied. For smaller data sets with small numbers of ties, it may be preferable to use Dunn's
 nonparametric method \code{method} = "dunn."
}
\references{
Glantz, S.A., 2005, Primer of biostatistics: McGraw Hill, New York, 520 p.\cr
Helsel, D.R., and Hirsch, R.M., 2002, Statistical methods in water resources:
 U.S. Geological Survey Techniques of Water-Resources Investigations, book 4, chap. A3, 522 p. 
Higgins, J.J., 2004, Introduction to modern nonparametric statistics:
 Pacific Grove, Calif., Brooks/Cole, 384 p.\cr
Yandell, B.S., 1997, Practical data analysis for designed experiments:
 London, United Kingdom, Chapman & Hall, 437 p.\cr
}
\keyword{nonparametric}
\keyword{htest}
