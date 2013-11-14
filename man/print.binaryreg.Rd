\name{print.binaryreg}
\alias{print.binaryreg}
\title{
Print Objects
}
\description{
Print the results of a logistic regression diagnotic analysis (\code{binaryReg}).
}
\usage{
\method{print}{binaryreg}(x, digits = 4, ...)
}
\arguments{
  \item{x}{
an object of class "binaryreg" from \code{mbinaryReg}.
}
  \item{digits}{
the number of significant digits to print numeric data.
}
  \item{\dots}{
not used for method, required for other methods.
}
}
\details{
The original regression model should be the output from \code{glm} with
\code{binomial} as the value for the \code{family} argument.
}
\value{
The object \code{x} is returned invisibly.{comp2 }{Description of 'comp2'}
}

\note{
The printed output contains the original call; a summary of the deviance residuals;
the coefficients with their estimates, standard errors, z scores, and attained p-values;
a summary of the deviance comparison between the model and null model (no explanatory varaibles),
also known as the G-squared or -2 log-likelihood; the response profile, which matches the value of the
response variable to 0 or 1 and the number of observations in each group; the le Cessie-van Houwelingen
and Hosmer-Lemeshow goodness of fit tests; three assessments of the predictinve power: 
the MaFadden R-squared and the adjusted R-squared, the classification
table, the concordance index, and the area under the reciever operating characteristics curve;
and selected test criteria with observations that exceed one of more of the criteria.
}

