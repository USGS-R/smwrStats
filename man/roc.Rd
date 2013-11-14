\name{roc}
\alias{roc}
\title{
Receiver Operator Characteristics for Logistic Regression
}
\description{
Compute the receiver operator characteristics (ROC) for a logistic regression model.
}
\usage{
roc(object)
}
\arguments{
  \item{object}{
an object of class "glm."
}
}
\value{
An object of class "roc" haveing these components:
  \item{c.val}{the area under the ROC curve}
  \item{table}{a 3-column matrix of the fitted values, sensitivity and specificity}
}
\references{
Gonen, M., 2006, Receiver Operating Characteristics (ROC) Curves: SAS Users Group International
Conference, Paper 210-31 18 p., available at http://www2.sas.com/proceedings/sugi31/210-31.pdf, 
last accessed Octocer 26, 2011.
}

\note{
Objects of class "roc" have \code{print} and \code{plot} methods.
}

\keyword{regression}

