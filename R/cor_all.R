#' Compute Cross Correlations
#' 
#' Compute correlations among numeric data
#' 
#' The null hypothesis is that the data are not correlated with one another.
#' The alternate hypothesis is that they are correlated with one another. This
#' is a two-sided test. For other options, see \code{cor.all}.
#' 
#' If \code{method} is "pearson," then the correlation is based on Pearson's
#' product moment correlation coefficient. If \code{method} is "kendall," then
#' Kendall's tau is used to estimate a rank-based measure of association. If
#' \code{method} is "spearman", then Spearman's rho is used to estimate the
#' correlation of the ranks of the data. The last two methods may be used if
#' the data do not necessarily come from a bivariate normal distribution.
#' 
#' If \code{na.method} is "fail," then \code{cor.all} stops if there are any
#' missing numeric values. If it is "omit," then all rows with any missing
#' values is removed before the correlations are computed. That option will
#' always produce a correlation matrix that is positive definite. If
#' \code{na.method} is "pairwise," then missing values are removed from each
#' pairwise correlation.
#' 
#' If \code{distribution} is "normal," then the assumption for \code{method} =
#' "pearson" is that the data are bivariate normal. If \code{distribution} is
#' "lognormal," then the assumption for \code{method} = "pearson" is that the
#' data are bivariate log-normal and all data are natural log-transformed. If
#' \code{distribution} is "log1p," then the assumption for \code{method} =
#' "pearson" is that the data are bivariate log-normal after adding 1 and all
#' data are transformed using the \code{log1p} function. The data are
#' transformed for any \code{method}, but only produce a different result for
#' \code{method} = "pearson."
#' 
#' @param data any rectangular object such as a data.frame or matrix.
#' @param method a character string indicating which correlation coefficient is
#' to be used. One of "pearson," "kendall," or "spearman," can be abbreviated.
#' @param na.method a character string indicating which method to use for
#' missing values. One of "fail," "omit," "pairwise," can be abbreviated.
#' @param distribution a character string indicating the assumed distribution
#' of the data. One of "normal," "lognormal," or "log1p", which can be
#' abbreviated.
#' @return An object of class "cor.all," which has these components:
#' \item{estimates}{ a matrix of the correlations between each pair of numeric
#' variables in \code{data} } \item{p.values}{ a matrix of the attained
#' p-values between each pair of numeric variables in \code{data} }
#' \item{counts}{ a matrix of observations in each pair of numeric variables in
#' \code{data} } \item{alternative}{ a character string indicating the
#' alternative hypothesis, always "two.sided" } \item{na.method}{ a character
#' string indicating the method to handle missing values } \item{method}{ a
#' character string describing the method to compute the correlations }
#' \item{data.name}{ the name of the data set, derived from \code{data} }
#' \item{data}{ a data frame of the numeric variables } \item{call.method}{ a
#' character string indicating the method to compute the correlations }
#' \item{distribution}{ a character string indicating the distribution
#' assumption of the data }
#' @note The \code{print}, \code{plot}, and \code{summary} methods are
#' available for an object of class "cor.all."
#' @seealso \code{\link{cor.test}}, \code{\link{plot.cor.all}},
#' \code{\link{summary.cor.all}}
#' @references Conover, W.J., 1980, Practical nonparametric statistics (2d
#' ed.): New York, Wiley, 512 p.\cr
#' 
#' Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in water resources:
#' U.S.  Geological Survey Techniques of Water-Resources Investigations, book
#' 4, chap. A3, 522 p.
#' @keywords htest nonparametric
#' @examples
#' 
#' library(smwrData)
#' data(TNLoads)
#' cor.all(TNLoads[, 1:5])
#' cor.all(TNLoads, method="spearman")
#' 
#' @export cor.all
cor.all <- function(data, method="pearson", na.method="pairwise",
                    distribution="normal") {
	# Coding history:
	#    2004Oct14 DLLorenz Initial coding for library.
	#    2005Jul14 DLLorenz date fix
	#    2006Feb03 DLLorenz Bug Fix when not all data were numeric
	#    2008Apr21 DLLorenz Continuity Adjustment to Kendall's method
	#    2011Jul28 DLLorenz Conversion to R
	#    2011Oct25 DLLorenz Update for package
	#    2012Mar08 DLLorenz Added summary function
	#    2012Jun19 DLLorenz Tweeks in summary and other methods
	#    2012Aug28 DLLorenz Change from diagPlot to plot
	#    2012Dec21 DLLorenz Change gaussian to normal
	#    2013Apr01 DLLorenz Bug fix in printed output
	#    2013Apr09 DLLorenz Added setGD to plot
	#    2013Apr31 DLLorenz Switched which to y-axis and bug fix
	#    2014May27 DLLorenz Added log1p distribution
	#    2014Dec22 DLLorenz Roxygen headers
  ##
  ## Compute correlations between all numeric columns in data
  data.name <- deparse(substitute(data))
  data <- as.data.frame(data) # force to be a data. frame
  nums <- unlist(lapply(data, is.numeric))
  nams <- names(data)[nums]
  nums <- which(nums)
  N <- length(nams)
  if(N < 2)
    stop("Insufficient numeric data to compute correlations")
  ## Set up for handling NAs
  na.methods <- c("fail", "omit", "pairwise")
  na.meth <- pmatch(na.method, na.methods, nomatch=0)
  data <- data[, nums] # extract the data to be analyzed
  if(na.meth == 1) { # check if any NAs
    if(any(is.na(unlist(data))))
      stop("There are missing values in the numeric data")
  }
  else if(na.meth == 2)
    data <- na.omit(data)
  else if(na.meth == 0)
    stop("cor.all: Invalid na.method specified")
  ## Apply log-transform if requested
  distribution <- match.arg(distribution, c("normal", "lognormal", "log1p"))
  if(distribution == "lognormal") {
  	data <- as.data.frame(lapply(data, log))
  } else if(distribution == "log1p") 
  	data <- as.data.frame(lapply(data, log1p))
  cor.grid <- expand.grid(1:N, 1:N)
  cor.grid <- cor.grid[cor.grid[,1] < cor.grid[,2],]
  cor.stat <- matrix(1, ncol=N, nrow=N)
  cor.pval <- matrix(0, ncol=N, nrow=N)
  ## Construct matrix of counts of samples
  cor.count <-as.matrix( !is.na(data))
  cor.count <- t(cor.count) %*% cor.count
  ## Protect against warning messages for spearman and kendall options
  warn <- options("warn")
  options(warn=-1)
  method=match.arg(method, c("pearson", "spearman", "kendall"))
  for(i in 1:nrow(cor.grid)) {
    Xs <- cor.grid[i,1]
    Ys <- cor.grid[i,2]
    temp <- try(cor.test(data[[Xs]], data[[Ys]], alternative="two.sided",
                         method=method, continuity=TRUE))
    if(class(temp) != "Error") {
      cor.pval[Xs, Ys] <- cor.pval[Ys, Xs] <- temp$p.value
      cor.stat[Xs, Ys] <- cor.stat[Ys, Xs] <- temp$estimate
      meth <- temp$method
    } else {
      cor.pval[Xs, Ys] <- cor.pval[Ys, Xs] <- NA
      cor.stat[Xs, Ys] <- cor.stat[Ys, Xs] <- NA
    }
  }
  options(warn)
  dimnames(cor.stat) <- list(nams, nams)
  dimnames(cor.pval) <- list(nams, nams)
  retval <- list(estimates = cor.stat, p.values = cor.pval, counts = cor.count,
                 alternative = "two.sided", na.method = na.methods[na.meth],
                 method = meth, data.name = data.name, data=data,
                 call.method=method, distribution=distribution)
  oldClass(retval) <- "cor.all"
  return(retval)
}

