#' General tools for hydrologic data and trend analysis.
#' 
#' Useful tools for the analysis of hydrologic data, not inlcuding censored
#' water-quality data.
#' 
#' \tabular{ll}{ 
#'Package: \tab smwrStats\cr 
#'Type: \tab Package\cr 
#'License: \tab CC0\cr 
#'}
#' Regression applications:\cr allReg\cr binaryReg\cr hosmerLemeshow.test\cr
#' leCessie.test\cr multReg\cr press\cr rmse\cr selBestWave\cr
#' roc\cr vif\cr predictDuan\cr predictFerguson\cr predictMVUE\cr
#' senSlope\cr seasonalPeak\cr seasonalWave\cr
#' 
#' Record extension applications:\cr move.1\cr jackknifeMove.1\cr move.2\cr
#' optimBoxCox\cr
#' 
#' Trend applications:\cr curvi\cr kensen.test\cr regken\cr seaken\cr 
#' trends\cr serial.test\cr regken\cr trends\cr
#' 
#' Summary statistics:\cr cor.all\cr printCor\cr percentile\cr
#' quantile.numeric\cr qtiles.CI\cr skew\cr sumStats\cr timeWeightedMean\cr
#' multicomp.test\cr 
#' 
#' Hypothesis tests:\cr grubbs.test\cr ppcc.test\cr serial.test\cr
#' 
#' @name smwrStats-package
#' @aliases smwrStats-package smwrStats
#' @docType package
#' @author Dave Lorenz \cr
#' 
#' @references Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in
#' water resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' 
#' Lorenz, D.L., in preparation. smwrStats---an R package for the analysis of
#' hydrologic data, version 0.7.5.
#' @keywords package
#' @import smwrBase
#' @import smwrGraphs
NULL
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This information is preliminary or provisional and
is subject to revision. It is being provided to meet
the need for timely best science. The information
has not received final approval by the U.S. Geological
Survey (USGS) and is provided on the condition that
neither the USGS nor the U.S. Government shall be held
liable for any damages resulting from the authorized
or unauthorized use of the information.

****Orphaned Package****
This package is looking for a new maintainer. For more information, 
see: https://owi.usgs.gov/R/packages.html#orphan")
}


