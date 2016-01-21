#' General tools for hydrologic data and trend analysis.
#' 
#' Useful tools for the analysis of hydrologic data, not inlcuding censored
#' water-quality data.
#' 
#' \tabular{ll}{ 
#'Package: \tab smwrStats\cr 
#'Type: \tab Package\cr 
#'Version: \tab 0.7.5\cr 
#'Date: \tab 2016-01-21\cr 
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
#' @author Dave Lorenz <lorenz@@usgs.gov>\cr
#' 
#' Maintainer: Dave Lorenz <lorenz@@usgs.gov>
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
	packageStartupMessage("Although this software program has been used by the U.S. Geological Survey (USGS), no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.")
}


