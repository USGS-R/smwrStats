#' Select the "Best" Seasonal Wave
#' 
#' Selects the "best" parameters for a seasonal wave fit.
#' 
#' For logistic regression, use \code{Regression}=\code{glm},
#' \code{Test}=\code{deviance}, \code{family=binomial(link="logit")}.
#' 
#' @param formula the formula describing the model without a
#' \code{seasonalWave} term.
#' @param data the data.frame that contians the veriable specified in the
#' formula.
#' @param dec.time a character string of the name of the column in \code{data}
#' in decimal time format.
#' @param wave.list an object of class "seasonalPeak" (confirmed) describing
#' the timing of the peak and potential candidate models.
#' @param exhaustive logical; if TRUE, then do a fairly complete search for the
#' timing of the peak, otherwise accept the timing specified in
#' \code{wave.list}.
#' @param Regression the regression function.
#' @param Test the function to perform the comparison test among all of the
#' candidate models.
#' @param \dots any additional arguments to \code{Regression}.
#' @return A 4- or 7-column matrix of the "best" models. The columns are the
#' timing of the peak (\code{Cmax}), the primary peak loading (\code{Loading}),
#' the half life (\code{Hlife}), and the test score (\code{Test}). If the model
#' is a two-peak model, then additional columns (\code{la}, \code{lo}, and
#' \code{w}) describeing the second peak are included.
#' @references Vecchia, A.V., Martin, J.D., and Gilliom, R.J., 2008, Modeling
#' variability and trends in pesticide concentrations in streams: Journal of
#' the American Water Resources Association, v. 44, no. 5, p.1308--1324.
#' @keywords models regression
#' @examples
#' 
#' ## See the SeasonalWave demo
#' 
#' @export selBestWave
selBestWave <- function(formula, data, dec.time, wave.list, exhaustive=FALSE,
                        Regression=lm, Test=AIC, ...) {
	# Coding history:
	#    2007Aug28 DLLorenz Initial coding.
	#    2007Aug29 DLLorenz Code tweaks and revision
	#    2007Sep14 DLLorenz Modified to use any regression technique
	#    2011Aug11 DLLorenz Conversion to R
	#    2014Dec29 DLLorenz Conversion to roxygen header
  ## 
  ## get candidate models, hlives, and cmax (WaveNum)
  NumPeaks <- attr(wave.list, "NumPeaks")
  if(is.null(NumPeaks)) {
    NumPeaks <- 1L
    warning("Number of peaks not set in wave.list, set to 1 peak")
  }
  Models <- attr(wave.list, "Models")
  if(is.null(Models)) {
    Models <- 1:6
    warning("Models not set in wave.list, set loading periods 1 to 6 months")
  }
  Models <- data.frame(Loading=Models)
  Hlives <- attr(wave.list, "hlife")
  if(is.null(Hlives)) {
    Hlives <- 1:4
    warning("Half lives not set in wave.list, set to 1 to 4 months")
  }
  Seconds <- attr(wave.list, "Second.peak")
  if(!is.null(Seconds))
    Models <- merge(Models, Seconds, all=TRUE)
  WaveNum <- round(as.double(wave.list), 4)
  ## Convert formula to character string to facilitate updates
  formula.str <- as.character(formula)
  AICmat <- matrix(0, nrow=nrow(Models), ncol=length(Hlives))
  dimnames(AICmat) <- list(NULL, as.character(Hlives))
  sp <- NULL
  for(i in seq(nrow(Models))) {
    imod <- Models$Loading[i]
    if(NumPeaks == 2L) {
      sp <- Models[i, -1]
      sp <- paste("list(la=", sp$la, ",lo=", sp$lo, ",w=", sp$w, ")", sep="")
    }
    for(j in seq(along=Hlives)) { # one for each half-life
      jlif <- Hlives[j]
      WaveTerm <- paste("seasonalWave(", dec.time, ",", WaveNum, ",", imod,
                        ",", jlif, ",", sp, ")", sep="")
      formula.model <- paste(formula.str[2], '~', formula.str[3], " + ",
                             WaveTerm, sep='')
      modelOut <- Regression(as.formula(formula.model), data=data,
                           na.action=na.exclude, ...)
      AICmat[i,j] <- Test(modelOut)
    }
  }
  if(exhaustive) { # Do a not very efficient search
    Sel <- whichRowCol(AICmat == min(AICmat))
    Model <- Models$Loading[Sel[1,1]]
    Hlife <- Hlives[Sel[1,2]]
    if(NumPeaks == 2L) {
      sp <- Models[Sel[1,1], -1]
      sp <- paste("list(la=", sp$la, ",lo=", sp$lo, ",w=", sp$w, ")", sep="")
    }
    else
      sp <- NULL
    ## Coarse grid search for best peak within .2 for broad peak
    ##  and .1 for narrow peak (2 peaks)
    if(NumPeaks == 2)
      Step=0.03
    else
      Step <- 0.005*pmax(5, Model + Hlife)
    AICck <- double(9)
    for(k in seq(-4, 4)) {
      Time <- WaveNum + k * Step
      if(Time < 0)
        Time <- Time + 1
      else if(Time > 1)
        Time <- Time - 1
      WaveTerm <- paste("seasonalWave(", dec.time, ",", Time,
                        ",", Model, ",", Hlife, ",", sp, ")", sep="")
      formula.model <- paste(formula.str[2], '~', formula.str[3], " + ",
                             WaveTerm, sep='')
      modelOut <- Regression(as.formula(formula.model), data=data,
                             na.action=na.exclude, ...)
      AICck[k+5] <- Test(modelOut)
    }
    ## Fine grid search with expanded model and half life included
    k <- which(AICck == min(AICck)) - 5
    WaveNum <- WaveNum + k * Step
    Step=Step/2.5
    ## Create the fine time steps
    WaveNum <- WaveNum + Step * seq(-3,3)
    WaveNum[WaveNum < 0] <- WaveNum[WaveNum < 0] + 1
    WaveNum[WaveNum > 1] <- WaveNum[WaveNum > 1] - 1
    Mods <- seq(max(Model - 1, min(Models$Loading)),
                min(Model + 1, max(Models$Loading)))
    Models <- Models[Models$Loading %in% Mods, , drop=FALSE]
    Hlives <- seq(max(Hlife - 1, 1), min(Hlife + 1, 4))
    AICmat <- matrix(0, nrow=nrow(Models)*7, ncol=length(Hlives))
    for(i in seq(nrow(Models))) {
      imod <- Models$Loading[i]
      if(NumPeaks == 2L) {
        sp <- Models[i, -1]
        sp <- paste("list(la=", sp$la, ",lo=", sp$lo, ",w=", sp$w, ")", sep="")
      }
      for(j in seq(along=Hlives)) { # one for each half-life
        jlif <- Hlives[j]
        for(k in seq(7)) { # the fine time grid
          WaveTerm <- paste("seasonalWave(", dec.time, ",", WaveNum[k], ",",
                            imod, ",", jlif, ",", sp, ")", sep="")
          formula.model <- paste(formula.str[2], '~', formula.str[3], " + ",
                                 WaveTerm, sep='')
          modelOut <- Regression(as.formula(formula.model), data=data,
                                 na.action=na.exclude, ...)
          AICmat[(i-1)*7 + k,j] <- Test(modelOut)
        }
      }
    }
    Best <- whichRowCol(AICmat < min(AICmat) + 2)
    BstM <- (Best[,1] - 1) %/% 7 + 1
    BstC <- (Best[,1] - 1) %% 7 + 1
    retval <- cbind(Cmax=WaveNum[BstC], Models[BstM,, drop=FALSE],
                    Hlife=Hlives[Best[,2]], Test=AICmat[Best])
  }
  else { # not exhaustive
    Best <- whichRowCol(AICmat < min(AICmat) + 2)
    retval <- cbind(Cmax=rep(WaveNum, nrow(Best)),
                    Models[Best[,1],, drop=FALSE],
                    Hlife=Hlives[Best[,2]], Test=AICmat[Best])
  }
  retval <- retval[order(retval[,"Test"]),]
  rownames(retval) <- as.character(seq(nrow(retval))) # Make pretty & sensible
  return(retval)
}
