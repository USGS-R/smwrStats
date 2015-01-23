#' Diagnostic Plots
#' 
#' Produce a series of diagnostic plots for statistical analyses.
#' 
#' @details For objects of class "ancovaReg" and "multReg," the argument \code{which}
#' can be a character string, "All" or any of a sequence of numbers from 1 to
#' 8. If it is "All," then all plots are produced. For class "multReg," can
#' also be the name of an explanatory variable so that a residual dependence
#' plot is created for a single variable.\cr Numeric values for \code{which}:
#' \enumerate{ \item Fitted with separate factor levels vs. Observed \item
#' Fitted vs. Residual \item S-L plot \item A correlogram if dates are
#' available in the model or in the data set \item Q-normal and Tukey boxplots
#' for each factor level \item Influence plot \item Outliers plot \item
#' Residual dependence plots for each explanatory variable }
#' 
#' For objects of class "binaryreg," the argument \code{which} can be a
#' character string, "All" or any of a sequence of numbers from 1 to 5. If it
#' is "All," then all plots are produced. By default, the le Cessie-van
#' Houwelingen diagnotic plot is not shown as it can be difficult to interpret.
#' Numeric values for \code{which}: \enumerate{ \item Le Cessie-van Houwelingen
#' overall fit \item Overall fit \item Classification plot \item ROC plot \item
#' Lin-Wei-Ying partial residual plots }
#' 
#' For objects of class "cor.all," the argument \code{which} must be either
#' "All," "Lower," or the name of one of the variables. If \code{which} is
#' "All", then the full scatter plot matrix is shown if there are 4 or fewer
#' variables, otherwise indivdual paired plots are shown. If \code{which} is
#' "Lower", then the lower part of the scatter plot matrix is shown if there
#' are 5 or fewer variables, otherwise indivdual paired plots are shown. If
#' \code{which} is the name of a variable, then only the scatter plots of that
#' variable and all others are shown.\cr
#' 
#' For objects of class "lecessie," the argument \code{which} must be either
#' "All," "First," the name of one or more of the variables, or a number
#' indicating which variable to plot.  If \code{which} is "All", then the
#' fitted vs. Residual and all partial residual plots are created. If
#' \code{which} is "First", then only the fitted vs. Residual plot is created.
#' If \code{which} is one or more of the variable names, then those partial
#' residual plots are created. If \code{which} is numeric, then the fitted vs.
#' Residual (1) sequentially numbered partial residual plots are created.\cr
#' 
#' For objects of class "move.1" or "move.2," the argument \code{which} can be a character
#' string, "All" or any of a sequence of numbers from 1 to 3. If it is "All,"
#' then all plots are produced.  Numeric values for \code{which}: \enumerate{
#' \item x on y and y on x \item The LOC regression line with ellipse \item Q-Q
#' plot of x and y }
#' 
#' For objects of class "roc," the argument \code{which} can be a character
#' string, "All" or 1. The only graph is the receiver operating characteristics
#' curve for a logistic regression.\cr
#' 
#' For objects of class "senSlope," the argument \code{which} can be a
#' character string, "All" or 1. The only graph avaialble is a scatter plot of
#' the data with the regression line in green and a smoothed line in cyan.
#' 
#' @aliases plot.Stats plot.ancovaReg plot.binaryreg plot.cor.all plot.lecessie
#' plot.move.1 plot.move.2 plot.multReg plot.roc plot.senSlope
#' @param x the object to be plotted.
#' @param which a character string or sequence of integers indicating which
#' diagnostic plots to produce. See \bold{Details}.
#' @param set.up set up the graphics page? Set to \code{FALSE} if the graphics
#' page has been set up with a call to \code{setPage}.
#' @param span the smoothing parameter for \code{loess.smooth}.
#' @param \dots not used, required for method function. These diagnostic plots
#' have fixed graphics output.
#' @param bandw the bandwidth for kernel smoothing for the Hosmer-Lemeshow
#' plot.
#' @return The object \code{x} is returned invisibly.
#' @seealso \code{\link{ancovaReg}}, \code{\link{binaryReg}},
#' \code{\link{cor.all}}, \code{\link{leCessie.test}}, \code{\link{move.1}},
#' \code{\link{move.2}}, \code{\link{multReg}}, \code{\link{roc}}, \code{\link{senSlope}},
#' \code{\link{loess.smooth}}, \code{\link{setPage}}
#' @references le Cessie, S. and van Houwelingen, H.C., 1995, Testing the fit
#' of a regression model via score tests in random effects models: Biometrics,
#' v. 51, p. 600--614.\cr Lin, D.Y., Wei, L.J., and Ying, Z., 2002,
#' Model-checking techniques based on cumulative residuals: Biometrics, v. 58,
#' p. 1--12.
#' @keywords hplot
#' @rdname plot.Stats
#' @export
#' @method plot ancovaReg
plot.ancovaReg <- function(x, which='All', set.up=TRUE, span=0.8, ...) {
	## Arguments:
	##  x (ancovareg object) the object to plot
	##  which (character or integer) which graphs to plot
	##  set.up (logical) call setPage?
	##  span (numeric scalar) smoothing parameter for loess.smooth
	## 
	## Identify which plots to do:
	## 1 Fitted - Actual
	## 2 Fitted - Residual
	## 3 S-L; colors by factor.var
	## 4 correlogram (only if dates are available)
	## 5 Q - normal and boxplot by factor.var
	## 6 Influence
	## 7 Outliers
	## 8 Residual dependence plots; colors by factor.var
	##
	## Set up graphics page if requested
	if(set.up)
		setGD("ANCOVA")
	## Set up to do all plots
	doPlot <- rep(TRUE, 8) 
	if(is.numeric(which)) {
		if(min(which) > 0) # select which to plot
			doPlot[seq(8)[-which]] <- FALSE
		else # which not to plot
			doPlot[-which] <- FALSE
	}
	## Anything else produces all plots
	## Final plot (8) residuals vs predictors
	## Extract weights and other data needed later
	Fits <- na.omit(fitted(x$object)) # needed for later plots
	wt.showCD <- if(is.null(x$object$weights))
		rep(1,length(Fits))
	else 
		x$object$weights
	xpred <- x$x.fr # Need the frame to preserve the factor
	if(doPlot[8L]) {
		xnames <- dimnames(xpred)[[2L]]
		xnames <- xnames[xnames != x$factor.var]
		## Residual dependence plots
		if(length(xnames) > 0L) {
			for(i in xnames) {
				## Need to protect against things like fourier and quadratic, which result in matrix entries
				if(is.matrix(xpred[, i])) {
					for(j in seq(ncol(xpred[,i]))) {
						xtopl <- xpred[, i][, j]
						AA <- colorPlot(xtopl, x$diagstats$stud.res,
														color=xpred[[x$factor.var]],
														Plot=list(what='points', size=0.05),
														xtitle=paste(i, colnames(xpred[, i])[j], sep=""), 
														ytitle='Studentized Residual',
														margin=c(NA, NA, 2.4, .5))
						if(span > 0) {
							smo <- loess.smooth(xtopl, x$diagstats$stud.res, family='sym', span=span)
							addXY(smo$x, smo$y)
						}
						refLine(horizontal=0, Plot=list(what='lines', width='standard',
																						type='dashed'))
						## The p-value of the second order fit on the residuals exactly matches
						## the p-value of adding the second order term to the regression
						nl.p <- summary(lm(x$diagstats$resids ~ poly(xtopl, 2),
															 weights=wt.showCD, model=TRUE), FALSE)$coefficients[3L,4L]
						addTitle(Main=paste("Second order polynomial test for linearity: p=",
																round(nl.p,4), sep=''), Bold=FALSE)
						## Must be after Title
						addExplanation(AA, "ul")
					}
				} else {
					AA <- colorPlot(xpred[,i], x$diagstats$stud.res,
													color=xpred[[x$factor.var]],
													Plot=list(what='points', size=0.05),
													xtitle=i, ytitle='Studentized Residual',
													margin=c(NA, NA, 2.4, .5))
					if(span > 0) {
						smo <- loess.smooth(xpred[,i], x$diagstats$stud.res, family='sym', span=span)
						addXY(smo$x, smo$y)
					}
					refLine(horizontal=0, Plot=list(what='lines', width='standard',
																					type='dashed'))
					## The p-value of the second order fit on the residuals exactly matches
					## the p-value of adding the second order term to the regression
					nl.p <- summary(lm(x$diagstats$resids ~ poly(xpred[,i], 2),
														 weights=wt.showCD, model=TRUE), FALSE)$coefficients[3L,4L]
					addTitle(Main=paste("Second order polynomial test for linearity: p=",
															round(nl.p,4), sep=''), Bold=FALSE)
					## Must be after title
					addExplanation(AA, "ul")
				}
			}
		}
	} # end of Residual dependence plots
	## Studentized vs fitted
	if(doPlot[7L]) {
		tval <- qt(.975, x$object$df.residual-1)
		ylim=max(tval, abs(range(x$diagstats$stud.res))) * 1.05
		xyPlot(Fits, x$diagstats$stud.res,
					 Plot=list(what='points', size=0.04),
					 yaxis.range=c(-ylim, ylim),
					 xtitle='Fitted', ytitle='Studentized Residual')
		refLine(horizontal=0, Plot=list(what='lines', width='standard', type='dashed'))
		refLine(horizontal=c(-tval, tval),
						Plot=list(what='lines', width='color', type='dashed', color='red'))
	}
	## Show the effect that each flagged point has on the regression
	## (an influence plot)
	Act <- Fits + residuals(x$object) # needed later
	if(doPlot[6L]) {
		xyPlot(Fits, Act,
					 Plot=list(what='points', size=0.02),
					 xtitle="Fitted", ytitle="Actual")
		refLine(coefficients=c(0,1), Plot=list(what='lines', width='standard', type='dashed'))
		if(is.null(x$object$na.action))
			fg.showCD <- which(x$flagobs)
		else 
			fg.showCD <- which(naresid(x$object$na.action, x$flagobs))
		for(j in seq(along=fg.showCD)) {
			i <- fg.showCD[j]
			## generate a random color
			Col <- sprintf('#%02x%02x%02x', as.integer(runif(1, 10, 245)),
										 as.integer(runif(1, 10, 245)),
										 as.integer(runif(1, 10, 245)))
			refLine(coefficients=lsfit(Fits[-i], Act[-i], wt.showCD[-i])$coef,
							Plot=list(what='lines', color=Col, width='color'))
			addXY(Fits[i], Act[i],
						Plot=list(what='points', color=Col, size=0.04))
			labelPoints(Fits[i], Act[i], dimnames(x$x)[[1L]][i],
									dir='NE', size=12, color=Col)
		}
	} # end of Influence plot
	## Q-normal of standardized residuals (H&H criterion 4)
	if(doPlot[5L]) {
		## First is box plot by factor level
		boxPlot(x$diagstats$stnd.res, group=x$x.fr[[x$factor.var]],
						Box=list(type="tukey"))
		probPlot(x$diagstats$stnd.res,
						 yaxis.log=FALSE, ylabels=7,
						 ytitle='Standardized Residual',
						 margin=c(NA, NA, 2.4, NA), mean=0, sd=1)
		PPCC <- ppcc.test(x$diagstats$stnd.res)$p.value
		addTitle(Main=paste("PPCC test for normality: p=", round(PPCC,4), sep=''), Bold=FALSE)
	}
	## for the next plots use Pearson residuals, which are weighted
	Res <- residuals(x$object, type="pearson")
	## If possible, plot a correlogram--requires finding 1 datelike column in
	## the data
	if(doPlot[4L] && !is.null(x$object$call$data)) {
		data <- x$object$call$data
		if(is.name(data))
			data <- try(eval(data))
		if(class(data) != "try-error") {
			anyDate <- which(sapply(data, isDateLike))
			if(length(anyDate) == 1) { # if more than 1, 
				Date <- dectime(data[[anyDate]])
				if(!is.null(skips <- x$object$na.action)) {
					if(attr(skips, "class") == "omit")
						Date <- Date[-skips] # Must remove if class is omit 
				}
				corGram(Date, Res)
			}
		}
	}
	## Add details of call on regression model to next plots
	Mod <- format(x$object$call$formula)
	## 3rd plot, S-L
	RSE <- rmse(x$object)
	if(doPlot[3L]) {
		xyPlot(Fits, sqrt(abs(Res)),
					 Plot=list(what='points', size=0.05),
					 xtitle="Fitted",
					 ytitle=as.expression(substitute(sqrt(abs(YL)),
					 																list(YL = as.name("Residuals")))),
					 margin=c(NA, NA, 2.4, NA))
		if(span > 0) {
			smo <- loess.smooth(Fits, sqrt(abs(Res)), span=span)
			addXY(smo$x, smo$y)
		}
		refLine(horizontal=0.82218*sqrt(RSE), Plot=list(what='lines', width='standard', type='dashed'))
		Woodings <- cor.test(Fits, abs(Res), method='s', exact=FALSE)
		addTitle(Main=paste("Woodings test for heteroscedasticity: p=",
												round(Woodings$p.value,4), sep=''), Bold=FALSE)
	} # end of S-L 
	## 2nd plot response vs. fit
	if(doPlot[2L]) {
		xyPlot(Fits, Res,
					 Plot=list(what='points', size=0.05),
					 xtitle="Fitted",
					 ytitle="Residuals")
		if(span > 0) {
			smo <- loess.smooth(Fits, Res, span=span)
			addXY(smo$x, smo$y)
		}
		refLine(horizontal=0, Plot=list(what="lines", width="standard", type="dashed"))
	}
	## First plot is actual vs partial fitted, with regression details
	if(doPlot[1L]) {
		## Construct fitted values using only the continuous variables
		Levels <- levels(xpred[[x$factor.var]])
		## This plot requires the input data set
		data <- x$object$call$data
		if(is.null(data)) {
			warning("Plot # 1 requires the original data")
			return(invisible(x))
		}
		else if(is.name(data)) {
			Par.fr <- try(eval(data))
			if(class(data) == "try-error") {
				warning("Plot # 1 requires the original data")
				return(invisible(x))
			}
			## We must protect against finding a symbol with the same name as an object in this function!
			## Why here and not for doPlot4?
			if(is.name(Par.fr)) {
				Par.fr <- as.character(Par.fr)
				Par.fr <- get(Par.fr, pos=1L)
			}
		}
		else # Must be a dataset
			Par.fr <- data
		## This makes a dummy data set for factor = first level
		Par.fr[[x$factor.var]] <- Levels[1L]
		## Make predictions and remove missings to match data
		Par.fits <- predict(x$object, Par.fr)
		if(!is.null(toDel <- x$object$na.action))
			Par.fits <- Par.fits[-toDel]
		AA <- colorPlot(Par.fits, Act, color=xpred[[x$factor.var]],
										Plot=list(what="points", size=0.05),
										xtitle=paste("Partial Fitted: ", x$factor.var, " dropped",
																 sep=""),
										ytitle="Response")
		Colors <- setColor(Levels) # The default for colorPlot
		for(i in seq(along=Levels)) {
			Par <- lm(Act ~ Par.fits, subset=xpred[[x$factor.var]] == Levels[i])
			refLine(coefficients=coefficients(Par),
							Plot=list(what="lines", width="color", color=Colors[i]),
							xrange=range(Par.fits[xpred[[x$factor.var]] == Levels[i]]))
		}
		addExplanation(AA, where="lr")
		## Add some details, regression eqn and RSE
		Eqn <- x$object$coef
		names(Eqn)[1L] <- ""
		Eqn <- paste(as.character(round(Eqn, 3)), names(Eqn), sep=' ')
		Eqn <- paste(Eqn, collapse=' + ')
		Resp <- as.expression(substitute(hat(R), 
																		 list(R=as.name(deparse(x$object$call$formula[[2L]])))))
		Model <- expression()
		RSE <- signif(RSE, 3)
		legend("topleft", legend=c(
			as.expression(substitute(hat(R) == EQN, 
															 list(R=as.name(deparse(x$object$call$formula[[2L]])),
															 		 EQN=Eqn))),
			paste("Residual Standard Error: ", RSE, sep='')), bty='n')
	}
	invisible(x)
}

#' @rdname plot.Stats
#' @export
#' @method plot binaryreg
plot.binaryreg <- function(x, which=2:5, set.up=TRUE, bandw=0.3, ...) {
	## Arguments:
	##  x (binaryreg object) the object to plot
	##  which (character or integer) which graphs to plot
	##  bandw (numeric scalar) bandwidth for kernel smoothing for H-L graph
	##
	## Set up graphics page
	if(set.up) {
		setGD("LOGISTIC")
	}
	## Set up to do all 5 plots
	doPlot <- rep(TRUE,5) 
	if(is.numeric(which)) {
		if(min(which) > 0) # select which to plot
			doPlot[seq(5)[-which]] <- FALSE
		else # which not to plot
			doPlot[-which] <- FALSE
	}
	## Anything else produces all plots
	## 
	## Last are the L-W-Y partial plots
	if(doPlot[5L]) {
		x.resid <- resid(x$object, 'response')
		x.coefs <- coef(x$object)
		x.mat <- model.matrix(x$object)
		for(i in seq(2, length(x.coefs))) {
			xs <- x.mat[,i]
			## Skip if binary
			if(length(unique(xs)) < 3) next
			rtoplot <- x.resid[order(xs)]
			xs <- sort(xs)
			rtoplot <-  cumsum(rtoplot)/sqrt(length(xs))
			xyPlot(xs, rtoplot, Plot=list(what='stairstep'),
						 xtitle=names(x.coefs[i]),
						 ytitle='L-W-Y Cumulative Residuals', yaxis.range=c(-.6, .6))
			refLine(horizontal=0)
		}
	}
	## ROC is plot # 4
	if(doPlot[4L])
		plot(x$roc, set.up=FALSE)
	## Third plot--classification graph
	if(doPlot[3L] && !is.null(x$PctCorrect)) {
		pd.0 <- density(x$object$fitted.values[x$object$y == 0], bw=0.25, kern='tri')
		pd.1 <- density(x$object$fitted.values[x$object$y == 1], bw=0.25, kern='tri')
		ylim <- pretty(c(0, max(pd.0$y, pd.1$y)))
		ylim <- ylim[c(1L, length(ylim))] # pick first and last for range
		AA <- xyPlot(pd.0$x, pd.0$y, Plot=list(name="specificity", 
																					 color='darkblue', what="lines", width="color"),
								 xaxis.range=c(0,1), yaxis.range=ylim,
								 xtitle='Predicted', ytitle='Relative Density', margin=c(NA, NA, 1.6, NA))
		AA <- addXY(pd.1$x, pd.1$y, Plot=list(name="sensitivity", color='darkgreen', width="color"), current=AA)
		addExplanation(AA, where='ur', title='')
		refLine(vertical=0.5)
		addTitle(paste("specificity:", round(x$PctCorrect[["0"]], 3),
									 "  sensitivity:", round(x$PctCorrect[["1"]], 3), sep=' '), Bold=FALSE)
	}
	## Second Plot, overall fit with H-L
	if(doPlot[2L] && !is.null(x$Hosmer)) {
		fits <- fitted(x$object)
		y <- x$object$y
		xyPlot(fits, y, Plot=list(what='points', filled=FALSE),
					 xtitle='Fitted Values', ytitle='Observed values',
					 xlabels=list(labels=5, extend.pct=5, extend.range=FALSE),
					 ylabels=list(labels=5, extend.pct=5, extend.range=FALSE),
					 margin=c(NA, NA, 1.6, NA))
		p2.smo <- ksmooth(fits, y, bandwidth=bandw, kernel='normal')
		addXY(range(fits), range(fits), Plot=list(what='lines'))
		addXY(p2.smo$x, p2.smo$y, Plot=list(what='lines'))
		HLx <- x$Hosmer$estimate[,2] / x$Hosmer$estimate[,1]
		HLy <- x$Hosmer$estimate[,3] / x$Hosmer$estimate[,1]
		addXY(HLx, HLy, Plot=list(what='points', filled=TRUE))
		addTitle(paste("Hosmer-Lemeshow Test, p-value =", round(x$Hosmer$p.value,4)), Bold=FALSE)
	}
	## First plot (last on page), leCessie
	if(doPlot[1L] && !is.null(x$leCessie)) {
		plot(x$leCessie, which=1, set.up=FALSE)
		addTitle(paste("le Cessie-van Houwelingen Test, p-value =", round(x$leCessie$p.value,4)),
						 Bold=FALSE)
	}
	invisible(x)
}

#' @rdname plot.Stats
#' @export
#' @method plot cor.all
plot.cor.all <- function(x, which="All", set.up=TRUE, ...) {
	## Arguments:
	##  x (a cor.all object) the object ot be plotted
	##  which (character scalar) which plots to plot
	##    'All' means full matrix, 'Lower' means lower tri implemented only 
	##    when ncol is 5 or less, other is show correlations for variable
	##  ... (dots) not used, required for method funtion
	##
	## Use splomPlot() with a neet little panel function for ncol(x$data <=4)
	## and paired x-y plots n*(n-1)/2 otherwise. Use a 3.5x6.5 vertical figure layout.
	## Set up graphics page
	if(set.up) 
		setGD("CORRELATIONS")
	## Continue
	data <- x$data
	call.method <- x$call.method
	if(x$distribution == "lognormal") {
		names(data) <- paste("log", names(data), sep='.') # fix for plot only
		rownames(x$estimates) <- paste("log", rownames(x$estimates), sep='.')
		colnames(x$estimates) <- paste("log", colnames(x$estimates), sep='.')
		rownames(x$p.values) <- paste("log", rownames(x$p.values), sep='.')
		colnames(x$p.values) <- paste("log", colnames(x$p.values), sep='.')
	} else if(x$distribution == "log1p") {
		names(data) <- paste("log1p", names(data), sep='.') # fix for plot only
		rownames(x$estimates) <- paste("log1p", rownames(x$estimates), sep='.')
		colnames(x$estimates) <- paste("log1p", colnames(x$estimates), sep='.')
		rownames(x$p.values) <- paste("log1p", rownames(x$p.values), sep='.')
		colnames(x$p.values) <- paste("log1p", colnames(x$p.values), sep='.')
	}
	N <- ncol(data)
	ckN <- as.integer(min(par("fin")) / 1.6249)
	if(which == "All") {
		Nok <- N <= ckN
		show.all <- TRUE
	} else if(which == "Lower") {
		Nok <- N <= (ckN + 1L)
		show.all <- FALSE
	} else { # Show which on the y-axis
		titles <- names(data)
		est <- c(pearson="cor", spearman="rho", kendall="tau")[call.method]
		if(x$distribution == "lognormal") {
			which <- paste("log", which, sep='.')
		} else if(x$distribution == "log1p")
			which <- paste("log1p", which, sep='.')
		for(j in titles) {
			if(j != which) {
				AA.xy <- xyPlot(data[[j]], data[[which]], Plot=list(what="points"),
												ylabels=5, xlabels=5,
												xtitle=j, ytitle=which, margin=c(NA, NA, 1.6, 0.5))
				addTitle(paste(est, " = ", round(x$estimates[which,j], 3),
											 ", p.val = ", round(x$p.values[which,j], 3), sep=''))
				addSmooth(AA.xy$x, AA.xy$y, span=1.0, current=AA.xy)
			}
		}
		return(invisible(x))
	}
	## Proceed with splom/paired plots
	if(Nok) {
		AA <- setSplom(num.variables=N, show.all=show.all, ymargin=3.5)
		splomPlot(data, AA, Panel=function(x, y, cur) {
			test <- cor.test(x, y, method=call.method)
			est <- names(test$est)
			tbl <- rbind(paste(est, " = ", round(test$est, 3), sep=''),
									 paste("p.val = ", round(test$p.value, 3), sep=''))
			addTable(tbl, "ul")
			return(addSmooth(x, y, span=1.0, current=cur))})
	} else { # more than 5 columns
		AA <- setLayout(width=3.5, height=c(3.35, 3.15), , num.rows=2, xtop=1.2)
		needpage <- FALSE
		titles <- names(data)
		est <- c(pearson="cor", spearman="rho", kendall="tau")[call.method]
		for(j in seq(N, 2L, -1L)) {
			for(i in seq(j - 1L, 1L, -1L)) {
				if(needpage)
					plot.new() # Sets up a new page
				needpage <- TRUE
				AA.m <- setGraph(1, AA)
				AA.m[3L] <- 1.6 # Need to allocate space for a title
				AA.xy <- xyPlot(data[[i]], data[[j]], Plot=list(what="points"),
												ylabels=5, xlabels=5,
												xtitle=titles[i], ytitle=titles[j], margin=AA.m)
				addTitle(paste(est, " = ", round(x$estimates[i,j], 3),
											 ", p.val = ", round(x$p.values[i,j], 3), sep=''))
				addSmooth(AA.xy$x, AA.xy$y, span=1.0, current=AA.xy)
				## Now plot j,i
				AA.m <- setGraph(2L, AA)
				AA.xy <- xyPlot(data[[j]], data[[i]], Plot=list(what="points"),
												ylabels=5, xlabels=5,
												xtitle=titles[j], ytitle=titles[i], margin=AA.m)
				addSmooth(AA.xy$x, AA.xy$y, span=1.0, current=AA.xy)
			} # end of j
		} # End of i
	} # end of else (too many to splom)
	invisible(x)
}

#' @rdname plot.Stats
#' @export
#' @method plot lecessie
plot.lecessie <- function(x, which="All", set.up=TRUE, ...) {
	## Arguments:
	##  x (lecessie object) the object to be printed
	##  which (character or numeric) which plots to plot
	##  ... (dots) not used, required for method function
	##
	## Set up graphics page
	if(set.up) 
		setGD("LECESIE")
	target.mat <- model.matrix(x$target.object)
	target.names <- dimnames(target.mat)[[2L]]
	target.mat[, 1L] <- fitted(x$object)
	target.names[1L] <- "Fitted values"
	target.resid <- residuals(x$object, type="dev")
	if(tolower(which[1]) == "all")
		which <- seq(ncol(target.mat), 1L)
	else if(tolower(which[1]) == "first")
		which <- 1
	else if(is.character(which))
		which <- which(target.names %in% which)
	for(i in which) {
		xyPlot(target.mat[,i,drop=TRUE], x$smoothed.residuals,
					 Plot=list(what="points", size=0.05), 
					 ytitle="Smoothed residuals", xtitle=target.names[i],
					 margin=c(NA, NA, 1.6, NA)) # leave room for title
		if(i == 1)
			addXY(target.mat[,i,drop=TRUE], target.resid,
						Plot=list(what="points", filled=FALSE, size=0.08))
	}
	invisible(x)
}

#' @rdname plot.Stats
#' @export
#' @method plot move.1
plot.move.1 <- function(x, which="All", set.up=TRUE, span=0.8, ...) {
	##
	## Identify which plots to do:
	## 1 Q-Q plot
	## 2 x on y, y on x
	## 
	## Set up graphics page
	if(set.up) 
		setGD("MOVE.1")
	fin <- par("fin")
	## Set up to do all plots
	doPlot <- c(rep(TRUE, 2L) , FALSE)
	if(is.numeric(which)) {
		if(min(which) > 0) # select which to plot
			doPlot[seq(2L)[-which]] <- FALSE
		else # which not to plot
			doPlot[-which] <- FALSE
	}
	xname <- x$var.names[2L]
	yname <- x$var.names[1L]
	if(doPlot[3L]) {
		## This is done only by special request
		AA <- scalePlot(x$x, x$y, scale=x$coef[2L], Plot=list(what="points"),
										xtitle=xname, ytitle=yname)
		refLine(coefficients=x$coef, current=AA)
		cov <- x$R*(x$xstats[2L]*x$ystats[2L])
		cov <- matrix(c(x$xstats[2L]^2, cov, cov, x$ystats[2L]^2), 2)
		el <- cov2Ellipse(cov, c(x$xstats[1L], x$ystats[1L]))
		addXY(el$x, el$y, Plot=list(what="lines", color="blue"), current=AA)
		el <- cov2Ellipse(cov, c(x$xstats[1L], x$ystats[1L]), scale=2)
		addXY(el$x, el$y, Plot=list(what="lines", color="blue"), current=AA) 
	}
	if(doPlot[2L]) {
		## Set up for 2 graphs
		AA.lo <- setLayout(width=fin[1L], height=fin[2L]/c(2,2))
		AA.gr <- setGraph(1, AA.lo)
		AA <- xyPlot(x$x, x$y, Plot=list(what="points"), xtitle=xname, ytitle=yname)
		refLine(coefficients=lsfit(x$x, x$y)$coef, Plot=list(color="green"), current=AA)
		addSmooth(x$x, x$y, span=span, Plot=list(color="cyan"), current=AA)
		AA.gr <- setGraph(2, AA.lo)
		AA <- xyPlot(x$y, x$x, Plot=list(what="points"), xtitle=yname, ytitle=xname)
		refLine(coefficients=lsfit(x$y, x$x)$coef, Plot=list(color="green"), current=AA)
		addSmooth(x$y, x$x, span=span, Plot=list(color="cyan"), current=AA)
	}
	if(doPlot[1L]) {
		# Reset if 2 done
		setLayout(width=fin[1L] - .5, height=(fin[2L] - .5))
		qqPlot(x$x, x$y, xtitle=xname, ytitle=yname, Line1.1=list(what="none"))
	}
	invisible(x)
}

#' @rdname plot.Stats
#' @export
#' @method plot move.2
plot.move.2 <- function(x, which="All", set.up=TRUE, span=0.8, ...) {
	##
	## Identify which plots to do:
	## 1 Q-Q plot
	## 2 x on y, y on x
	## 
	## Set up graphics page
	if(set.up) 
		setGD("MOVE.2")
	fin <- par("fin")
	## Set up to do all plots
	doPlot <- c(rep(TRUE, 2L) , FALSE)
	if(is.numeric(which)) {
		if(min(which) > 0) # select which to plot
			doPlot[seq(2L)[-which]] <- FALSE
		else # which not to plot
			doPlot[-which] <- FALSE
	}
	xname <- x$var.names[2L]
	yname <- x$var.names[1L]
	if(doPlot[3L]) {
		## This is done only by special request
		AA <- scalePlot(x$x, x$y, scale=x$coef[2L], Plot=list(what="points"),
										xtitle=xname, ytitle=yname)
		refLine(coefficients=x$coef, current=AA)
		cov <- x$R*(x$xstats[2L]*x$ystats[2L])
		cov <- matrix(c(x$xstats[2L]^2, cov, cov, x$ystats[2L]^2), 2)
		el <- cov2Ellipse(cov, c(x$xstats[1L], x$ystats[1L]))
		addXY(el$x, el$y, Plot=list(what="lines", color="blue"), current=AA)
		el <- cov2Ellipse(cov, c(x$xstats[1L], x$ystats[1L]), scale=2)
		addXY(el$x, el$y, Plot=list(what="lines", color="blue"), current=AA) 
	}
	if(doPlot[2L]) {
		## Set up for 2 graphs
		AA.lo <- setLayout(width=fin[1L], height=fin[2L]/c(2,2))
		AA.gr <- setGraph(1, AA.lo)
		AA <- xyPlot(x$x, x$y, Plot=list(what="points"), xtitle=xname, ytitle=yname)
		refLine(coefficients=lsfit(x$x, x$y)$coef, Plot=list(color="green"), current=AA)
		addSmooth(x$x, x$y, span=span, Plot=list(color="cyan"), current=AA)
		AA.gr <- setGraph(2, AA.lo)
		AA <- xyPlot(x$y, x$x, Plot=list(what="points"), xtitle=yname, ytitle=xname)
		refLine(coefficients=lsfit(x$y, x$x)$coef, Plot=list(color="green"), current=AA)
		addSmooth(x$y, x$x, span=span, Plot=list(color="cyan"), current=AA)
	}
	if(doPlot[1L]) {
		# Reset if 2 done
		setLayout(width=fin[1L] - .5, height=(fin[2L] - .5))
		qqPlot(x$x, x$y, xtitle=xname, ytitle=yname, Line1.1=list(what="none"))
	}
	invisible(x)
}

#' @rdname plot.Stats
#' @export
#' @method plot multReg
plot.multReg <- function(x, which='All', set.up=TRUE, span=1.0, ...) {
	## Arguments:
	##  x (multreg object) the object to plot
	##  which (character or integer) which graphs to plot
	##  span (numeric scalar) smoothing parameter for loess.smooth
	##  set.up and  ... (dots) unused, required for method function
	## 
	## Identify which plots to do:
	## 1 Fitted - Actual
	## 2 Fitted - Residual
	## 3 S-L
	## 4 correlogram (only if dates are available)
	## 5 Q - normal
	## 6 Influence
	## 7 Outliers
	## 8 Residual dependence plots
	##
	## Set up graphics page
	if(set.up)
		setGD("MULTREG")
	## Set up to do all plots
	doPlot <- rep(TRUE, 8L)
	do8 <- FALSE
	if(is.numeric(which)) {
		if(min(which) > 0) # select which to plot
			doPlot[seq(8L)[-which]] <- FALSE
		else # which not to plot
			doPlot[-which] <- FALSE
	}
	else if(is.character(which) && which[1] != "All") {
		doPlot[1:7] <- FALSE
		xnames <- which
		do8 <- TRUE
	}
	## Anything else produces all plots
	## Final plot (8) residuals vs predictors
	## Extract weights and other data needed later
	Fits <- na.omit(fitted(x$object)) # needed for later plots, but without missing values
	wt.showCD <- if(is.null(x$object$weights))
		rep(1,length(Fits))
	else
		x$object$weights
	## Protect aginst intercept only fit
	if(diff(range(Fits)) < 1.e-15) {
		doPlot[-5L] <- FALSE
		doPlot[5L] <- TRUE
	}
	if(doPlot[8L]) {
		xpred <- x$x
		if(!do8) # get all explanatory variable names
			xnames <- dimnames(xpred)[[2L]]
		## Residual dependence plots
		if(do8 || length(xnames) > 1L) {
			for(i in xnames) {
				## Protect against possible binary predictors as created by factors
				if(length(unique(xpred[,i])) < 3L)
					warning(i, " appears to be created by a factor variable, not plotted")
				else {
					xyPlot(xpred[,i], x$diagstats$stud.res,
								 Plot=list(what="points", size=0.05),
								 xtitle=i, ytitle="Studentized Residual",
								 margin=c(NA, NA, 1.5, .5))
					if(span > 0) {
						smo <- loess.smooth(xpred[,i], x$diagstats$stud.res, family="sym", span=span)
						addXY(smo$x, smo$y)
					}
					refLine(horizontal=0, Plot=list(what="lines", width="standard", type="dashed"))
					## The p-value of the second order fit on the residuals exactly matches
					## the p-value of adding the second order term to the regression
					nl.p <- summary(lm(x$diagstats$resids ~ poly(xpred[,i], 2),
														 weights=wt.showCD, model=TRUE), FALSE)$coefficients[3,4]
					addTitle(Main=paste("Second order polynomial test for linearity: p=",
															round(nl.p, 4L), sep=""), Bold=FALSE)
				}
			}
		}
	} # end of Residual dependence plots
	## Studentized vs fitted
	if(doPlot[7L]) {
		tval <- qt(.975, x$object$df.residual - 1)
		ylim=max(tval, abs(range(x$diagstats$stud.res))) * 1.05
		xyPlot(Fits, x$diagstats$stud.res,
					 Plot=list(what="points", size=0.04),
					 yaxis.range=c(-ylim, ylim),
					 xtitle="Fitted", ytitle="Studentized Residual")
		refLine(horizontal=0, Plot=list(what="lines", width="standard", type="dashed"))
		refLine(horizontal=c(-tval, tval),
						Plot=list(what="lines", width="color", type="dashed", color="red"))
	}
	## Show the effect that each flagged point has on the regression
	## (an influence plot)
	Act <- Fits + na.omit(residuals(x$object)) # needed later
	if(doPlot[6L]) {
		xyPlot(Fits, Act,
					 Plot=list(what="points", size=0.02),
					 xtitle="Fitted", ytitle="Actual")
		refLine(coefficients=c(0,1), Plot=list(what="lines", width="standard", type="dashed"))
		if(is.null(x$object$na.action))
			fg.showCD <- which(x$flagobs)
		else 
			fg.showCD <- which(naresid(x$object$na.action, x$flagobs))
		for(j in seq(along=fg.showCD)) {
			i <- fg.showCD[j]
			## generate a random color
			Col <- sprintf("#%02x%02x%02x", as.integer(runif(1, 10, 245)),
										 as.integer(runif(1, 10, 245)),
										 as.integer(runif(1, 10, 245)))
			refLine(coefficients=lsfit(Fits[-i], Act[-i], wt.showCD[-i])$coef,
							Plot=list(what="lines", color=Col, width="color"))
			addXY(Fits[i], Act[i],
						Plot=list(what="points", color=Col, size=0.04))
			labelPoints(Fits[i], Act[i], dimnames(x$x)[[1L]][i],
									dir="NE", size=12, color=Col)
		}
	} # end of Influence plot
	## Q-normal of standardized residuals (H&H criterion 4)
	if(doPlot[5L]) {
		qqPlot(x$diagstats$stnd.res, Plot=list(size=0.05),
					 yaxis.log=FALSE, ylabels=7,
					 ytitle="Standardized Residual",
					 xtitle="Standard Normal Quantiles",
					 margin=c(NA, NA, 2.4, NA))
		refLine(coefficients=c(0,1))
		PPCC <- ppcc.test(x$diagstats$stnd.res)$p.value
		addTitle(Main=paste("PPCC test for normality: p=", round(PPCC,4), sep=""), Bold=FALSE)
	}
	## for the next plots use Pearson residuals, which are weighted
	Res <- na.omit(residuals(x$object, type="pearson"))
	## If possible, plot a correlogram--requires finding 1 datelike column in
	## the data
	if(doPlot[4L] && !is.null(x$object$call$data)) {
		data <- x$object$call$data
		if(is.name(data))
			data <- try(eval(data))
		if(class(data) != "try-error") {
			anyDate <- which(sapply(data, isDateLike))
			if(length(anyDate) == 1L) { # if more than 1, skip it
				Date <- dectime(data[rownames(data) %in% rownames(xpred),anyDate])
				corGram(Date, Res)
			}
		}
	}
	## Add details of call on regression model to next plots
	Mod <- format(x$object$call$formula)
	## 3rd plot, S-L
	RSE <- rmse(x$object)
	if(doPlot[3L]) {
		xyPlot(Fits, sqrt(abs(Res)),
					 Plot=list(what="points", size=0.05),
					 xtitle="Fitted",
					 ytitle=as.expression(substitute(sqrt(abs(YL)),
					 																list(YL = as.name("Residuals")))),
					 margin=c(NA, NA, 2.4, NA))
		if(span > 0) {
			smo <- loess.smooth(Fits, sqrt(abs(Res)), span=span)
			addXY(smo$x, smo$y)
		}
		## 0.82218 is the expected value of the sqrt(abs(x)) for a normal dist.:
		## integrate(function(x) sqrt(abs(x))*dnorm(x), -Inf, Inf)
		refLine(horizontal=0.82218*sqrt(RSE), Plot=list(what="lines", width="standard", type="dashed"))
		Woodings <- cor.test(Fits, abs(Res), method="s", exact=FALSE)
		addTitle(Main=paste("Woodings test for heteroscedasticity: p=",
												round(Woodings$p.value,4), sep=""), Bold=FALSE)
	} # end of S-L 
	## 2nd plot response vs. fit
	if(doPlot[2L]) {
		xyPlot(Fits, Res,
					 Plot=list(what="points", size=0.05),
					 xtitle="Fitted",
					 ytitle="Residuals")
		if(span > 0) {
			smo <- loess.smooth(Fits, Res, span=span)
			addXY(smo$x, smo$y)
		}
		refLine(horizontal=0, Plot=list(what="lines", width="standard", type="dashed"))
	}
	## First plot is actual vs fitted, with regression details
	if(doPlot[1L]) {
		xyPlot(Fits, Act,
					 Plot=list(what="points", size=0.05),
					 xtitle=paste("Fitted:", Mod, sep=" "),
					 ytitle="Response")
		if(span > 0) {
			smo <- loess.smooth(Fits, Act, span=span)
			addXY(smo$x, smo$y)
		}
		refLine(coefficients=c(0,1), Plot=list(what="lines", width="standard", type="dashed"))
		## Add some details, regression eqn and RSE
		Eqn <- x$object$coef
		names(Eqn)[1L] <- ""
		Eqn <- paste(as.character(round(Eqn, 3)), names(Eqn), sep=" ")
		Eqn <- paste(Eqn, collapse=" + ")
		Resp <- as.expression(substitute(hat(R), 
																		 list(R=as.name(deparse(x$object$call$formula[[2]])))))
		Model <- expression()
		RSE <- signif(RSE, 3)
		legend("topleft", legend=c(
			as.expression(substitute(hat(R) == EQN, 
															 list(R=as.name(deparse(x$object$call$formula[[2]])),
															 		 EQN=Eqn))),
			paste("Residual Standard Error: ", RSE, sep="")), bty="n")
	}
	invisible(x)
}

#' @rdname plot.Stats
#' @export
#' @method plot senSlope
plot.senSlope <- function(x, which="All", set.up=TRUE, span=0.8, ...) {
	## Coding history:
	##    2013Apr15 DLLorenz Original Coding
	##
	##
	## Identify which plots to do:
	## y on x with regression line
	##
	## Set up graphics page
	if(set.up) 
		setGD("SenSlope")
	## Set up to do all plots
	doPlot <- TRUE    
	if(is.numeric(which)) {    
		if(which[1L] == -1) # why is beyond me!
			doPlot <- FALSE    
		if(length(which) > 1 || which != 1)
			warning("Only one diagnostic plot for senSlope")
	}
	xname <- x$var.names[2L]
	yname <- x$var.names[1L]
	if(doPlot[1L]) {
		AA <- xyPlot(x$x, x$y, Plot=list(what="points"), xtitle=xname, ytitle=yname)
		refLine(coefficients=x$coefficients, Plot=list(color="green"), current=AA)
		addSmooth(x$x, x$y, span=span, Plot=list(color="cyan"), current=AA)
	}
	invisible(x)
}

#' @rdname plot.Stats
#' @export
#' @method plot roc
plot.roc <- function(x, which="All", set.up=TRUE, ...) {
	## Arguments:
	##  x (roc object) the object to plot
	##   ... (dots) unused, required for method function
	## 
	## Set up graphics page
	if(set.up) 
		setGD("ROC")
	## Set up to do all plots
	doPlot <- TRUE    
	if(is.numeric(which)) {    
		if(which[1L] == -1) # why is beyond me!
			doPlot <- FALSE    
		if(length(which) > 1 || which != 1)
			warning("Only one diagnostic plot for senSlope")
	}
	if(doPlot[1L]) {
		sens <- x$table$sens
		spec<- 1-x$table$spec
		fit <- x$table$fits
		## Need to start at origin
		xyPlot(c(1,spec,0), c(1,sens,0),
					 xaxis.range=c(0,1), yaxis.range=c(0,1),
					 ytitle='True positive rate (Sensitivity)',
					 xtitle='False positive rate (1-Specificity)',
					 margin=c(NA,NA, 1.6, NA))
		refLine(coefficients=c(0,1))
		## Select labels along the length of the curve and plot the predicted value
		curve.dist <- cumsum(diff(c(0, 1 - spec))^2 + diff(c(0, 1 - sens))^2)
		for(i in seq(.05, .95, by=0.05) * max(curve.dist)) {
			xdist <- abs(curve.dist - i)
			min.xdist <-  min(xdist)
			pick <- which(xdist == min.xdist)[1]
			xpick <- spec[pick]
			ypick <- sens[pick]
			if(fit[pick] > 0.008 && fit[pick] < 0.993) { # label it
				if(xpick > .04 && ypick < .96) {
					addAnnotation(xpick-0.01, ypick-0.01,
												format(round(fit[pick], 2)), angle=-45, justification='right')
				}
				else { # Put on the left side of the line
					addAnnotation(xpick+0.01, ypick-0.03,
												format(round(fit[pick], 2)), angle=-45, justification='left')
				}
			}
		}
		addTitle(Main=paste('ROC Analysis, Area under curve =',
												round(x$c.val,3), sep=' '))
	}
	invisible(x)
}

